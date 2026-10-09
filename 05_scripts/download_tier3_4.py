#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tier3–4 GEO 批量自动下载器（零依赖，仅标准库）

功能：
  * 解析 crc_gse_download_list_tier3_4.txt，抽取所有真实 GSE 编号（跳过注释/已用/脚本片段）。
  * 按 GEO 真实 FTP 目录规则自动定位 matrix 文件并下载到 <out-root>/<GSE>/。
  * 断点续传、失败重试、并发（默认 3，礼貌 NCBI）、下载清单 + 失败列表。
  * RNA-seq / scRNA-seq 系列（行注释含 RNA/scRNA）额外尝试拉 supplementary RAW 包（非致命）。
  * --dry-run 只解析并打印计划 URL，不联网（用于在本机先验证）。
  * 实时进度显示：解析阶段计数、下载阶段进度条 + 每文件状态 + 大文件内进度（标准库实现）。

用法（在"可访问 GEO"的环境里运行，Windows 下用 D:/AIwork/GEO 这类路径）：
  python download_tier3_4.py                      # 默认读同目录清单，下到 D:/AIwork/GEO
  python download_tier3_4.py --out-root D:/AIwork/GEO --workers 3
  python download_tier3_4.py --only-tier 3        # 只下 Tier3
  python download_tier3_4.py --dry-run            # 仅打印计划
  python download_tier3_4.py --list mylist.txt --exclude GSE999 GSE888

随后转换 + 检测（复用已修复流水线）：
  python -m geo_matrix_to_sop --root D:/AIwork/GEO --out D:/AIwork/GEO_sop
  python analyze_tier.py --indirs D:/AIwork/GEO_sop [--indirs E:/CRC_GEO_2026_combined_sop] \
      --out D:/AIwork/GEO_dedup --name "Tier3+4" --pos1 GSE41258 --pos2 GSE68468
"""
import argparse
import json
import os
import re
import sys
import time
import queue
import threading
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed

# ----------------------------- 配置 -----------------------------
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/"
DEFAULT_LIST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "crc_gse_download_list_tier3_4.txt")
EXCLUDE_DEFAULT = {  # Tier1+Tier2 已用，避免重复下载
    "GSE110224","GSE115261","GSE14297","GSE17536","GSE17538","GSE18392",
    "GSE20916","GSE33113","GSE35982","GSE37364","GSE37892","GSE38832",
    "GSE39582","GSE41258","GSE41568","GSE41655","GSE41657","GSE44076",
    "GSE49355","GSE68468",
}
USER_AGENT = "geo-redundancy-sop-downloader/1.0 (academic; contact: user)"
TIMEOUT = 60          # 文件下载超时（秒）
RESOLVE_TIMEOUT = 20  # 目录桶解析(HEAD)超时（秒，失败更快暴露）
CHUNK = 1 << 16       # 64KB
MSG_Q = queue.Queue()  # 下载阶段的日志队列（由状态线程统一刷屏，避免并发错行）


# ------------------------- 输出辅助 -------------------------
def eprint(msg):
    """在非状态刷新阶段直接打印（带 flush，避免缓冲吞输出）。"""
    print(msg, flush=True)


def qlog(msg):
    """下载阶段：把消息投递给日志队列，由状态线程统一打印。"""
    MSG_Q.put(msg)


def _bar(done, total, width=24):
    if total <= 0:
        return "·" * width
    filled = int(width * done / total)
    return "█" * filled + "░" * (width - filled)


# ------------------------- GEO 目录规则 -------------------------
def bucket_candidates(gse):
    """返回可能的 GEO 目录桶（按真实 FTP 规则 + 容错）。
    真实规则（已对 4/5/6 位 accession 实测验证）：把 GSE 数字的最后 3 位
    替换成 'nnn'，即 bucket = 'GSE' + digits[:-3] + 'nnn'。
      4 位: GSE4107  -> GSE4nnn
      5 位: GSE32323 -> GSE32nnn ;  GSE24514 -> GSE24nnn
      6 位: GSE156451-> GSE156nnn ; GSE110223-> GSE110nnn
    """
    digits = gse[3:]
    n = len(digits)
    primary = ("GSE" + digits[:-3] + "nnn") if n > 3 else ("GSE" + digits[0] + "nnn")
    cands = [primary]
    # 容错：额外尝试"保留前 2 位"和"保留前 3 位" + nnn
    if n >= 2:
        cands.append("GSE" + digits[:2] + "nnn")
    if n >= 3:
        cands.append("GSE" + digits[:3] + "nnn")
    # 去重保序
    seen, out = set(), []
    for c in cands:
        if c not in seen:
            seen.add(c); out.append(c)
    return out


def matrix_url(gse, bucket):
    return f"{BASE}{bucket}/{gse}/matrix/{gse}_series_matrix.txt.gz"


def supp_url(gse, bucket):
    # GEO supplementary 常见两种：整体 tar 或分目录；先试 tar
    return f"{BASE}{bucket}/{gse}/suppl/{gse}_RAW.tar"


# ------------------------- 清单解析 -------------------------
def parse_list(path, only_tier=None, exclude=None):
    exclude = set(exclude or [])
    entries = []  # list of dict: {gse, tier, rna}
    tier = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            # 段落标题（Tier 3 / Tier 4）
            m = re.match(r"\s*#\s*=*\s*Tier\s*([34])\b", line)
            if m:
                tier = int(m.group(1)); continue
            # 真实条目行：兼容两种写法
            #   "GSExxxx   # 注释"   （原格式）
            #   "GSExxxx"           （裸编号，如 gen_tier4_list.py 生成的清单）
            m = re.match(r"\s*(GSE\d+)\s*(?:#(.*))?\s*$", line)
            if not m:
                continue
            gse = m.group(1)
            if gse in exclude:
                continue
            comment = (m.group(2) or "").lower()
            rna = ("rna" in comment) or ("scrna" in comment)
            if only_tier and tier != only_tier:
                continue
            entries.append({"gse": gse, "tier": tier, "rna": rna})
    # 去重（同 GSE 只保留一次，优先标记为 rna 的）
    by = {}
    for e in entries:
        if e["gse"] not in by or (e["rna"] and not by[e["gse"]]["rna"]):
            by[e["gse"]] = e
    return list(by.values())


# ------------------------- 网络工具 -------------------------
def _req(url, method="GET", range_start=None):
    headers = {"User-Agent": USER_AGENT}
    if range_start:
        headers["Range"] = f"bytes={range_start}-"
    return urllib.request.Request(url, headers=headers, method=method)


def resolve_bucket(gse):
    for b in bucket_candidates(gse):
        url = matrix_url(gse, b)
        try:
            with urllib.request.urlopen(_req(url, "HEAD"), timeout=RESOLVE_TIMEOUT) as r:
                if r.status in (200, 206):
                    return b, url
        except urllib.error.HTTPError as e:
            if e.code == 404:
                continue
        except Exception:
            continue
    return None, None


def _stream(url, dest, cb=None):
    """流式下载到 dest，支持断点续传。cb(downloaded_bytes, final=False) 用于进度回调。"""
    start = 0
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        start = os.path.getsize(dest)
    req = _req(url, "GET", range_start=start if start else None)
    downloaded = start
    last_log = start
    tries = 3
    while tries > 0:
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
                mode = "ab" if start and r.status == 206 else "wb"
                with open(dest, mode) as f:
                    while True:
                        buf = r.read(CHUNK)
                        if not buf:
                            break
                        f.write(buf)
                        downloaded += len(buf)
                        if cb and downloaded - last_log >= 2 * 1024 * 1024:
                            cb(downloaded)
                            last_log = downloaded
                if cb:
                    cb(downloaded, final=True)
                return
        except urllib.error.HTTPError as e:
            if e.code == 416:  # 范围不满足，重新完整下载
                start = 0
                req = _req(url, "GET")
                downloaded = 0
                last_log = 0
            else:
                tries -= 1
                time.sleep(2)
        except Exception:
            tries -= 1
            time.sleep(2)
    raise RuntimeError("重试后仍失败")


def download_one(gse, rna, out_root, no_rna_supp=False, cb=None):
    bucket, url = resolve_bucket(gse)
    res = {"gse": gse, "rna": rna, "ok": False, "matrix": None,
           "supp": None, "error": "", "size": 0}
    if not url:
        res["error"] = "无法解析 GEO 目录桶(全部 404)"
        return res
    dest_dir = os.path.join(out_root, gse)
    os.makedirs(dest_dir, exist_ok=True)
    mpath = os.path.join(dest_dir, gse + "_series_matrix.txt.gz")
    # ---- 下载 matrix（带续传） ----
    try:
        _stream(url, mpath, cb=cb)
        res["matrix"] = mpath
        res["size"] = os.path.getsize(mpath) if os.path.exists(mpath) else 0
        res["ok"] = True
    except Exception as e:
        res["error"] = f"matrix 下载失败: {e}"
        return res
    # ---- RNA-seq：额外拉 supplementary（非致命） ----
    if rna and not no_rna_supp:
        surl = supp_url(gse, bucket)
        try:
            with urllib.request.urlopen(_req(surl, "HEAD"), timeout=RESOLVE_TIMEOUT) as r:
                if r.status in (200, 206):
                    spath = os.path.join(dest_dir, gse + "_RAW.tar")
                    _stream(surl, spath)
                    res["supp"] = spath
        except Exception as e:
            res["supp"] = f"(跳过) {e}"
    return res


# ------------------------- 状态刷新线程 -------------------------
def status_loop(st, stop):
    """统一刷新进度条并打印队列中的消息，避免并发写屏错行。"""
    while True:
        # 1) 打印队列中等待的消息（清掉状态行再换行，避免覆盖）
        while True:
            try:
                item = MSG_Q.get_nowait()
            except queue.Empty:
                break
            sys.stdout.write("\r\033[K" + item + "\n")
            sys.stdout.flush()
        # 2) 刷新进度条
        with st["lock"]:
            done = st["done"]; total = st["total"]
            mb = st["bytes"] / (1024 * 1024); act = len(st["active"])
        pct = (done * 100 // total) if total else 0
        line = (f"进度 [{_bar(done, total)}] {done}/{total} ({pct}%)"
                f"  活跃 {act}  已下 {mb:.1f} MB")
        sys.stdout.write("\r\033[K" + line)
        sys.stdout.flush()
        # 3) 结束条件
        if stop.is_set():
            # 再清空一次队列后退出
            while True:
                try:
                    item = MSG_Q.get_nowait()
                except queue.Empty:
                    break
                sys.stdout.write("\r\033[K" + item + "\n")
                sys.stdout.flush()
            sys.stdout.write("\r\033[K")
            sys.stdout.flush()
            break
        time.sleep(0.4)


# ------------------------- 主流程 -------------------------
def main():
    ap = argparse.ArgumentParser(description="Tier3–4 GEO 批量下载器（带进度显示）")
    ap.add_argument("--list", default=DEFAULT_LIST, help="下载清单 .txt")
    ap.add_argument("--out-root", default="D:/AIwork/GEO", help="下载根目录")
    ap.add_argument("--workers", type=int, default=3, help="并发数(建议≤3, 礼貌 NCBI)")
    ap.add_argument("--only-tier", type=int, choices=[3, 4], help="只下指定 Tier")
    ap.add_argument("--exclude", nargs="*", default=(), help="额外排除的 GSE 号")
    ap.add_argument("--dry-run", action="store_true", help="只打印计划，不下载")
    ap.add_argument("--no-rna-supp", action="store_true", help="不拉 RNA-seq 的 supplementary")
    args = ap.parse_args()

    exclude = set(EXCLUDE_DEFAULT) | set(args.exclude)
    entries = parse_list(args.list, only_tier=args.only_tier, exclude=exclude)
    if not entries:
        sys.exit("[错误] 清单未解析到任何 GSE 条目")
    eprint(f"解析到 {len(entries)} 个待下载 GSE"
           + (f" (Tier{args.only_tier} 仅)" if args.only_tier else "")
           + f"，输出根: {args.out_root}")

    # 预解析 bucket
    plan = []
    if args.dry_run:
        eprint("---- 计划（仅本地推算，未联网）----")
        for i, e in enumerate(entries, 1):
            b = bucket_candidates(e["gse"])[0]
            url = matrix_url(e["gse"], b)
            tag = "[RNA]" if e["rna"] else "      "
            plan.append((e, b, url))
            eprint(f"  [{i:>3}/{len(entries)}] {e['gse']:12s} tier={e['tier']} {tag} -> {url}")
        eprint(f"[dry-run] 计划下载 {len(plan)} 个，未实际联网。去掉 --dry-run 执行真实下载。")
        return

    # 真实下载：先解析目录桶（串行、带计数，避免“静默卡住”观感）
    eprint("---- 解析 GEO 目录桶（定位 matrix 文件）----")
    for i, e in enumerate(entries, 1):
        b, url = resolve_bucket(e["gse"])
        if url:
            eprint(f"  [{i:>3}/{len(entries)}] {e['gse']:12s} 桶={b} OK")
        else:
            eprint(f"  [{i:>3}/{len(entries)}] {e['gse']:12s} 桶解析失败(全部 404) -> 将标记 FAIL")
        plan.append((e, b, url))

    os.makedirs(args.out_root, exist_ok=True)
    eprint(f"---- 开始下载 {len(plan)} 个 GSE -> {args.out_root}  (并发={args.workers}) ----")

    # 共享进度状态
    st = {"lock": threading.Lock(), "done": 0, "total": len(plan),
          "bytes": 0, "active": set()}
    stop = threading.Event()
    t = threading.Thread(target=status_loop, args=(st, stop), daemon=True)
    t.start()

    results = []
    no_rna = args.no_rna_supp

    def worker(e):
        with st["lock"]:
            st["active"].add(e["gse"])
        def cb(d, final=False):
            mb = d / (1024 * 1024)
            qlog(f"    · {e['gse']} 下载中 {mb:.1f} MB" + ("" if not final else " ✓"))
        try:
            res = download_one(e["gse"], (e["rna"] and not no_rna),
                               args.out_root, no_rna_supp=no_rna, cb=cb)
        finally:
            with st["lock"]:
                st["active"].discard(e["gse"])
                st["done"] += 1
                if res.get("matrix") and os.path.exists(res["matrix"]):
                    st["bytes"] += os.path.getsize(res["matrix"])
        # 每文件完成，单独一行汇总（经队列，避免与进度条错行）
        if res["ok"]:
            sup = ""
            if res["supp"] and isinstance(res["supp"], str) and res["supp"].startswith("("):
                sup = f"  supp跳过={res['supp']}"
            elif res["supp"]:
                sup = "  supp=OK"
            qlog(f"  [OK]   {res['gse']}: matrix={os.path.basename(res['matrix'])}"
                 f" ({res['size']/1024/1024:.1f} MB){sup}")
        else:
            qlog(f"  [FAIL] {res['gse']}: {res['error']}")
        return res

    with ThreadPoolExecutor(max_workers=max(1, min(args.workers, 4))) as ex:
        futs = {ex.submit(worker, e): e for e, _, _ in plan}
        for fut in as_completed(futs):
            results.append(fut.result())

    stop.set()
    t.join()

    # 清单输出
    ok = [r for r in results if r["ok"]]
    fail = [r for r in results if not r["ok"]]
    manifest = {
        "total": len(results), "ok": len(ok), "fail": len(fail),
        "ok_list": [r["gse"] for r in ok],
        "fail_list": [{"gse": r["gse"], "error": r["error"]} for r in fail],
    }
    mp = os.path.join(args.out_root, "download_manifest.json")
    with open(mp, "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    eprint("")
    if fail:
        fp = os.path.join(args.out_root, "download_failed.txt")
        with open(fp, "w", encoding="utf-8") as f:
            for r in fail:
                f.write(r["gse"] + "\n")
        eprint(f"[完成] 成功 {len(ok)}/{len(results)}；失败 {len(fail)} 个 -> {fp}")
        eprint(f"        成功清单 -> {mp}")
    else:
        eprint(f"[完成] 全部成功 {len(ok)}/{len(results)}；清单 -> {mp}")


if __name__ == "__main__":
    main()
