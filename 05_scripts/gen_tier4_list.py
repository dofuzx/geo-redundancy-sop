#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tier4（领域级）GEO 系列清单程序化生成器

用 NCBI E-utilities 枚举结直肠癌领域的 GPL570/GPL96 系列，剔除已用于 Tier1/2/3 的编号，
产出可直接喂给 download_tier3_4.py --list 的清单。

为什么必须程序化：手工抄写 200+ 个 GSE 编号极易出错；且领域名随时间变动，
程序化清单可重跑刷新，保证可复现。

用法：
    python gen_tier4_list.py                      # 默认产出 250 条到 --out
    python gen_tier4_list.py --target 300         # 目标条数
    python gen_tier4_list.py --platforms 570 96   # 指定平台
NCBI 礼貌限速：无 API key 时 3 req/s，脚本已内置 sleep；可用 --api-key 提速。
"""
import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# 已用系列目录（自动生成排除集）
USED_DIRS = [
    ("D:/AIwork/GEO", "_excluded_rna"),
    ("D:/AIwork/GEO_sop_v2", None),
    ("E:/CRC_GEO_2026T1_sop_v2", None),
    ("E:/CRC_GEO_2026T2_sop_v2", None),
    ("E:/CRC_GEO_2026T1_sop", None),
    ("E:/CRC_GEO_2026T2_sop", None),
]
# 之前已知下载失败/损坏的编号（也可放行重下，由 --allow-failed 控制）
KNOWN_FAILED = ["GSE13471", "GSE110225"]

KEYWORDS = [
    "colorectal", "colon cancer", "colon adenocarcinoma", "colonic",
    "rectal", "rectum", "colon carcinoma", "colorectal cancer",
    "colorectal carcinoma", "colitis-associated",
    # 补充召回词：部分系列标题使用别称，不补会漏掉大量候选
    "colon polyp", "colonic mucosa", "rectosigmoid", "large bowel",
    "colorectal liver metastasis", "colon adenoma", "rectal cancer",
    "colonic epithelium", "colonic biopsy", "colon", "serrated",
]

SEARCH_TERMS = [
    '"{kw}"[All Fields]',
]

_last = [0.0]


def _throttle(min_interval=0.34):
    gap = time.time() - _last[0]
    if gap < min_interval:
        time.sleep(min_interval - gap)
    _last[0] = time.time()


def _get(url, retries=5, timeout=120):
    """网络请求（含重试与退避）。NCBI 大批量 esummary 偶发读超时，必须重试。"""
    last = None
    for i in range(retries):
        try:
            _throttle()
            req = urllib.request.Request(url)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.load(r)
        except Exception as e:
            last = e
            print(f"    [重试 {i+1}/{retries}] {type(e).__name__}: {e}", flush=True)
            time.sleep(2.0 * (i + 1))
    raise last


def esearch(term, retmax=10000, api_key=None):
    q = {"db": "gds", "retmode": "json", "retmax": retmax, "term": term}
    if api_key:
        q["api_key"] = api_key
    return _get(f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(q))


def esummary(uids, api_key=None, batch=150):
    out = []
    for i in range(0, len(uids), batch):
        chunk = uids[i:i + batch]
        q = {"db": "gds", "retmode": "json", "retmax": batch, "id": ",".join(chunk)}
        if api_key:
            q["api_key"] = api_key
        try:
            j = _get(f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(q))
        except Exception as e:
            print(f"    [跳过] 批次 {i//batch + 1} 彻底失败: {type(e).__name__}", flush=True)
            continue
        res = j.get("result", {})
        for u in res.get("uids", []):
            out.append(res[u])
        print(f"    esummary 批次 {i//batch + 1}/{(len(uids)+batch-1)//batch}: +{len(out)} 累计", flush=True)
    return out


def collect_used():
    """从既用语料目录扫描出已用 GSE 编号"""
    used = set()
    for d, skip_sub in USED_DIRS:
        if not os.path.isdir(d):
            continue
        for fn in os.listdir(d):
            if fn.endswith(".csv"):
                used.add(fn[:-4].split("_GPL")[0])
    # 源数据目录也扫一遍（含 RNA-seq / 被剔除的）
    for root in ["D:/AIwork/GEO", "E:/CRC_GEO_2026T1", "E:/CRC_GEO_2026T2"]:
        if not os.path.isdir(root):
            continue
        for name in os.listdir(root):
            m = __import__("re").match(r"(?:.*_)?(GSE\d+)", name)
            if m:
                used.add(m.group(1))
    return used


def main():
    ap = argparse.ArgumentParser(description="生成 Tier4 领域级 GEO 下载清单")
    ap.add_argument("--target", type=int, default=250, help="目标条数(默认250)")
    ap.add_argument("--platforms", nargs="+", default=["570", "96"],
                    help="平台 GPL 号，默认 570 96")
    ap.add_argument("--min-samples", type=int, default=6, help="最小样本数过滤")
    ap.add_argument("--out", default="crc_gse_tier4_full.txt", help="清单输出文件")
    ap.add_argument("--detail", default="crc_gse_tier4_detail.tsv", help="明细输出文件")
    ap.add_argument("--api-key", default=None, help="NCBI API key(可选,提速)")
    ap.add_argument("--allow-failed", action="store_true", help="不排除已知失败编号")
    ap.add_argument("--require-title", action="store_true",
                    help="只保留标题/ssinfo 命中 CRC 关键词的高置信项(默认也保留中置信项以达到 200+)")
    args = ap.parse_args()

    used = collect_used()
    print(f"[信息] 已用系列 {len(used)} 个，将剔除（示例: {sorted(used)[:5]} ...）", flush=True)

    # 1) 多轮关键词搜索并合并 uid
    uid_set = set()
    for plat in args.platforms:
        for kw in KEYWORDS:
            term = f'GPL{plat}[Platform] AND "GSE"[Entry Type] AND "{kw}"[All Fields]'
            try:
                j = esearch(term, api_key=args.api_key)
            except Exception as e:
                print(f"  [警告] 检索失败 plat={plat} kw={kw}: {type(e).__name__}", flush=True)
                continue
            res = j.get("esearchresult", {})
            ids = res.get("idlist", [])
            new = set(ids) - uid_set
            uid_set |= set(ids)
            print(f"  GPL{plat:4s} / {kw:24s} count={res.get('count','?'):>6} 新增 {len(new)}", flush=True)

    print(f"[信息] 候选 uid 合计 {len(uid_set)}，开始取摘要...", flush=True)
    rows = esummary(sorted(uid_set), api_key=args.api_key)

    # 2) 过滤
    keep = {}
    excl_super = []
    for d in rows:
        acc = d.get("accession", "")
        if not acc.startswith("GSE"):
            continue
        if d.get("entrytype") != "GSE":
            continue
        plat = str(d.get("gpl", ""))
        if plat not in args.platforms:
            continue
        n = int(d.get("n_samples", 0) or 0)
        if n < args.min_samples:
            continue
        if acc in used:
            continue
        if (not args.allow_failed) and acc in KNOWN_FAILED:
            continue
        title = str(d.get("seriestitle") or d.get("title") or "")
        rel = str(d.get("relations") or "")
        # SuperSeries 会重复包含自己的 SubSeries，排除以免自我重复计数
        if "SuperSeries" in rel or "superseries" in title.lower():
            excl_super.append((acc, n, title))
            continue
        # 相关性分级（不硬性丢弃中置信项，否则凑不到领域级 200+ 规模）：
        #   high   = 标题/ssinfo 命中 CRC 关键词
        #   medium = 仅 esearch All Fields 命中（标题用了别称，需人工看一眼）
        ctx = (title + " " + str(d.get("ssinfo") or "")).lower()
        conf = "high" if any(k.lower() in ctx for k in KEYWORDS) else "medium"
        if conf == "medium" and args.require_title:
            continue
        if acc in keep:
            continue
        keep[acc] = {"gse": acc, "gpl": plat, "n_samples": n, "conf": conf,
                     "title": title.replace("\t", " ")[:120],
                     "pdat": str(d.get("pdat", "")),
                     "pubmed": ",".join(str(x) for x in (d.get("pubmedids") or [])[:3])}

    items = sorted(keep.values(), key=lambda x: -x["n_samples"])
    sel = items[:args.target]

    with open(args.out, "w", encoding="utf-8") as f:
        f.write("# Tier4 领域级 CRC GEO 下载清单\n")
        f.write(f"# 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# 平台: {','.join('GPL'+p for p in args.platforms)} | 最小样本数: {args.min_samples}\n")
        f.write(f"# 候选 {len(rows)} -> 通过过滤 {len(items)} -> 本清单输出 {len(sel)}\n")
        f.write("# 已剔除 Tier1/2/3 既用编号；剔除 SuperSeries（避免与其 SubSeries 自我重复）\n")
        for it in sel:
            f.write(it["gse"] + "\n")

    with open(args.detail, "w", encoding="utf-8") as f:
        f.write("GSE\tGPL\tn_samples\tpubmed\tdate\ttitle\n")
        for it in sel:
            f.write(f"{it['gse']}\tGPL{it['gpl']}\t{it['n_samples']}\t{it['pubmed']}\t{it['pdat']}\t{it['title']}\n")

    print(f"\n[完成] 候选 {len(rows)} -> 通过 {len(items)} -> 输出 {len(sel)} 条", flush=True)
    print(f"  清单: {args.out}")
    print(f"  明细: {args.detail}")
    print(f"  被排除的 SuperSeries {len(excl_super)} 个: "
          f"{', '.join(a for a,_,_ in excl_super[:8])}{' ...' if len(excl_super)>8 else ''}", flush=True)
    tot = sum(it["n_samples"] for it in sel)
    print(f"  清单样本总量估计: {tot}")


if __name__ == "__main__":
    main()
