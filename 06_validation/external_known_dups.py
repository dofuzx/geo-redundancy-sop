#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
工作流 1（可靠版）— 外部真实正控：基于 GEO/Gemma 样本 ID 重叠的已知重复对

为什么不用 Gemma 自由文本去重注记：
  Gemma 的 API 文本字段（description / curationNote / troubleDetails / auditEvents）
  里【没有】干净的“数据集 X 是 GSExxxx 的重复”成对标注。实测：295 个含单词
  “duplicate” 的数据集中，0 条在注记里出现 “duplicate/redundant + 引用 GSE”。
  Gemma 的去重是样本级内部策展，不对外暴露成对链接。故改为下面更可靠的方案。

本方案的独立真值来源：
  GEO/Gemma 的【样本 ID（GSM）】是第三方权威标识，与你的调参语料无关。
  若两个不同 GSE 共享大量完全相同的 GSM，则 GEO 自己就“确认”它们是同一批样本
  （重复提交 / 再分析），这就是【独立已知重复对】。

测试目标（直接回应审稿人）：
  检验你的【表达指纹层 Layer C/D】（非 trivial 的 GSM 匹配层）能否【独立】找回这些
  已知重复对——尤其当两 GSE【不同平台】（GSM 前缀都不同）时，GSM 匹配层根本帮不上，
  只能靠表达指纹，正测你被质疑的那一层。

两阶段：
  discover : 拉取独立 CRC/colon 队列（排除本地 Tire1-4，确保与调参无关），取每数据集
             GSM 集合 + 平台；算 GSM-Jaccard，挑出重叠对 -> external_dup_candidates.csv
             【你】可人工复核该 CSV（确认确为同批样本重复提交）。
  evaluate : 对候选对，优先读本地 series_matrix，否则回退 Gemma /data/processed；
             用 Layer C/D（pair_recall，纯表达指纹）算召回 -> external_dup_eval.csv

网络/代理：
  默认绕过系统代理直连（实测用户机器 HTTPS_PROXY 会间歇性掉线）；身处需代理内网时
  加 --use-system-proxy。
"""
import argparse
import csv
import json
import os
import re
import sys
import time

import numpy as np

from dedup_lib import (load_series_matrix, find_series_matrix, pair_recall,
                       column_centered_pearson, Series)

BASE = "https://gemma.msl.ubc.ca/rest/v2"
ROOT = "D:/AIwork/GEO"


# ---------------------------------------------------------------------------
# 网络会话
# ---------------------------------------------------------------------------
def _make_session(use_system_proxy=False):
    import requests
    s = requests.Session()
    if not use_system_proxy:
        s.trust_env = False
        s.proxies = {"http": None, "https": None}
    s.headers.update({"Accept": "application/json", "User-Agent": "dedup-ext/1.0"})
    return s


def _get(session, path, params=None, retries=2):
    import requests
    last = None
    for attempt in range(retries):
        try:
            r = session.get(BASE + path, params=params, timeout=(10, 30))
            r.raise_for_status()
            return r.json()
        except (requests.exceptions.ProxyError,
                requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.HTTPError) as e:
            last = e
            time.sleep(1.5 + attempt * 2.0)
    if last:
        raise last
    raise RuntimeError("network error")


def _local_accessions(root):
    """本地 Tire1-4 已有的 GSE 集合（用于排除、保证独立）。"""
    accs = set()
    for cur, _, fs in os.walk(root):
        for fn in fs:
            m = re.match(r"(GSE\d+)_series_matrix\.txt\.gz", fn)
            if m:
                accs.add(m.group(1))
    return accs


# ---------------------------------------------------------------------------
# discover：取 GSM 集合 + 平台，算重叠对
# ---------------------------------------------------------------------------
def discover(query="colorectal", out_csv="external_dup_candidates.csv",
             jaccard_min=0.5, overlap_min=20, limit=400,
             use_system_proxy=False, sleep=0.12):
    session = _make_session(use_system_proxy)
    local = _local_accessions(ROOT)
    print("[discover] 本地 Tire1-4 已排除 %d 个 GSE（保证与调参无关）" % len(local))

    # 1) 拉数据集列表（人类）
    accs = []
    seen = set()
    offset = 0
    while offset < limit:
        data = _get(session, "/datasets",
                    params={"query": query, "offset": offset, "limit": 100})
        chunk = data.get("data", [])
        if not chunk:
            break
        for d in chunk:
            a = d.get("accession")
            if not a or a in seen or d.get("taxonId") != 1 or a in local:
                continue
            seen.add(a)
            accs.append(a)
        offset += 100
        if len(chunk) < 100:
            break
    print("[discover] 独立候选数据集（人类, 排除本地）: %d 个" % len(accs))

    # 2) 取每数据集 GSM 集合 + 平台
    gsm_map = {}
    plat_map = {}
    for i, a in enumerate(accs):
        try:
            dd = _get(session, "/datasets/%s" % a)
            obj = dd.get("data", dd)
            if isinstance(obj, list):
                if not obj:          # 空 data（如未找到）-> 跳过该数据集
                    print("[discover] 跳过 %s: 返回空 data" % a)
                    continue
                obj = obj[0]
            plats = obj.get("platforms") or []
            plat = (plats[0].get("shortName") if plats and isinstance(plats[0], dict)
                    else (plats[0] if plats else ""))
            # samples：Gemma v2 用 cursor 分页（不接受 offset）
            gsms = set()
            cursor = None
            while True:
                params = {"limit": 100}
                if cursor is not None:
                    params["cursor"] = cursor
                sd = _get(session, "/datasets/%s/samples" % a, params=params)
                sdata = sd.get("data", [])
                if not sdata:
                    break
                for sm in sdata:
                    sn = sm.get("shortName") or ""
                    if sn.startswith("GSM"):
                        gsms.add(sn)
                nxt = sd.get("nextCursor")
                if not nxt:
                    break
                cursor = nxt
                time.sleep(sleep)
            gsm_map[a] = gsms
            plat_map[a] = plat
            if (i + 1) % 25 == 0:
                print("  ...已取 %d/%d（累计 GSM 集合 %d 个）" % (i + 1, len(accs), len(gsm_map)))
        except Exception as e:
            print("[discover] 跳过 %s: %s" % (a, e))
        time.sleep(sleep)

    # 3) 算重叠对
    accs2 = [a for a in accs if gsm_map.get(a)]
    rows = []
    for i in range(len(accs2)):
        for j in range(i + 1, len(accs2)):
            a, b = accs2[i], accs2[j]
            ga, gb = gsm_map[a], gsm_map[b]
            inter = ga & gb
            if len(inter) < overlap_min:
                continue
            union = ga | gb
            jac = len(inter) / len(union) if union else 0.0
            if jac < jaccard_min:
                continue
            rows.append({
                "accA": a, "accB": b,
                "platformA": plat_map.get(a, ""),
                "platformB": plat_map.get(b, ""),
                "same_platform": plat_map.get(a) == plat_map.get(b),
                "shared_samples": len(inter),
                "nA": len(ga), "nB": len(gb),
                "jaccard": round(jac, 3),
            })
    rows.sort(key=lambda r: (not r["same_platform"], -r["shared_samples"]))
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accA", "accB", "platformA", "platformB",
                                          "same_platform", "shared_samples",
                                          "nA", "nB", "jaccard"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    diff = sum(1 for r in rows if not r["same_platform"])
    print("[discover] 命中 GEO 确认重复对(独立) %d 个（其中不同平台 %d 个）-> %s"
          % (len(rows), diff, out_csv))
    print("[discover] 下一步：【你】复核 CSV，确认确为同批样本重复提交；然后运行 evaluate。")
    return rows


# ---------------------------------------------------------------------------
# evaluate：Layer C/D 召回（纯表达指纹，不依赖 GSM）
# ---------------------------------------------------------------------------
def _load_gemma_processed(session, acc):
    r = session.get("%s/datasets/%s/data/processed" % (BASE, acc), timeout=180)
    r.raise_for_status()
    lines = r.text.splitlines()
    header = lines[0].split("\t")
    samples = [h.strip() for h in header[1:]]
    genes = []; rows = []
    for ln in lines[1:]:
        p = ln.split("\t")
        genes.append(p[0].strip())
        rows.append([float(x) if x not in ("", "NA", "NaN") else np.nan
                     for x in p[1:]])
    return genes, samples, np.array(rows, dtype=np.float32)


def evaluate(cand_csv="external_dup_candidates.csv",
             out_csv="external_dup_eval.csv", use_system_proxy=False):
    session = _make_session(use_system_proxy)
    results = []
    with open(cand_csv, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            a, b = row["accA"], row["accB"]
            pa, pb = find_series_matrix(ROOT, a), find_series_matrix(ROOT, b)
            src = "local"
            try:
                if pa and pb:
                    sa = load_series_matrix(pa); sb = load_series_matrix(pb)
                else:
                    ga, sga, Ma = _load_gemma_processed(session, a)
                    gb, sgb, Mb = _load_gemma_processed(session, b)
                    sa = Series(ga, sga, Ma, a); sb = Series(gb, sgb, Mb, b)
                    src = "gemma_api"
                o = pair_recall(sa, sb, floor=0.97)
                o.update({"accA": a, "accB": b, "source": src,
                          "same_platform": row.get("same_platform"),
                          "shared_samples": row.get("shared_samples")})
                results.append(o)
                print("[evaluate] %s↔%s [%s] shared=%s samePlat=%s recall@0.97=%.3f"
                      % (a, b, src, row.get("shared_samples"),
                         row.get("same_platform"), o["recall"]))
            except Exception as e:
                print("[evaluate] %s↔%s 失败: %s" % (a, b, e))
                results.append({"accA": a, "accB": b, "source": src,
                                "error": str(e)})
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accA", "accB", "source",
                                          "same_platform", "shared_samples",
                                          "recall", "threshold", "cross_max",
                                          "mean_best_A2B", "nA", "nB", "error"])
        w.writeheader()
        for r in results:
            w.writerow(r)
    recs = [r["recall"] for r in results if "recall" in r]
    if recs:
        print("[evaluate] 总体召回(已知重复被 Layer C/D 找回) = %.3f (%d 对)"
              % (np.mean(recs), len(recs)))
    print("[evaluate] -> %s" % out_csv)


def main():
    ap = argparse.ArgumentParser(description="外部真实正控（GEO 样本 ID 重叠）")
    ap.add_argument("stage", choices=["discover", "evaluate"])
    ap.add_argument("--query", default="colorectal")
    ap.add_argument("--limit", type=int, default=400)
    ap.add_argument("--jaccard", type=float, default=0.5)
    ap.add_argument("--overlap", type=int, default=20)
    ap.add_argument("--cand", default="external_dup_candidates.csv")
    ap.add_argument("--out", default=None)
    ap.add_argument("--use-system-proxy", action="store_true")
    args = ap.parse_args()
    if args.stage == "discover":
        discover(query=args.query, out_csv=args.cand, jaccard_min=args.jaccard,
                 overlap_min=args.overlap, limit=args.limit,
                 use_system_proxy=args.use_system_proxy)
    else:
        evaluate(cand_csv=args.cand, out_csv=args.out or "external_dup_eval.csv",
                 use_system_proxy=args.use_system_proxy)


if __name__ == "__main__":
    main()
