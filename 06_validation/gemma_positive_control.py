#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
工作流 1 — 外部真实正控（回答"检测器能找回已知重复吗"）

来源：Gemma (https://gemma.msl.ubc.ca/rest/v2) —— 其策展团队会主动移除重复样本，
并在 description / curationNote / troubleDetails 里标注
"N samples were removed as duplicates of dataset X"。
这些是【外部权威、且与本稿件调参语料无关的】已知隐蔽重复，正可证伪"只在自己数据上超参"。

两阶段：
  discover : 拉取 human + colorectal/cancer 数据集（仅翻页列表，列表项已含 description /
            curationNote / troubleDetails / taxonId，无需逐条请求），解析去重注记，
            抽取被引用的 GSE -> candidates.csv
            【你】打开 CSV 人工锁定"不同 GSM 的隐蔽重复对" -> locked_pairs.csv
  evaluate : 对锁定对，优先从本地 D:/AIwork/GEO 读 series_matrix 算召回；
             本地没有的，回退到 Gemma /data/processed 算召回。报告 recall + 跨系列 r 分布。

关于网络/代理（已实测用户环境）：
  - Gemma 的 /datasets 不接受 `taxon=1` 这种参数（返回 400），必须只用 query= 然后在
    客户端按 taxonId==1 过滤。
  - 本机若设置了会掉线的系统代理（HTTPS_PROXY），请求会间歇性 ProxyError。
    脚本默认【绕过系统代理】直连；若你确实身处需要代理才能出网的内网，加 --use-system-proxy。

依赖：requests（仅 discover & Gemma 回退需要联网）。
"""
import argparse
import csv
import os
import re
import sys
import time

import numpy as np

from dedup_lib import (load_series_matrix, find_series_matrix, pair_recall,
                       column_centered_pearson)

BASE = "https://gemma.msl.ubc.ca/rest/v2"
ROOT = "D:/AIwork/GEO"

# 去重 / 重叠 注记的正则（大小写不敏感）
DUP_PATTERNS = [
    r"removed as duplicates? of",
    r"duplicate(?:s|d)? (?:of|to) (?:dataset|study|series|accession)?\s*[GSE]",
    r"same as (?:the )?(?:study|dataset|series)?\s*[GSE]",
    r"overlapp(?:s|ing|ed)? with [GSE]",
    r"redundant (?:with|to) [GSE]",
]
GSE_RE = re.compile(r"GSE\d+")


# ---------------------------------------------------------------------------
# 网络：会话 + 重试 + 代理开关
# ---------------------------------------------------------------------------
def _make_session(use_system_proxy=False):
    import requests
    s = requests.Session()
    if not use_system_proxy:
        # 绕过系统代理（实测用户机器的 HTTPS_PROXY 会间歇性掉线）
        s.trust_env = False
        s.proxies = {"http": None, "https": None}
    s.headers.update({"Accept": "application/json",
                      "User-Agent": "dedup-validation/1.0"})
    return s


def _get(session, path, params=None, retries=5):
    import requests
    last = None
    for attempt in range(retries):
        try:
            r = session.get(BASE + path, params=params, timeout=60)
            r.raise_for_status()
            return r.json()
        except (requests.exceptions.ProxyError,
                requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.HTTPError) as e:
            last = e
            # HTTPError(如 500) 也重试；400 通常是参数错，但这里统一重试几次更稳
            time.sleep(1.5 + attempt * 2.0)
    if last:
        raise last
    raise RuntimeError("unknown network error")


# ---------------------------------------------------------------------------
# discover
# ---------------------------------------------------------------------------
def discover(query="colorectal", out_csv="candidates.csv", limit=400,
             use_system_proxy=False, page_size=100, sleep=0.1):
    session = _make_session(use_system_proxy)
    print("[discover] 检索 human + %r（翻页列表，客户端过滤 taxonId==1）..." % query)
    rows = []
    seen = set()
    offset = 0
    fetched = 0
    while fetched < limit:
        data = _get(session, "/datasets",
                    params={"query": query, "offset": offset, "limit": page_size})
        chunk = data.get("data", [])
        if not chunk:
            break
        for d in chunk:
            acc = d.get("accession") or d.get("shortName")
            if not acc or acc in seen:
                continue
            seen.add(acc)
            if d.get("taxonId") != 1:        # 只保留人类
                continue
            desc = d.get("description") or ""
            note = d.get("curationNote") or ""
            trbl = d.get("troubleDetails") or ""
            blob = " ".join([desc, note, trbl])
            hit = [p for p in DUP_PATTERNS if re.search(p, blob, re.I)]
            if not hit:
                continue
            refs = sorted(set(GSE_RE.findall(blob)))
            refs = [g for g in refs if g != acc]   # 去掉自己
            rows.append({
                "accession": acc,
                "name": (d.get("name") or "")[:120],
                "bioAssayCount": d.get("bioAssayCount"),
                "dup_pattern": ";".join(hit),
                "referenced_gse": ";".join(refs),
                "note_excerpt": blob[:500],
            })
        fetched += len(chunk)
        offset += page_size
        if len(chunk) < page_size:
            break
        time.sleep(sleep)
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "name", "bioAssayCount",
                                          "dup_pattern", "referenced_gse",
                                          "note_excerpt"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("[discover] 命中去重/重叠注记 %d 条（人类）-> 写入 %s" % (len(rows), out_csv))
    print("[discover] 下一步：【你】打开 CSV，人工锁定真正“不同 GSM 的隐蔽重复对”，")
    print("           整理成 locked_pairs.csv，列: accA,accB（每行一个已知重复数据集对）。")
    return rows


# ---------------------------------------------------------------------------
# evaluate
# ---------------------------------------------------------------------------
def _load_gemma_processed(session, acc):
    """回退：从 Gemma /data/processed 拉 log2 表达矩阵 (TSV)。返回 (genes, samples, M)。"""
    r = session.get("%s/datasets/%s/data/processed" % (BASE, acc), timeout=120)
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
    M = np.array(rows, dtype=np.float32)
    return genes, samples, M


def evaluate(locked_csv="locked_pairs.csv", out_csv="gemma_eval.csv",
             use_system_proxy=False):
    session = _make_session(use_system_proxy)
    results = []
    with open(locked_csv, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            a = row.get("accA", "").strip()
            b = row.get("accB", row.get("b", "")).strip()
            if not a or not b:
                continue
            pa = find_series_matrix(ROOT, a)
            pb = find_series_matrix(ROOT, b)
            src = "local"
            try:
                if pa and pb:
                    sa = load_series_matrix(pa)
                    sb = load_series_matrix(pb)
                else:
                    ga, sga, Ma = _load_gemma_processed(session, a)
                    gb, sgb, Mb = _load_gemma_processed(session, b)
                    from dedup_lib import Series
                    sa = Series(ga, sga, Ma, a)
                    sb = Series(gb, sgb, Mb, b)
                    src = "gemma_api"
                o = pair_recall(sa, sb, floor=0.97)
                o.update({"accA": a, "accB": b, "source": src})
                results.append(o)
                print("[evaluate] %s↔%s [%s] recall@0.97=%.3f cross_max=%.4f"
                      % (a, b, src, o["recall"], o["cross_max"]))
            except Exception as e:
                print("[evaluate] %s↔%s 失败: %s" % (a, b, e))
                results.append({"accA": a, "accB": b, "source": src,
                                "error": str(e)})
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accA", "accB", "source", "recall",
                                          "threshold", "cross_max", "mean_best_A2B",
                                          "nA", "nB", "error"])
        w.writeheader()
        for r in results:
            w.writerow(r)
    recs = [r["recall"] for r in results if "recall" in r]
    if recs:
        print("[evaluate] 总体召回率(已知隐蔽重复被找回比例) = %.3f (%d 对)"
              % (np.mean(recs), len(recs)))
    print("[evaluate] 结果 -> %s" % out_csv)


def main():
    ap = argparse.ArgumentParser(description="Gemma 外部真实正控")
    ap.add_argument("stage", choices=["discover", "evaluate"])
    ap.add_argument("--query", default="colorectal")
    ap.add_argument("--limit", type=int, default=400)
    ap.add_argument("--locked", default="locked_pairs.csv")
    ap.add_argument("--out", default=None)
    ap.add_argument("--use-system-proxy", action="store_true",
                   help="默认绕过系统代理直连；仅在你身处需代理才能出网的内网时加此开关")
    args = ap.parse_args()
    if args.stage == "discover":
        discover(query=args.query, out_csv="candidates.csv", limit=args.limit,
                 use_system_proxy=args.use_system_proxy)
    else:
        evaluate(locked_csv=args.locked, out_csv=args.out or "gemma_eval.csv",
                 use_system_proxy=args.use_system_proxy)


if __name__ == "__main__":
    main()
