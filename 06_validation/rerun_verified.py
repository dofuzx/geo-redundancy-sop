#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rerun_verified.py — 在本机 CSV 语料上重跑补充包的两个工作流，核验/更正其数值。

背景：补充包 threshold_robustness.py 的 pair_link_report 存在索引 bug
（C[:nA, nB:] 在非方阵时错位成子块），已修复；本脚本用修复后的内核 +
本机 CSV 数据重跑，产出权威数值：
  - 三个锚点对（2 个阳性 + 1 个特异性）的链路复现/门限特异性
  - 4 个阴性对照的 FPR
  - spike-in 扰动基准（GSE203024, 100 样本, seed=42）
"""
import csv
import os
import sys

import numpy as np

sys.path.insert(0, r"G:/YYG/supplement")
import dedup_lib as dl  # noqa: E402
from corpus_sweep import read_matrix  # noqa: E402

Series = dl.Series
CSV_DIRS = [r"G:/YYG/data/GEO_Tire3_sop", r"G:/YYG/data/GEO_Tire4_sop",
            r"J:/CRC_GEO_2026T1_sop_v2", r"J:/CRC_GEO_2026T2_sop_v2",
            r"J:/CRC_GEO_2026T1_sop_v2/_excluded_structural"]
OUT = r"G:/YYG/supplement/sweep_output"


def find_series_matrix_csv(root, gse):
    import re
    g = re.sub(r"_GPL\d+$", "", gse)
    for d in CSV_DIRS:
        p = os.path.join(d, g + ".csv")
        if os.path.exists(p):
            return p
    return None


def build_corpus_index_csv(root):
    idx = {}
    for d in CSV_DIRS:
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn.endswith(".csv"):
                m = __import__("re").search(r"(GSE\d+)", fn)
                if m:
                    idx.setdefault(m.group(1), os.path.join(d, fn))
    return idx


def load_series_matrix_csv(path):
    genes, samples, M = read_matrix(path)
    return Series(genes=genes, samples=samples, M=M, path=path)


# ---- 打补丁：dedup_lib 的加载器换成本机 CSV 版 ----
dl.load_series_matrix = load_series_matrix_csv
dl.find_series_matrix = find_series_matrix_csv
dl.build_corpus_index = build_corpus_index_csv

FLOOR, Q = 0.97, 99.9
POS = [("GSE41258", "GSE68468"), ("GSE12251", "GSE23597")]
NEG = ("GSE29621", "GSE17536")
FPR_PAIRS = [("GSE106325", "GSE59196"), ("GSE20916", "GSE71571"),
             ("GSE5816", "GSE23878"), ("GSE41326", "GSE11440")]


def pair_report(ga, gb):
    sa = load_series_matrix_csv(find_series_matrix_csv(None, ga))
    sb = load_series_matrix_csv(find_series_matrix_csv(None, gb))
    return dl.pair_link_report(sa, sb, floor=FLOOR, q=Q)


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    print("=== (A) 链路复现（阳性锚点）===")
    for ga, gb in POS:
        r = pair_report(ga, gb)
        r["pair"] = "%s<->%s" % (ga, gb)
        r["section"] = "link_recovery"
        rows.append(r)
        print("  %s: nA=%d nB=%d common=%d bA=%.4f bB=%.4f T=%.4f peak=%.4f naive=%d gated=%d recall=%.3f"
              % (r["pair"], r["nA"], r["nB"], r["common_genes"], r["baseline_A"], r["baseline_B"],
                 r["threshold"], r["cross_max"], r["naive_pairs_above_floor"],
                 r["gated_pairs"], r["recall_at_T"]))

    print("=== (B) 门限特异性（假阳性塌缩对）===")
    r = pair_report(*NEG)
    r["pair"] = "%s<->%s" % NEG
    r["section"] = "specificity"
    rows.append(r)
    print("  %s: naive=%d gated=%d T=%.4f peak=%.4f  (稿件: 2825 -> 66)"
          % (r["pair"], r["naive_pairs_above_floor"], r["gated_pairs"], r["threshold"], r["cross_max"]))

    print("=== (C) 阴性对照 FPR ===")
    fpr_rows = []
    for ga, gb in FPR_PAIRS:
        sa = load_series_matrix_csv(find_series_matrix_csv(None, ga))
        sb = load_series_matrix_csv(find_series_matrix_csv(None, gb))
        o = dl.pair_fpr(sa, sb, floor=FLOOR, q=Q)
        o["pair"] = "%s<->%s" % (ga, gb)
        fpr_rows.append(o)
        print("  %s: cross_max=%.4f recall@T=%.3f fp=%s" % (o["pair"], o["cross_max"], o["recall_at_T"], "no" if o["recall_at_T"]==0 else "YES"))
    fpr = sum(1 for x in fpr_rows if x["recall_at_T"] > 0) / len(fpr_rows)
    print("  FPR = %.3f" % fpr)

    with open(os.path.join(OUT, "threshold_robustness_verified.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["section", "pair", "nA", "nB", "common_genes", "baseline_A", "baseline_B",
                    "threshold", "cross_max_peak_r", "naive_pairs_above_floor", "gated_pairs", "recall_at_T"])
        for r in rows:
            w.writerow([r["section"], r["pair"], r["nA"], r["nB"], r["common_genes"], r["baseline_A"],
                        r["baseline_B"], r["threshold"], r["cross_max"],
                        r["naive_pairs_above_floor"], r["gated_pairs"], r["recall_at_T"]])
        w.writerow(["fpr_negative", "", "", "", "", "", "", "", "", "", "", ""])
        for o in fpr_rows:
            w.writerow(["fpr_negative", o["pair"], "", "", "", "", "", "", o["cross_max"], "", "", o["recall_at_T"]])
        w.writerow(["FPR", "", "", "", "", "", "", "", "", "", "", fpr])

    print("=== (D) spike-in 扰动基准（GSE203024, 100 样本, seed=42）===")
    import spike_in_benchmark as sib
    gse, srows = sib.run_real(None, "GSE203024", 100, 42)
    with open(os.path.join(OUT, "spike_in_results_verified.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["perturbation", "level", "recall", "threshold", "mean_best_r", "min_best_r"])
        w.writeheader()
        for r in srows:
            w.writerow(r)
    sib.plot(gse, srows, os.path.join(OUT, "spike_in_recall_verified.png"))
    print("完成 -> %s" % OUT)


if __name__ == "__main__":
    main()
