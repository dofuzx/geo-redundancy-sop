#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_positive_control.py — 用与稿件检测器 1:1 一致的内核，在真实 GEO 数据上
复算 GSE41258 ↔ GSE68468 与 GSE12251 ↔ GSE23597 两个阳性对照的全部引用数值。

输出（供 V4 稿件 3.5 / 3.6 节与 Figure 5 图注使用）：
  - 双方样本数、共同基因数
  - 互向最佳匹配 r ≥ 0.95 的样本数、峰值 r、均值
  - 绝对 r≥0.97 对数、门限后对数、门限值
"""
import sys
import numpy as np

sys.path.insert(0, r"G:/YYG/supplement")
from corpus_sweep import read_matrix  # noqa: E402

PAIRS = [
    ("GSE41258", "GSE68468", r"J:/CRC_GEO_2026T1_sop_v2"),
    ("GSE12251", "GSE23597", None),
]


def find(g):
    import os
    for d in [r"G:/YYG/data/GEO_Tire3_sop", r"G:/YYG/data/GEO_Tire4_sop",
              r"J:/CRC_GEO_2026T1_sop_v2", r"J:/CRC_GEO_2026T2_sop_v2"]:
        p = d + "/" + g + ".csv"
        if os.path.exists(p):
            return p
    raise FileNotFoundError(g)


def centre(M):
    col_mean = np.nanmean(M, axis=0, keepdims=True)
    Mc = M - col_mean
    np.nan_to_num(Mc, copy=False)
    nrm = np.sqrt((Mc ** 2).sum(axis=0, keepdims=True))
    Mn = Mc
    Mn /= np.where(nrm == 0, 1, nrm)
    return Mn


def analyse(g1, g2):
    genes1, s1, M1 = read_matrix(find(g1))
    genes2, s2, M2 = read_matrix(find(g2))
    common = sorted(set(genes1) & set(genes2))
    i1 = {g: i for i, g in enumerate(genes1)}
    i2 = {g: i for i, g in enumerate(genes2)}
    r1 = [i1[g] for g in common]
    r2 = [i2[g] for g in common]
    A = centre(np.ascontiguousarray(M1[r1, :]))
    B = centre(np.ascontiguousarray(M2[r2, :]))
    C = np.clip(A.T @ B, -1, 1)                       # nA × nB
    bA = float(np.percentile(C[np.triu_indices(C.shape[0], 1)] if False else
                             np.clip(A.T @ A, -1, 1)[np.triu_indices(len(s1), 1)], 99.9))
    bB = float(np.percentile(np.clip(B.T @ B, -1, 1)[np.triu_indices(len(s2), 1)], 99.9))
    gate = max(0.97, bA, bB)
    naive97 = int(np.count_nonzero(C >= 0.97))
    gated = int(np.count_nonzero(C >= gate))
    # 互向最佳匹配
    a2b = C.max(axis=1)          # A 每样本在 B 的最佳匹配
    b2a = C.max(axis=0)          # B 每样本在 A 的最佳匹配
    print("=" * 70)
    print("%s (%d samples) <-> %s (%d samples) | common genes %d" % (g1, len(s1), g2, len(s2), len(common)))
    print("  intra p99.9: A=%.4f  B=%.4f  -> gate=%.4f" % (bA, bB, gate))
    print("  peak cross r = %.4f" % C.max())
    print("  naive pairs r>=0.97: %d | gated pairs: %d" % (naive97, gated))
    print("  A->B best match: >=0.95: %d/%d | >=0.97: %d | mean(best)=%.4f"
          % (int((a2b >= 0.95).sum()), len(s1), int((a2b >= 0.97).sum()), a2b.mean()))
    print("  B->A best match: >=0.95: %d/%d | >=0.97: %d | mean(best)=%.4f"
          % (int((b2a >= 0.95).sum()), len(s2), int((b2a >= 0.97).sum()), b2a.mean()))
    print("  mean r over all cross pairs = %.4f" % C.mean())
    print("  pairs r>=0.95: %d | r>=0.9999: %d" % (int((C >= 0.95).sum()), int((C >= 0.9999).sum())))


if __name__ == "__main__":
    for g1, g2, _ in PAIRS:
        analyse(g1, g2)
