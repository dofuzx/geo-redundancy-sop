#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""make_figure_S1.py — 用 corpus_sweep.py 的真实扫描结果重画 Figure S1 热图。

输入: sweep_output/sweep_matrix.csv (tau x percentile -> redundancy rate %)
输出: Figure_S1_threshold_sensitivity.png (300 dpi, 浅色主题)
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = r"G:/YYG/supplement/sweep_output"
PNG = os.path.join(OUT, "Figure_S1_threshold_sensitivity.png")

rows = list(csv.DictReader(open(os.path.join(OUT, "sweep_matrix.csv"), encoding="utf-8")))
taus = [float(r["tau"]) for r in rows]
qs = [c for c in rows[0].keys() if c != "tau"]
M = np.array([[float(r[q]) for q in qs] for r in rows])          # tau x q

fig, ax = plt.subplots(figsize=(7.4, 4.3), dpi=300)
im = ax.imshow(M, cmap="YlOrRd", aspect="auto", vmin=7.5, vmax=13.0)
ax.set_xticks(range(len(qs)))
ax.set_xticklabels([q.replace("p", "p") for q in qs], fontsize=10)
ax.set_yticks(range(len(taus)))
ax.set_yticklabels(["%.2f" % t for t in taus], fontsize=10)
ax.set_xlabel("Intra-series baseline percentile (q)", fontsize=11)
ax.set_ylabel("Correlation threshold (\u03c4)", fontsize=11)
ax.set_title("Corpus redundancy rate (%) across the sensitivity grid\n"
             "319 CRC series / 25,460 samples; baseline-aware gate = max(\u03c4, p$q$ baseline$_A$, p$q$ baseline$_B$)",
             fontsize=11)

for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v = M[i, j]
        ax.text(j, i, "%.2f" % v, ha="center", va="center", fontsize=10,
                color="white" if v > 10.8 else "#7f1d1d",
                fontweight="bold" if (i == 2 and j == 2) else "normal")

# 标注操作点 (tau=0.97, q=99.9)
oi, oj = taus.index(0.97), qs.index("p99.9")
ax.add_patch(plt.Rectangle((oj - 0.5, oi - 0.5), 1, 1, fill=False,
                           edgecolor="#1d4ed8", lw=2.2))
ax.text(oj, oi - 0.62, "operating point\n8.52%", ha="center", va="bottom",
        fontsize=8.5, color="#1d4ed8", fontweight="bold")

cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
cb.set_label("Redundancy rate (%)", fontsize=10)
cb.ax.tick_params(labelsize=9)
ax.set_xticks(np.arange(-0.5, len(qs), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(taus), 1), minor=True)
ax.grid(which="minor", color="white", lw=1.4)
ax.tick_params(which="minor", length=0)
fig.tight_layout()
fig.savefig(PNG, dpi=300, facecolor="white")
print("saved", PNG)
