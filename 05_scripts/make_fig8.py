#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""实验B·Fig8：RNA-seq 扩展结果图（双面板）。
(a) 绝对阈值 vs 基线感知（伪相关削减）——与 Fig3 对应的平台迁移版
(b) 组内基线分布（RNA-seq vs microarray 对照）——展示判据自适应价值
"""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RJ = sys.argv[1] if len(sys.argv) > 1 else r"D:/AIwork/GEO_RNAseq/out/dedup_report_rnaseq.json"
MJJ = sys.argv[2] if len(sys.argv) > 2 else r"D:/AIwork/GEO_T1234_final/dedup_report.json"
OUT = sys.argv[3] if len(sys.argv) > 3 else r"D:/AIwork/GEO_RNAseq/out/Fig8_rnaseq_extension.png"

r = json.load(open(RJ, encoding="utf-8"))
s = r["summary"]
pair_stats = r["layerC_pair_stats"]

fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

# ---------- (a) 绝对 vs 基线感知（TOP 10 系列对） ----------
ax = axes[0]
top = pair_stats[:10]
if top:
    labels = [f"{p['g1'][3:]}↔{p['g2'][3:]}" for p in top]
    x = np.arange(len(top))
    ax.bar(x - 0.2, [p["edges_abs_0.97"] for p in top], width=0.4,
           color="#c0392b", edgecolor="black", linewidth=0.4, label="absolute ρ≥0.97")
    ax.bar(x + 0.2, [p["edges_baseline_aware"] for p in top], width=0.4,
           color="#27ae60", edgecolor="black", linewidth=0.4, label="baseline-aware")
    ax.set_yscale("symlog")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("cross-series sample-pairs")
ax.legend(fontsize=8)
ax.set_title(f"(a) RNA-seq: absolute vs baseline-aware\n"
             f"({s['edges_absolute_0.97_total']} → {s['edges_baseline_aware_total']} pairs)", fontsize=10)

# ---------- (b) 组内基线：RNA-seq vs microarray ----------
ax = axes[1]
try:
    m = json.load(open(MJJ, encoding="utf-8"))
    if "series_baselines" in m:
        ma_base = [v for v in m["series_baselines"].values() if v and v > 0]
    else:   # microarray 报告无逐系列基线，从系列对统计提取
        q = "baseline_p%.4g_g1" % m["summary"]["baseline_q"]
        q2 = "baseline_p%.4g_g2" % m["summary"]["baseline_q"]
        seen = {}
        for p in m["layerC_pair_stats"]:
            if p.get(q): seen[p["g1"]] = p[q]
            if p.get(q2): seen[p["g2"]] = p[q2]
        ma_base = list(seen.values())
except Exception:
    ma_base = []
rna_base = [v for v in r.get("series_baselines", {}).values() if v and v > 0]
bp = ax.boxplot([ma_base, rna_base], tick_labels=["microarray\n(GPL570+)", "RNA-seq\n(this study)"],
                widths=0.5, patch_artist=True, showfliers=False)
for patch, c in zip(bp["boxes"], ["#95a5a6", "#2980b9"]):
    patch.set_facecolor(c)
    patch.set_alpha(0.7)
ax.set_ylabel("intra-series p99.9 baseline (per series)")
med_r = float(np.median(rna_base)) if rna_base else 0
med_m = float(np.median(ma_base)) if ma_base else 0
ax.axhline(0.97, ls="--", color="gray", lw=1)
ax.text(1.42, 0.972, "floor 0.97", fontsize=8, color="gray", ha="right")
ax.set_title(f"(b) Intra-series baselines (p99.9)\nmedian: microarray {med_m:.3f} / RNA-seq {med_r:.3f}", fontsize=10)

fig.suptitle("Experiment B — SOP transfers to RNA-seq with the same safeguards", fontsize=11, y=1.00)
fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight")
print("saved:", OUT)
print(f"summary: series={s['n_series_used']} samples={s['n_samples']} "
      f"abs={s['edges_absolute_0.97_total']} ada={s['edges_baseline_aware_total']} "
      f"rate={s['redundancy_rate']} baselines(med RNA)={med_r:.4f} (med MA)={med_m:.4f}")
