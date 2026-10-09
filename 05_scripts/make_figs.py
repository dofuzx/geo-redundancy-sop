# -*- coding: utf-8 -*-
"""Generate publication-quality figures for the CRC GEO redundancy manuscript."""
import json, csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

OUT = "manuscript_figs"
os.makedirs(OUT, exist_ok=True)
REP = "D:/AIwork/GEO_T1234_final"
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman"],
    "font.size": 11,
    "axes.linewidth": 1.0,
    "figure.dpi": 300,
})

# ---- load data ----
d = json.load(open(f"{REP}/dedup_report.json", encoding="utf-8"))
bm = d["best_match_cross"]
n_series = d["summary"]["n_series"]
n_samples = d["summary"]["total_samples"]
n_remove = d["summary"]["samples_to_remove"]
red_rate = d["summary"]["redundancy_rate"]
n_dup_pairs = len(d["layerC_expression_duplicates"])
n_groups = len(d["duplicate_groups"])
n_near = sum(1 for x in d["layerC_expression_duplicates"] if x["pearson"] >= 0.9999)

# sample -> gse from clean_sample_list.csv
s2g = {}
with open(f"{REP}/clean_sample_list.csv", encoding="utf-8") as f:
    r = csv.DictReader(f)
    for row in r:
        s2g[row["sample"].strip().strip('"')] = row["gse"]

# platform composition from report md (verified final run)
plat = {
    "GPL570":  (246, 19806),
    "GPL96":   (36,  2266),
    "GPL10558":(2,   330),
    "GPL4133": (2,   146),
}
known_series = sum(v[0] for v in plat.values())
known_samples = sum(v[1] for v in plat.values())
others_series = n_series - known_series
others_samples = n_samples - known_samples

# calibration numbers (from final dedup_report.md, verified)
raw_edges = 34497
retained_edges = n_dup_pairs  # 2956

# =====================================================================
# FIG 1 — SOP workflow (5 layers)
# =====================================================================
fig, ax = plt.subplots(figsize=(9.2, 6.4))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc, ec="#333333", tc="white", fs=10, bold=True):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.12",
                                fc=fc, ec=ec, lw=1.4))
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", color=tc,
            fontsize=fs, fontweight="bold" if bold else "normal", wrap=True)

def arrow(x1, y1, x2, y2, style="-|>", color="#444444", lw=1.6):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                 mutation_scale=14, color=color, lw=lw))

box(0.3, 8.4, 3.0, 1.2, "GEO series input\n(expression matrix + metadata)", "#5B6770")
arrow(3.3, 9.0, 4.0, 9.0)
box(4.0, 8.4, 5.7, 1.2, "Four complementary detection layers", "#2C3E50", fs=11)

# layer boxes
ly = [6.7, 5.2, 3.7, 2.2]
cols = ["#1F77B4", "#9467BD", "#D9822B", "#2CA02C"]
labels = [
    "Layer A — GSM exact overlap\n(un-concealed duplicates)",
    "Layer B — metadata fuzzy Jaccard\n(near-duplicate samples)",
    "Layer C — expression Pearson r\n(extreme-gene fingerprint pre-screen)",
    "Layer D — series-level Jaccard\n(‘face-swapped’ re-accessioning)",
]
for y, c, t in zip(ly, cols, labels):
    box(4.0, y, 5.7, 1.25, t, c, fs=9.5)
    arrow(3.3, y+0.62, 4.0, y+0.62, color="#888888", lw=1.2)

# baseline-aware gate
box(0.3, 2.0, 3.0, 2.4, "Baseline-aware gate (L5)\n• r ≥ max(per-series p99.9, 0.97)\n• probe-space integrity filter\n(cut pseudocorrelation)",
    "#C0392B", fs=9)
arrow(3.3, 3.2, 4.0, 4.0, color="#C0392B", lw=1.6)
arrow(3.3, 3.2, 4.0, 3.0, color="#C0392B", lw=1.6)

# AI semantic review
box(0.3, 0.3, 3.0, 1.3, "AI / LLM semantic review\n(final retain / remove)", "#7F8C8D", fs=9.5)
arrow(3.3, 0.95, 4.0, 1.4, color="#7F8C8D", lw=1.4)

# output
box(4.0, 0.3, 5.7, 1.3, "Outputs: dedup_report + clean_sample_list\n(KEEP / REMOVE per sample)", "#148F77", fs=9.5)
arrow(6.7, 1.6, 6.7, 2.0, color="#888888", lw=1.2)

ax.text(5.0, 9.85, "Reproducible SOP for cross-GSE sample redundancy & ‘face-swapping’ detection",
        ha="center", fontsize=11.5, fontweight="bold", color="#2C3E50")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig1_SOP_workflow.png", bbox_inches="tight")
print("Fig1 done")

# =====================================================================
# FIG 2 — Corpus composition
# =====================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.2))
# (a) series per platform
names = list(plat.keys()) + ["Other\nplatforms"]
sv = [plat[k][0] for k in plat] + [others_series]
colors = ["#1F77B4", "#9467BD", "#D9822B", "#2CA02C", "#95A5A6"]
y = np.arange(len(names))
ax1.barh(y, sv, color=colors)
ax1.set_yticks(y); ax1.set_yticklabels(names, fontsize=10)
ax1.invert_yaxis()
for i, v in enumerate(sv):
    ax1.text(v+2, i, str(v), va="center", fontsize=9)
ax1.set_xlabel("Number of GEO series")
ax1.set_title(f"(a) Corpus by platform  (total {n_series} series)", fontsize=11)
ax1.set_xlim(0, max(sv)*1.15)

# (b) redundancy donut
removed = n_remove
kept = n_samples - removed
ax2.pie([kept, removed], labels=["Retained\n%d" % kept, "To remove\n%d" % removed],
        colors=["#BDC3C7", "#C0392B"], startangle=90, counterclock=False,
        wedgeprops=dict(width=0.42, edgecolor="white"),
        textprops=dict(fontsize=10))
ax2.text(0, 0, f"{red_rate*100:.1f}%\nredundancy", ha="center", va="center",
         fontsize=13, fontweight="bold", color="#C0392B")
ax2.set_title(f"(b) Redundancy in {n_samples:,} samples", fontsize=11)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig2_corpus_composition.png", bbox_inches="tight")
print("Fig2 done")

# =====================================================================
# FIG 3 — Baseline-aware calibration
# =====================================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.3))
# (a) total pairs
ax1.bar(["Absolute\nr≥0.97", "Baseline-aware\nretained"], [raw_edges, retained_edges],
        color=["#E67E22", "#2C3E50"])
for i, v in enumerate([raw_edges, retained_edges]):
    ax1.text(i, v+600, f"{v:,}", ha="center", fontsize=10, fontweight="bold")
ax1.set_ylabel("Cross-series duplicate sample-pairs")
ax1.set_title(f"(a) Pseudocorrelation removed: {100*(1-retained_edges/raw_edges):.1f}%",
              fontsize=10.5)
ax1.set_ylim(0, raw_edges*1.15)

# (b) representative pairs
pairs = [
    ("GSE29621↔\nGSE17536\n(pseudocorrelation)", 2825, 66),
    ("GSE161158↔\nGSE17536", 6397, 184),
    ("GSE102079↔\nGSE112790", 1109, 176),
    ("GSE21510↔\nGSE27854", 824, 104),
    ("GSE50831↔\nGSE50832", 1861, 225),
]
xs = np.arange(len(pairs))
w = 0.38
ax2.bar(xs-w/2, [p[1] for p in pairs], w, label="Absolute r≥0.97", color="#E67E22")
ax2.bar(xs+w/2, [p[2] for p in pairs], w, label="Baseline-aware", color="#2C3E50")
ax2.set_xticks(xs); ax2.set_xticklabels([p[0] for p in pairs], fontsize=7.6)
ax2.set_ylabel("Duplicate pairs")
ax2.legend(fontsize=8.5, loc="upper right")
ax2.set_title("(b) Selected series-pairs: cut vs retained", fontsize=10.5)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig3_baseline_calibration.png", bbox_inches="tight")
print("Fig3 done")

# =====================================================================
# FIG 4 — Juntendo face-swap evidence chain
# =====================================================================
fig, ax = plt.subplots(figsize=(9.4, 3.6))
ax.set_xlim(0, 10); ax.set_ylim(0, 4); ax.axis("off")
nodes = [("GSE18105", "2009", "PMID 20162577"),
         ("GSE22598", "2010", "PMID 21922135"),
         ("GSE32323", "2011", "PMID 22399497")]
xpos = [1.0, 4.5, 8.0]
for x, (g, yr, pmid) in zip(xpos, nodes):
    ax.add_patch(FancyBboxPatch((x-0.85, 1.4), 1.7, 1.4, boxstyle="round,pad=0.05",
                 fc="#2C3E50", ec="#1A252F", lw=1.6))
    ax.text(x, 2.45, g, ha="center", fontsize=12, fontweight="bold", color="white")
    ax.text(x, 2.0, yr, ha="center", fontsize=10, color="#D6DBDF")
    ax.text(x, 1.62, pmid, ha="center", fontsize=7.5, color="#AEB6BF")
arrow = dict(arrowstyle="-|>", color="#C0392B", lw=2.4, mutation_scale=16)
ax.annotate("", xy=(xpos[1]-0.85, 2.1), xytext=(xpos[0]+0.85, 2.1), arrowprops=arrow)
ax.annotate("", xy=(xpos[2]-0.85, 2.1), xytext=(xpos[1]+0.85, 2.1), arrowprops=arrow)
ax.text(2.75, 2.55, "same 34 arrays\n(r = 1.0)", ha="center", fontsize=8.5, color="#C0392B")
ax.text(6.25, 2.55, "re-accessioned\nNEW GSM → face-swap", ha="center", fontsize=8.5, color="#C0392B")
ax.text(5.0, 0.7, "Juntendo University (submitter K. Mogushi) — one cohort re-published three times;\n"
        "GSE32323 conceals replication by assigning entirely new GSM identifiers (Layer A blind, Layer C catches).",
        ha="center", fontsize=8.8, color="#34495E")
ax.set_title("Evidence chain: L2 ‘face-swap’ re-accessioning", fontsize=11.5, fontweight="bold")
fig.tight_layout()
fig.savefig(f"{OUT}/Fig4_juntendo_faceswap.png", bbox_inches="tight")
print("Fig4 done")

# =====================================================================
# FIG 5 — Positive control r distribution (GSE41258 vs GSE68468)
# =====================================================================
gsa = set(s for s, g in s2g.items() if g == "GSE41258")
gsb = set(s for s, g in s2g.items() if g == "GSE68468")
r_a = [v[0] for k, v in bm.items() if k.strip('"') in gsa and v[1] == "GSE68468"]
r_b = [v[0] for k, v in bm.items() if k.strip('"') in gsb and v[1] == "GSE41258"]
fig, ax = plt.subplots(figsize=(7.6, 4.2))
bins = np.linspace(0.6, 1.0, 33)
ax.hist(r_a, bins=bins, alpha=0.6, color="#1F77B4", label=f"GSE41258 → GSE68468 (n={len(r_a)})")
ax.hist(r_b, bins=bins, alpha=0.6, color="#C0392B", label=f"GSE68468 → GSE41258 (n={len(r_b)})")
ax.axvline(0.95, color="#2C3E50", ls="--", lw=1.3)
ax.axvline(0.981, color="#27AE60", ls=":", lw=1.6)
ax.text(0.951, ax.get_ylim()[1]*0.9, " r=0.95", color="#2C3E50", fontsize=9)
ax.text(0.84, ax.get_ylim()[1]*0.75, "peak r=0.981", color="#27AE60", fontsize=9)
ax.set_xlabel("Best cross-series Pearson r")
ax.set_ylabel("Number of samples")
ax.set_title("Positive control: GSE41258 ↔ GSE68468 (verified re-accession)", fontsize=10.5)
ax.legend(fontsize=8.5)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig5_positive_control.png", bbox_inches="tight")
print("Fig5 done; r_a=%d r_b=%d" % (len(r_a), len(r_b)))

# =====================================================================
# FIG 6 — Probe-space integrity filter
# =====================================================================
excl = d["excluded_by_probe_integrity"]
zero_gene = [e["gse"] for e in excl if e["gene_count"] == 0]
ratio_bad = [e["gse"] for e in excl if e["gene_count"] > 0]
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.6, 4.2))
# (a) GSE20916 collapse example
ax1.bar(["Full GPL570\nspace", "With GSE20916\n(probe subset)"], [54675, 27697],
        color=["#27AE60", "#E74C3C"])
for i, v in enumerate([54675, 27697]):
    ax1.text(i, v+800, f"{v:,}", ha="center", fontsize=9, fontweight="bold")
ax1.set_ylabel("Common gene features")
ax1.set_title("(a) Probe-subset collapse example\nGSE20916: 27697/54675 genes → 88% of\ncomparison space lost", fontsize=9.5)
# (b) excluded series by reason
cats = {"No expression\nmatrix (0 genes)": len(zero_gene), "Probe ratio\n< 0.9": len(ratio_bad)}
xs = np.arange(len(cats))
ax2.bar(xs, list(cats.values()), color=["#8E44AD", "#E67E22"])
for i, v in enumerate(cats.values()):
    ax2.text(i, v+0.2, str(v), ha="center", fontsize=10, fontweight="bold")
ax2.set_xticks(xs); ax2.set_xticklabels(list(cats.keys()), fontsize=9)
ax2.set_ylabel("Series excluded")
ax2.set_title(f"(b) Probe-integrity exclusions (n={len(excl)})\nprevents space collapse & false inflation",
              fontsize=9.5)
fig.tight_layout()
fig.savefig(f"{OUT}/Fig6_probe_integrity.png", bbox_inches="tight")
print("Fig6 done; excl=%d" % len(excl))

print("ALL FIGURES WRITTEN TO", OUT)
