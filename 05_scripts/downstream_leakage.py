#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Experiment A (downstream impact): redundant cross-GSE samples inflate
classifier AUC via data leakage.

Builds a tumor-vs-normal LogisticRegression on the GPL570 subset of the
CRC corpus, then compares four evaluation conditions:
  RAW  + random split      (all samples, leakage allowed)
  RAW  + grouped split     (by GSE, within-series leakage removed)
  CLEAN + random split     (REMOVE samples excluded, cross-GSE redundancy gone)
  CLEAN + grouped split    (cleanest baseline)

Also estimates the train/test leakage rate of redundant (REMOVE) samples.
Writes downstream_leakage_results.json and Fig7_leakage_auc.png.
"""
import os, glob, csv, json, re
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold, GroupShuffleSplit
from collections import defaultdict

SOP_DIRS = [
    "D:/AIwork/GEO_Tire3_sop",
    "D:/AIwork/GEO_Tire4_sop",
    "E:/CRC_GEO_2026T1_sop_v2",
    "E:/CRC_GEO_2026T2_sop_v2",
]
CLEAN_CSV = "D:/AIwork/GEO_T1234_final/clean_sample_list.csv"
DEDUP_JSON = "D:/AIwork/GEO_T1234_final/dedup_report.json"
K = 1000
OUT_JSON = "D:/AIwork/GEO/paper/07_downstream/downstream_leakage_results.json"
FIG = "D:/AIwork/GEO/paper/07_downstream/Fig7_leakage_auc.png"
FIG_FIGS = "D:/AIwork/GEO/paper/01_manuscript/manuscript_figs/Fig7_leakage_auc.png"

def norm(s):
    return s.strip().strip('"').strip()

# ---------- action map (sample -> KEEP/REMOVE) ----------
action = {}
with open(CLEAN_CSV, encoding='utf-8') as f:
    r = csv.reader(f)
    next(r)
    for row in r:
        if len(row) < 3:
            continue
        action[norm(row[1])] = row[2].strip()
print("action map size:", len(action))

# ---------- label heuristic ----------
normal_re = re.compile(r'normal|healthy|control|adjacent|mucosa|non-?tumor|benign|polyp|uninvolved|nonmalign|non-malignant|noncancer|non-cancer')
tumor_re  = re.compile(r'tumor|tumour|carcinoma|adenocarcinoma|cancer|neoplas|malignant|lesion|colorectal cancer|tumoral')
def label_of(txt):
    n = bool(normal_re.search(txt)); t = bool(tumor_re.search(txt))
    if n and not t: return 0
    if t and not n: return 1
    return None

def meta_labels(metap):
    out = {}
    with open(metap, encoding='utf-8', errors='replace') as f:
        rd = csv.reader(f, delimiter='\t')
        header = next(rd)
        i_gsm = header.index('gsm')
        i_src = header.index('Sample_source_name_ch1') if 'Sample_source_name_ch1' in header else None
        i_ch  = header.index('Sample_characteristics_ch1') if 'Sample_characteristics_ch1' in header else None
        i_ti  = header.index('Sample_title') if 'Sample_title' in header else None
        for row in rd:
            if len(row) <= max(i_gsm, i_src or 0, i_ch or 0, i_ti or 0):
                continue
            sid = norm(row[i_gsm])
            blob = " ".join(row[i] for i in (i_src, i_ch, i_ti) if i is not None and i < len(row)).lower()
            out[sid] = label_of(blob)
    return out

def platform_of(metap):
    with open(metap, encoding='utf-8', errors='replace') as f:
        rd = csv.reader(f, delimiter='\t')
        header = next(rd)
        if 'platform' not in header:
            return None
        row = next(rd)
        return norm(row[header.index('platform')])

def load_matrix(csvp):
    try:
        with open(csvp, encoding='utf-8', errors='replace') as f:
            header = f.readline().rstrip('\n').split(',')
        sids = [norm(h) for h in header[1:]]
        arr = np.loadtxt(csvp, delimiter=',', skiprows=1, dtype=np.float32,
                         usecols=range(1, len(header)))          # genes x nsamples
        if arr.size == 0:
            return None, None, None
        gids = np.loadtxt(csvp, delimiter=',', skiprows=1, dtype=str, usecols=0)
        gids = np.array([x.strip('"').strip() for x in gids.astype(str)])
    except Exception as e:
        print("  load fail:", os.path.basename(csvp), e)
        return None, None, None
    if gids.size == 0:
        return None, None, None
    return sids, gids, arr

# ---------- collect GPL570 (GPL570) series ----------
gpl570 = []
for d in SOP_DIRS:
    for csvp in glob.glob(os.path.join(d, "GSE*.csv")):
        metap = csvp[:-4] + "_meta.tsv"
        if os.path.exists(metap) and platform_of(metap) == "GPL570":
            gpl570.append((csvp, metap))
print("GPL570 series:", len(gpl570))

# ---------- build MAJORITY gene space across GPL570 (>=95% of series) ----------
def load_gids(csvp):
    try:
        g = np.loadtxt(csvp, delimiter=',', skiprows=1, dtype=str, usecols=0)
        return np.array([x.strip('"').strip() for x in g.astype(str)])
    except Exception:
        return np.array([])

counts = defaultdict(int)
n_series = 0
for ci, (csvp, metap) in enumerate(gpl570):
    g = load_gids(csvp)
    if g is None or g.size == 0:
        continue
    n_series += 1
    for x in set(g.tolist()):
        counts[x] += 1
thr = 0.95 * n_series
master = np.array(sorted([gid for gid, c in counts.items() if c >= thr]))
gene_pos = {g: i for i, g in enumerate(master)}
G = len(master)
mean = np.zeros(G, np.float64)
M2 = np.zeros(G, np.float64)
count = 0
print("GPL570 series with readable gene ids:", n_series, " majority gene space (>=95%):", G)

def load_aligned(csvp):
    sids, gids, arr = load_matrix(csvp)
    if arr is None or gids is None:
        return None, None
    local = {g: i for i, g in enumerate(gids)}
    ns = arr.shape[1]
    out = np.empty((G, ns), dtype=np.float32)
    col_mean = arr.mean(axis=0)
    for mi, g in enumerate(master):
        li = local.get(g, -1)
        out[mi] = arr[li] if li >= 0 else col_mean
    return sids, out

# ---------- pass1 variance (Welford) ----------
for ci, (csvp, metap) in enumerate(gpl570):
    sids, arr = load_aligned(csvp)
    if arr is None:
        print("  skip (gene set mismatch):", csvp)
        continue
    ns = arr.shape[1]
    for j in range(ns):
        x = arr[:, j].astype(np.float64)
        count += 1
        delta = x - mean
        mean += delta / count
        M2 += delta * (x - mean)
    del arr
    if ci % 40 == 0:
        print("  pass1", ci, "samples seen", count)
var = M2 / (count - 1)
top_idx = np.argsort(var)[::-1][:K]
print("variance pass done; total GPL570 samples:", count, " top-K:", K)

# ---------- pass2: extract features + labels ----------
X_list, Y, GSEl, ACT, SID = [], [], [], [], []
for ci, (csvp, metap) in enumerate(gpl570):
    gse = os.path.basename(csvp)[:-4]
    sids, arr = load_aligned(csvp)
    if arr is None:
        continue
    feat = arr[top_idx, :].T.astype(np.float32)        # ns x K
    labels = meta_labels(metap)
    for j, sm in enumerate(sids):
        lab = labels.get(sm)
        if lab is None:
            continue
        X_list.append(feat[j])
        Y.append(lab)
        GSEl.append(gse)
        ACT.append(action.get(sm, 'KEEP'))
        SID.append(sm)
    del arr
    if ci % 40 == 0:
        print("  pass2", ci, "labeled so far", len(Y))
X = np.array(X_list, dtype=np.float32)
Y = np.array(Y, dtype=np.int64)
GSEl = np.array(GSEl)
ACT = np.array(ACT)
SID = np.array(SID)
print("labeled samples (RAW):", len(Y), " classes:", np.bincount(Y))
print("series used:", len(np.unique(GSEl)))

# ---------- cross-validation ----------
pipe = make_pipeline(StandardScaler(),
                     LogisticRegression(max_iter=8000, C=1.0))
def score(Xs, ys, groups, cv):
    return cross_val_score(pipe, Xs, ys, groups=groups, cv=cv, scoring='roc_auc')

mask_clean = ACT == 'KEEP'
auc_raw_rand   = score(X, Y, None, StratifiedKFold(10, shuffle=True, random_state=42))
auc_raw_group  = score(X, Y, GSEl, GroupShuffleSplit(n_splits=10, random_state=42))
auc_clean_rand = score(X[mask_clean], Y[mask_clean], None, StratifiedKFold(10, shuffle=True, random_state=42))
auc_clean_group= score(X[mask_clean], Y[mask_clean], GSEl[mask_clean], GroupShuffleSplit(n_splits=10, random_state=42))

# ---------- leakage rate simulation ----------
dedup = json.load(open(DEDUP_JSON, encoding='utf-8'))
labeled_set = set(SID.tolist())
sample_to_group = {}
group_indices = defaultdict(set)
gi = 0
for grp in dedup['duplicate_groups']:
    members = [norm(m[1]) for m in grp]
    members_in = [m for m in members if m in labeled_set]
    if len(members_in) >= 2:
        for m in members_in:
            sample_to_group[m] = gi
        gi += 1
for i, s in enumerate(SID.tolist()):
    if s in sample_to_group:
        group_indices[sample_to_group[s]].add(i)

redundant_idx = [i for i in range(len(SID)) if ACT[i] == 'REMOVE' and SID[i] in sample_to_group]
leak_hits = 0
leak_tot = 0
for rs in range(50):
    skf = StratifiedKFold(10, shuffle=True, random_state=rs)
    for tr, te in skf.split(X, Y):
        tr_set = set(tr.tolist())
        for i in redundant_idx:
            if i in te:
                gs = group_indices[sample_to_group[SID[i]]]
                if gs & tr_set:
                    leak_hits += 1
                leak_tot += 1

leak_rate = leak_hits / leak_tot if leak_tot else 0.0

res = {
    "platform_subset": "GPL570",
    "n_labeled_raw": int(len(Y)),
    "n_labeled_clean": int(mask_clean.sum()),
    "n_series_used": int(len(np.unique(GSEl))),
    "class_counts_raw": {str(k): int(v) for k, v in enumerate(np.bincount(Y))},
    "top_genes_K": K,
    "AUC_RAW_random_mean": float(auc_raw_rand.mean()),
    "AUC_RAW_random_std": float(auc_raw_rand.std()),
    "AUC_RAW_group_mean": float(auc_raw_group.mean()),
    "AUC_RAW_group_std": float(auc_raw_group.std()),
    "AUC_CLEAN_random_mean": float(auc_clean_rand.mean()),
    "AUC_CLEAN_random_std": float(auc_clean_rand.std()),
    "AUC_CLEAN_group_mean": float(auc_clean_group.mean()),
    "AUC_CLEAN_group_std": float(auc_clean_group.std()),
    "delta_AUC_raw_random_minus_clean_group": float(auc_raw_rand.mean() - auc_clean_group.mean()),
    "n_redundant_labeled": int(len(redundant_idx)),
    "leakage_rate_redundant_samples": float(leak_rate),
}
print(json.dumps(res, indent=2, ensure_ascii=False))
os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
json.dump(res, open(OUT_JSON, 'w', encoding='utf-8'), indent=2, ensure_ascii=False)
print("saved:", OUT_JSON)

# ---------- figure ----------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
labels = ['RAW\nrandom', 'RAW\nseries-grouped', 'CLEAN\nrandom', 'CLEAN\nseries-grouped']
means = [auc_raw_rand.mean(), auc_raw_group.mean(), auc_clean_rand.mean(), auc_clean_group.mean()]
stds  = [auc_raw_rand.std(), auc_raw_group.std(), auc_clean_rand.std(), auc_clean_group.std()]
colors = ['#c0392b', '#e67e22', '#27ae60', '#16a085']
fig, ax = plt.subplots(figsize=(7.2, 4.6))
bars = ax.bar(labels, means, yerr=stds, capsize=6, color=colors, edgecolor='black', linewidth=0.6)
ax.set_ylim(min(means) - 0.05, 1.005)
ax.set_ylabel('Cross-validated AUC (tumor vs normal)')
ax.set_title('Experiment A — redundant samples inflate classifier AUC (GPL570, K=%d)' % K)
for b, m in zip(bars, means):
    ax.text(b.get_x() + b.get_width()/2, m + 0.008, '%.3f' % m, ha='center', fontsize=9)
ax.axhline(means[3], ls='--', color='gray', lw=1)
ax.text(3.4, means[3] - 0.004, 'cleanest', color='gray', fontsize=8, ha='right')
fig.tight_layout()
os.makedirs(os.path.dirname(FIG), exist_ok=True)
fig.savefig(FIG, dpi=300)
fig.savefig(FIG_FIGS, dpi=300)
print("saved figure:", FIG)
