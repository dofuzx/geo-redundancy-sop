#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""实验B·核心：RNA-seq 版 GEO 跨系列冗余检测（四层框架的平台迁移版）。

相对 analyze_tier.py (microarray) 的适配:
  * Layer C: Pearson(按样本中心化) → Spearman（先按列求秩, 秩上做 Pearson;
    Spearman 对任意单调变换不变, 故 log/归一化差异不影响结果）。
  * 基因空间完整性 → 命名空间感知的完整性策略:
      - 基因 ID 归一化（ENSG 去版本号、符号大写）;
      - 按 ID 命名空间分三组: ENSG / SYMBOL / OTHER;
      - 组内构建"多数基因空间"(出现在 >= --space-frac 系列中的基因);
      - 系列须覆盖组空间 >= --gene-integrity, 否则判为残缺剔除;
      - 组空间内缺失基因用该样本观测值列均值插补;
      - Layer C 仅在同名空间组内比较（跨命名空间不做 ID 映射, 如实报告）。
  * 新增计数数据 QC: 每样本基因检出率 >= --detect-ratio x 组中位数;
    基因需在 >= --min-gene-frac 组内样本中检出。
  * 基线感知判据完全复用: 交叉 rho 须 >= max(双方组内 p99.9 基线, 绝对下限)。
  * Layer A (GSM 重叠) / Layer D (标题元数据 Jaccard) / 连通分量 / 报告 schema 同构。

注意: 不同系列的样本表头可能重名（TCGA 条码、S1/S2 等）, 故样本全局键为 (GSE, 样本名)。
"""
import argparse, csv, json, os, re
from collections import defaultdict, Counter
import numpy as np
from scipy.stats import rankdata

R_THRESH = 0.97
BASE_Q = 99.9
BASE_FLOOR = 0.97
TOKEN = re.compile(r"[\s,.;:()\[\]/\"'<>-]+")
ENSG_RE = re.compile(r"^ENSG\d+", re.I)


def tokenize(text):
    if not text:
        return set()
    return set(t for t in TOKEN.split(str(text).lower()) if len(t) > 3)


def read_meta(path):
    meta = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            sid = (row.get("sample_id") or "").strip()
            if sid:
                meta[sid] = row
    return meta


def norm_gid(x):
    x = str(x).strip().strip('"').strip("'")
    if ENSG_RE.match(x):
        return re.split(r"\.", x)[0].upper()
    return x.upper()


def series_gene_ids(path):
    """返回 (ids 归一化列表, 命名空间)。"""
    ensg = sym = other = 0
    ids = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        r = csv.reader(f)
        next(r, None)
        for row in r:
            if not row or not row[0].strip():
                continue
            raw = row[0].strip()
            n = norm_gid(raw)
            ids.append(n)
            if n.startswith("ENSG"):
                ensg += 1
            elif re.fullmatch(r"[A-Z][A-Z0-9.-]{0,14}", n):
                sym += 1
            else:
                other += 1
    tot = ensg + sym + other or 1
    ns = "ENSG" if ensg > 0.5 * tot else ("SYMBOL" if sym > 0.5 * tot else "OTHER")
    return ids, ns


def classify_group(glist, ids_by_g, space_frac):
    """按命名空间构建多数基因空间 + 覆盖率过滤。返回 (keep, space, dropped)。"""
    cnt = Counter()
    for g in glist:
        cnt.update(set(ids_by_g[g]))
    thr = space_frac * len(glist)
    space = sorted(g for g, c in cnt.items() if c >= thr)
    keep, dropped = [], []
    for g in glist:
        s = set(ids_by_g[g])
        cov = len(s & set(space)) / len(space) if space else 0.0
        if cov >= GENE_INTEGRITY:
            keep.append(g)
        else:
            dropped.append((g, len(s), f"组空间覆盖率 {cov:.2f} < {GENE_INTEGRITY}"))
    return keep, space, dropped


def load_group_matrix(glist, sopdir, space):
    """加载组内多数基因空间矩阵 (G×N, float32), 缺失=NaN(后插补)。"""
    g2i = {g: i for i, g in enumerate(space)}
    sample_list, series_at = [], []
    col_off = {}
    for g in glist:
        with open(os.path.join(sopdir, g + ".csv"), newline="", encoding="utf-8-sig") as f:
            r = csv.reader(f)
            h = next(r)
            cols = [x.strip() for x in h[1:]]
            col_off[g] = len(sample_list)
            sample_list.extend(cols)
            series_at.extend([g] * len(cols))
    M = np.full((len(space), len(sample_list)), np.nan, dtype=np.float32)
    for g in glist:
        base = col_off[g]
        with open(os.path.join(sopdir, g + ".csv"), newline="", encoding="utf-8-sig") as f:
            r = csv.reader(f)
            h = next(r)
            ncol = len(h)
            for row in r:
                if len(row) < 2:
                    continue
                gi = g2i.get(norm_gid(row[0]))
                if gi is None:
                    continue
                for k, v in enumerate(row[1:ncol]):
                    if v == "":
                        continue
                    try:
                        M[gi, base + k] = float(v)
                    except ValueError:
                        pass
    return M, sample_list, series_at


def impute_colmean(M):
    """NaN → 该样本(列)在组空间内的观测均值。"""
    with np.errstate(invalid="ignore"):
        cm = np.nanmean(M, axis=0, keepdims=True)
    cm = np.where(np.isfinite(cm), cm, 0.0)
    idx = ~np.isfinite(M)
    M[idx] = np.broadcast_to(cm, M.shape)[idx]
    return M


def run_group(ns_name, glist, sopdir, space, args, acc, series):
    """跑一个命名空间组, 结果累加进 acc。"""
    if len(glist) < 2 or not space:
        print(f"[{ns_name}] 系列 {len(glist)}, 空间 {len(space)} — 组内不可比较, 跳过", flush=True)
        return
    M, sample_list, series_at = load_group_matrix(glist, sopdir, space)
    M = impute_colmean(M)
    # ---- 列级防御: 近常数列(注释列残留, 如染色体索引)剔除 ----
    Ms = np.sort(M, axis=0)
    dup_frac = (Ms[1:] == Ms[:-1]).mean(axis=0) if M.shape[0] > 1 else np.zeros(M.shape[1])
    del Ms
    bad = dup_frac > 0.98
    n_bad = int(bad.sum())
    if n_bad:
        print(f"[{ns_name} 防御] 剔除 {n_bad} 个近常数列(疑似注释列残留)", flush=True)
        acc["const_cols_dropped"] = acc.get("const_cols_dropped", 0) + n_bad
        M = M[:, ~bad]
        sample_list = [s for s, b in zip(sample_list, bad) if not b]
        series_at = [g for g, b in zip(series_at, bad) if not b]
    sa = np.array(series_at)
    print(f"[{ns_name}] 矩阵 {M.shape[0]} 基因 x {M.shape[1]} 样本 ({len(glist)} 系列)", flush=True)

    # ---- 样本级 QC: 基因检出率 ----
    detected = np.count_nonzero(M > 0, axis=0)
    med = np.median(detected) or 1
    drop_samples = detected < args.detect_ratio * med
    n_drop = int(drop_samples.sum())
    if n_drop:
        print(f"[{ns_name} QC] 检出率过滤: 剔除 {n_drop} 样本 (<{args.detect_ratio}x中位数={int(args.detect_ratio*med)})", flush=True)
    for j in np.where(drop_samples)[0]:
        acc["qc_dropped"].append([sample_list[j], series_at[j], int(detected[j])])
    M = M[:, ~drop_samples]
    sample_list = [s for s, d in zip(sample_list, drop_samples) if not d]
    series_at = [g for g, d in zip(series_at, drop_samples) if not d]
    sa = np.array(series_at)
    n = len(sample_list)
    acc["n_samples"] += n
    acc["n_series"] += len(glist)
    acc["group_detail"].append({
        "namespace": ns_name, "n_series": len(glist), "n_samples": n,
        "n_genes_space": len(space), "qc_samples_dropped": n_drop,
    })

    # ---- 基因级 QC: 检出频率 ----
    gene_det = np.count_nonzero(M > 0, axis=1)
    keep_genes = gene_det >= args.min_gene_frac * max(1, M.shape[1])
    M = M[keep_genes, :]
    genes_qc = int(keep_genes.sum())
    acc["n_genes_qc"] += genes_qc
    print(f"[{ns_name} QC] 基因检出过滤: {genes_qc}/{len(space)} 基因进入比对空间", flush=True)

    # ---- Spearman: 列秩 → 中心化/单位化 ----
    Mr = rankdata(M, axis=0, method="average").astype(np.float32)
    del M
    Mr -= Mr.mean(axis=0, keepdims=True)
    nrm = np.sqrt((Mr ** 2).sum(axis=0, keepdims=True))
    Mr /= np.where(nrm == 0, 1, nrm)
    del nrm

    cols_of = {g: np.where(sa == g)[0] for g in glist}
    print(f"[{ns_name} 扫描] {len(glist)} 系列 / {n} 样本 / {genes_qc} 基因", flush=True)

    # ---- 组内基线 ----
    baseline = {}
    for g in glist:
        idx = cols_of[g]
        if len(idx) < 4:
            baseline[g] = 0.0
            continue
        Cw = np.clip(Mr[:, idx].T @ Mr[:, idx], -1, 1)
        iu = np.triu_indices(len(idx), k=1)
        v = Cw[iu]
        baseline[g] = float(np.percentile(v, BASE_Q)) if len(v) else 0.0
        del Cw, v
    acc["baseline"].update(baseline)

    # ---- 分块扫描（仅组内交叉）----
    g_abs = g_ada = 0
    for ai, g1 in enumerate(glist):
        idx1 = cols_of[g1]
        if len(idx1) == 0:
            continue
        block = np.clip(Mr[:, idx1].T @ Mr, -1, 1)     # n1 × N(组内)
        for bi in range(len(idx1)):
            i = idx1[bi]
            vals = block[bi][sa != g1]
            if not vals.size:
                continue
            jb = int(np.argmax(vals))
            colsb = np.where(sa != g1)[0]
            acc["bestmatch"][f"{g1}|{sample_list[i]}"] = [
                round(float(vals[jb]), 4), series_at[colsb[jb]], sample_list[colsb[jb]]]
        for bj in range(ai + 1, len(glist)):
            g2 = glist[bj]
            idx2 = cols_of[g2]
            v = block[:, idx2].ravel()
            n_abs = int(np.count_nonzero(v >= R_THRESH))
            T = R_THRESH
            if args.baseline_mode != "off":
                T = max(T, baseline.get(g1, 0.0), baseline.get(g2, 0.0), BASE_FLOOR)
            sel = np.where(v >= T)[0]
            n_ada = int(len(sel))
            g_abs += n_abs
            g_ada += n_ada
            if n_ada or n_abs:
                acc["pair_stats"].append({
                    "g1": g1, "g2": g2, "namespace": ns_name,
                    "n_cross_pairs": int(v.size),
                    f"baseline_p{BASE_Q:g}_g1": round(baseline.get(g1, 0.0), 4),
                    f"baseline_p{BASE_Q:g}_g2": round(baseline.get(g2, 0.0), 4),
                    "threshold_used": round(float(T), 4),
                    "edges_abs_0.97": n_abs, "edges_baseline_aware": n_ada,
                    "cross_max": round(float(v.max()), 4) if v.size else None,
                    "n_r>=0.9999": int(np.count_nonzero(v >= 0.9999)),
                    "reduced_by": round(1 - n_ada / n_abs, 4) if n_abs else None,
                })
            for kk in sel:
                i = idx1[kk // len(idx2)]
                j = idx2[kk % len(idx2)]
                acc["edges"].append({"g1": g1, "s1": sample_list[i], "g2": g2,
                                     "s2": sample_list[j],
                                     "spearman": round(float(v[kk]), 4)})
        if (ai + 1) % 50 == 0 or ai + 1 == len(glist):
            print(f"  [{ns_name}] 扫描 {ai+1}/{len(glist)}  绝对边={g_abs}  感知边={g_ada}", flush=True)
        del block
    acc["abs_edges"] += g_abs
    acc["ada_edges"] += g_ada


def main():
    ap = argparse.ArgumentParser(description="RNA-seq 版 GEO 跨系列冗余检测（实验B）")
    ap.add_argument("--sopdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", default="CRC-RNAseq")
    ap.add_argument("--rthresh", type=float, default=0.97)
    ap.add_argument("--baseline-mode", choices=["off", "adapt"], default="adapt")
    ap.add_argument("--baseline-q", type=float, default=99.9)
    ap.add_argument("--baseline-floor", type=float, default=0.97)
    ap.add_argument("--gene-integrity", type=float, default=0.6,
                    help="系列对组多数基因空间的覆盖率须 >= 该值")
    ap.add_argument("--space-frac", type=float, default=0.5,
                    help="基因出现在组内 >= 该比例系列中才进入组空间")
    ap.add_argument("--detect-ratio", type=float, default=0.9,
                    help="每样本基因检出数须 >= 该比例 x 组中位数")
    ap.add_argument("--min-gene-frac", type=float, default=0.05,
                    help="基因至少在 该比例 组内样本中检出才进入比对空间")
    ap.add_argument("--enum-tsv", default=r"D:/AIwork/GEO_RNAseq/rnaseq_enum.tsv",
                    help="系列标题元数据(Layer D 用)")
    args = ap.parse_args()
    global R_THRESH, BASE_Q, BASE_FLOOR, GENE_INTEGRITY
    R_THRESH, BASE_Q, BASE_FLOOR = args.rthresh, args.baseline_q, args.baseline_floor
    GENE_INTEGRITY = args.gene_integrity
    os.makedirs(args.out, exist_ok=True)

    # ---------- 系列清单 ----------
    series = {}
    for fn in sorted(os.listdir(args.sopdir)):
        if fn.endswith(".csv"):
            g = fn[:-4]
            if os.path.exists(os.path.join(args.sopdir, g + "_meta.tsv")):
                series[g] = read_meta(os.path.join(args.sopdir, g + "_meta.tsv"))
    gses = sorted(series.keys())
    print(f"[载入] {len(gses)} 个 RNA-seq 系列", flush=True)

    # ---------- 归一化 ID + 命名空间分组 ----------
    ids_by_g, ns_of = {}, {}
    dropped_empty = []
    for g in gses:
        ids, ns = series_gene_ids(os.path.join(args.sopdir, g + ".csv"))
        if not ids:
            dropped_empty.append((g, 0, "空矩阵"))
            continue
        ids_by_g[g] = ids
        ns_of[g] = ns
    print(f"[命名空间] ENSG={sum(1 for v in ns_of.values() if v=='ENSG')} "
          f"SYMBOL={sum(1 for v in ns_of.values() if v=='SYMBOL')} "
          f"OTHER={sum(1 for v in ns_of.values() if v=='OTHER')}", flush=True)

    acc = {"edges": [], "pair_stats": [], "bestmatch": {}, "baseline": {},
           "qc_dropped": [], "n_samples": 0, "n_series": 0,
           "n_genes_qc": 0, "abs_edges": 0, "ada_edges": 0, "group_detail": []}
    excluded = [{"gse": g, "gene_count": c, "reason": why} for g, c, why in dropped_empty]
    for ns_name in ["ENSG", "SYMBOL", "OTHER"]:
        glist = [g for g in gses if ns_of.get(g) == ns_name]
        keep, space, dropped = classify_group(glist, ids_by_g, args.space_frac)
        for g, c, why in dropped:
            excluded.append({"gse": g, "gene_count": c, "reason": why, "namespace": ns_name})
        print(f"[{ns_name}] 保留 {len(keep)}/{len(glist)}, 组空间 {len(space)} 基因, 剔除 {len(dropped)}", flush=True)
        acc["group_detail"].append({
            "namespace": ns_name, "n_series_input": len(glist),
            "n_series_kept": len(keep), "n_genes_space_raw": len(space),
            "n_excluded": len(dropped)})
        run_group(ns_name, keep, args.sopdir, space, args, acc, series)

    keep_all = sorted({e[k] for e in acc["edges"] for k in ("g1", "g2")} |
                      set(acc["baseline"].keys()))
    c_edges, pair_stats = acc["edges"], acc["pair_stats"]
    n = acc["n_samples"]

    # ---------- Layer A ----------
    gsm_map = defaultdict(list)
    for g in gses:
        for s in series[g]:
            gsm_map[s].append(g)
    layerA = {g: [{"gsm": s, "also_in": [x for x in sorted(set(v)) if x != g]}]
              for s, v in gsm_map.items() if len(set(v)) > 1 for g in set(v)}

    # ---------- Layer D（标题元数据 Jaccard）----------
    titles = {}
    try:
        with open(args.enum_tsv, encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                titles[row["gse"]] = row.get("title") or ""
    except OSError:
        pass
    toks = {g: tokenize(titles.get(g, "")) for g in gses}
    layerD = []
    for i in range(len(gses)):
        for j in range(i + 1, len(gses)):
            g1, g2 = gses[i], gses[j]
            tt1, tt2 = toks[g1], toks[g2]
            jc = len(tt1 & tt2) / len(tt1 | tt2) if (tt1 | tt2) else 0
            if jc >= 0.7:
                shared = set(series[g1]) & set(series[g2])
                layerD.append({"g1": g1, "g2": g2, "meta_jaccard": round(jc, 3),
                               "sample_name_overlap": sorted(shared)[:20],
                               "verdict": "FACE_SWAP_CANDIDATE" if not shared else "OVERLAP_BUT_SHARED_NAME"})

    # ---------- 连通分量 / 冗余率 ----------
    # 枚举所有保留系列的样本表头（扣除 QC 剔除样本），构建全局 (gse, sample) 并查集
    qc_dropped_set = {(g, s) for s, g, _ in acc["qc_dropped"]}
    all_samples = set()
    for g in keep_all:
        with open(os.path.join(args.sopdir, g + ".csv"), newline="", encoding="utf-8-sig") as f:
            r = csv.reader(f)
            h = next(r)
            for s in h[1:]:
                if (g, s.strip()) not in qc_dropped_set:
                    all_samples.add((g, s.strip()))
    parent = {x: x for x in all_samples}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for e in c_edges:
        a, b = (e["g1"], e["s1"]), (e["g2"], e["s2"])
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    gd = defaultdict(list)
    for x in parent:
        gd[find(x)].append(x)
    dup_groups = [sorted(f"{g}|{s}" for g, s in v) for v in gd.values() if len(v) > 1]
    remove_set = set()
    for grp in dup_groups:
        for x in grp[1:]:
            remove_set.add(x)
    total_samples = len(all_samples)
    redundancy_rate = len(remove_set) / total_samples if total_samples else 0.0

    pair_stats.sort(key=lambda x: -x["edges_baseline_aware"])
    base_vals = [v for v in acc["baseline"].values() if v > 0]
    report = {
        "summary": {
            "n_series_input": len(gses), "n_series_used": acc["n_series"],
            "n_samples": total_samples,
            "n_genes_qc": acc["n_genes_qc"],
            "qc_samples_dropped": len(acc["qc_dropped"]), "qc_sample_detail": acc["qc_dropped"][:200],
            "duplicate_sample_pairs": len(c_edges), "duplicate_groups": len(dup_groups),
            "samples_to_remove": len(remove_set), "redundancy_rate": round(redundancy_rate, 4),
            "r_threshold": R_THRESH, "baseline_mode": args.baseline_mode,
            "baseline_q": BASE_Q, "baseline_floor": BASE_FLOOR,
            "gene_integrity": args.gene_integrity, "space_frac": args.space_frac,
            "detect_ratio": args.detect_ratio, "min_gene_frac": args.min_gene_frac,
            "edges_absolute_0.97_total": acc["abs_edges"],
            "edges_baseline_aware_total": acc["ada_edges"],
            "pseudo_correlations_removed": acc["abs_edges"] - acc["ada_edges"],
            "near_identity_edges_r>=0.9999": sum(1 for e in c_edges if e["spearman"] >= 0.9999),
            "median_intra_series_baseline": round(float(np.median(base_vals)), 4) if base_vals else 0.0,
            "namespace_groups": acc["group_detail"],
            "n_excluded_by_integrity": len(excluded),
        },
        "excluded_by_gene_integrity": excluded,
        "layerA_gsm_overlap": layerA,
        "layerC_expression_duplicates": c_edges,
        "layerC_pair_stats": pair_stats,
        "layerD_faceswap": layerD,
        "duplicate_groups": dup_groups,
        "best_match_cross": acc["bestmatch"],
        "recommended_remove": sorted(remove_set),
        "series_baselines": {g: round(v, 4) for g, v in acc["baseline"].items()},
    }
    with open(os.path.join(args.out, "dedup_report_rnaseq.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    # ---------- Markdown 报告 ----------
    L = [f"# RNA-seq 扩展检测报告 ({args.name})", ""]
    L.append(f"- 输入系列 {len(gses)} → 命名空间分组后可用 **{acc['n_series']}**；样本 **{total_samples}**")
    L.append(f"- 相关度量: **Spearman**（列秩 + Pearson，对单调变换不变）；比较仅在同名空间组内进行")
    L.append(f"- 判据: 基线感知 T = max(双方组内p{BASE_Q:g}, {BASE_FLOOR})；组内基线中位数 = {report['summary']['median_intra_series_baseline']}")
    L.append(f"- 绝对阈值 rho>={R_THRESH} 原始边: **{acc['abs_edges']}**")
    L.append(f"- 基线感知后重复边: **{acc['ada_edges']}**（滤除伪相关 {acc['abs_edges'] - acc['ada_edges']} 对）")
    L.append(f"- 重复组 **{len(dup_groups)}**；建议剔除 **{len(remove_set)}**；**冗余率 {redundancy_rate*100:.2f}%**")
    L.append(f"- 近乎完全一致 (rho>=0.9999): **{sum(1 for e in c_edges if e['spearman'] >= 0.9999)}** 对")
    L.append("")
    L.append("## 命名空间分组")
    for gd_ in acc["group_detail"]:
        L.append(f"- **{gd_['namespace']}**: 输入 {gd_.get('n_series_input', gd_.get('n_series','-'))} 系列 → "
                 f"保留 {gd_.get('n_series_kept', gd_.get('n_series','-'))}, "
                 f"组空间 {gd_.get('n_genes_space_raw','-')} 基因, 比对样本 {gd_.get('n_samples', 0)}, "
                 f"QC 剔除样本 {gd_.get('qc_samples_dropped', 0)}")
    L.append("")
    if excluded:
        L.append(f"## 基因空间完整性过滤（剔除 {len(excluded)}）")
        for x in excluded[:40]:
            L.append(f"- 剔除 `{x['gse']}`（{x['gene_count']} 基因）：{x['reason']}")
        if len(excluded) > 40:
            L.append(f"- ... 共 {len(excluded)} 个")
        L.append("")
    if pair_stats:
        L.append("## TOP 系列对（基线感知 vs 绝对阈值）")
        L.append("")
        L.append("| 系列A | 系列B | 空间 | 交叉对数 | 基线A | 基线B | 阈值 | 绝对边 | 感知边 | max rho |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for p in pair_stats[:20]:
            L.append(f"| {p['g1']} | {p['g2']} | {p['namespace']} | {p['n_cross_pairs']} | "
                     f"{p['baseline_p%.4g_g1' % BASE_Q]} | {p['baseline_p%.4g_g2' % BASE_Q]} | "
                     f"{p['threshold_used']} | {p['edges_abs_0.97']} | {p['edges_baseline_aware']} | {p['cross_max']} |")
        L.append("")
    L.append("## Layer C — 表达层重复（Spearman, TOP 40）")
    if c_edges:
        for e in sorted(c_edges, key=lambda x: -x["spearman"])[:40]:
            L.append(f"- `{e['g1']}/{e['s1']}` ≈ `{e['g2']}/{e['s2']}` (rho={e['spearman']})")
        if len(c_edges) > 40:
            L.append(f"- ... 共 {len(c_edges)} 对")
    else:
        L.append("- 未检出跨 GSE 表达重复。")
    L.append("")
    if layerD:
        L.append("## Layer D — 变脸候选（标题 Jaccard>=0.7）")
        for d in layerD:
            L.append(f"- `{d['g1']}` ↔ `{d['g2']}` J={d['meta_jaccard']} → {d['verdict']}")
        L.append("")
    with open(os.path.join(args.out, "dedup_report_rnaseq.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    with open(os.path.join(args.out, "clean_sample_list_rnaseq.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gse", "sample", "action"])
        for g, s in sorted(all_samples):
            w.writerow([g, s, "REMOVE" if f"{g}|{s}" in remove_set else "KEEP"])
    print(f"[完成] 系列 {acc['n_series']} / 样本 {total_samples} / 重复对 {len(c_edges)} / 冗余率 {redundancy_rate*100:.2f}%", flush=True)
    print(f"[完成] 报告 -> {args.out}", flush=True)


if __name__ == "__main__":
    main()
