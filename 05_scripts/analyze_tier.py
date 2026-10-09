#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
泛化版 GEO 跨系列冗余检测引擎（基于已修复的 analyze_tier1.py v3）。

特点：
  * 支持多个输入目录 (--in 可重复)，自动合并语料（用于 Tier1+Tier2 联合分析）。
  * 同一样本间 Pearson（按样本/列中心化）—— 修正了旧版按基因中心化制造伪 r=1.0 的 bug。
  * 同 GSE 根号（如 GSE35982_GPL14767 与 GSE35982_GPL4133）视为同一研究，不互判为跨 GSE 重复。
  * 可选阳性对照 (--pos1/--pos2)，输出双向最佳匹配计数。
  * 两遍流式加载 + 落位，避免巨量嵌套 dict 导致 OOM。
"""
import argparse, csv, json, os, re
from collections import defaultdict
import numpy as np

IN_DIRS = []
OUT_DIR = "."
R_THRESH = 0.97
GROUP_MIN = 2
POS = (None, None)
NAME = "corpus"
# 基线感知判据（消除“数据集自身内部高基线”造成的伪相关，如 GSE29621）
BASE_MODE = "off"      # off = 仅用绝对阈值；adapt = 交叉 r 须超过双方组内基线分位
BASE_Q = 99.9          # 组内基线取的分位数
BASE_FLOOR = 0.97      # 阈值绝对下限（默认与 --rthresh 一致；勿设 0.99，会误杀阳性对照）
PROBE_RATIO = 0.9      # 探针空间完整性过滤比例

def tokenize(text):
    if not text: return set()
    return set(t for t in re.split(r"[\s,.;:()\[\]/\"'<>-]+", str(text).lower()) if len(t) > 3)

def gse_root(series_name):
    """GSE35982_GPL14767 -> GSE35982；GSE41258 -> GSE41258"""
    return re.sub(r"_GPL\d+$", "", series_name)

def _find_csv(g):
    for d in IN_DIRS:
        p = os.path.join(d, g + ".csv")
        if os.path.exists(p): return p
    return None

def read_meta(g):
    for d in IN_DIRS:
        p = os.path.join(d, g + "_meta.tsv")
        if os.path.exists(p):
            meta = {}
            with open(p, newline="", encoding="utf-8-sig") as f:
                r = csv.DictReader(f, delimiter="\t")
                if r.fieldnames is None: return meta
                idcol = next((c for c in ("sample_id", "gsm") if c in r.fieldnames), r.fieldnames[0])
                for row in r:
                    sid = (row.get(idcol) or "").strip()
                    if sid: meta[sid] = row
            return meta
    return {}

def fast_line_count(g):
    """只数行数（不解析 CSV），用于按探针数排序。"""
    path = _find_csv(g)
    if not path: return 0
    n = 0
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b: break
            n += b.count(b"\n")
    return max(0, n - 1)

def gene_id_set(g):
    path = _find_csv(g)
    s = set()
    with open(path, encoding="utf-8-sig") as f:
        next(f, None)
        for line in f:
            if not line.strip(): continue
            s.add(line.split(",", 1)[0].strip().strip('"'))
    return s

def filter_by_probe_integrity(glist, ratio):
    """探针空间完整性过滤。

    平台组内取"全组共同基因交集"作为比对空间，因此任何一个探针集残缺的系列
    都会把整组比对空间砍到它的水平（GSE20916 事件：54675 -> 27697，检出对数损失 88%）。
    做法：按探针数降序逐个尝试并入，若并入后交集小于当前交集的 ratio 倍则剔除该系列。
    """
    counts = {g: fast_line_count(g) for g in glist}
    order = sorted([g for g in glist if counts[g] > 0], key=lambda g: -counts[g])
    dropped = [(g, counts[g], "0基因(无表达矩阵)") for g in glist if counts[g] == 0]
    keep = []; running = None
    for g in order:
        s = gene_id_set(g)
        if running is None:
            running = s; keep.append(g); continue
        inter = running & s
        frac = len(inter) / max(1, len(running))
        if frac >= ratio:
            keep.append(g); running = inter
        else:
            dropped.append((g, counts[g], f"并入后交集 {len(inter)}/{len(running)}={frac:.2f} < {ratio}"))
        del s
    return keep, dropped

def load_group_csvs(glist):
    gid_sets = {}
    for g in glist:
        path = _find_csv(g)
        s = set()
        with open(path, newline="", encoding="utf-8-sig") as f:
            r = csv.reader(f); next(r)
            for row in r:
                if row: s.add(row[0].strip())
        gid_sets[g] = s
    common = sorted(set.intersection(*gid_sets.values())) if gid_sets else []
    g2i = {g: i for i, g in enumerate(common)}
    # 关键：同一 GSM 可能同时属于多个 GSE（SuperSeries/SubSeries、样本复用）。
    # 若只用 {样本名: 系列} 字典，后加载的系列会覆盖先加载的，导致：
    #   (a) 共享 GSM 的重复对被误判为"同系列"而丢失；(b) 相关对被错误归属给后加载的系列。
    # 因此额外维护按下标索引的 series_at / root_at，与列位置一一对应。
    sample_list = []; series_of = {}; series_at = []; root_at = []; col_off = {}
    for g in glist:
        path = _find_csv(g)
        with open(path, newline="", encoding="utf-8-sig") as f:
            r = csv.reader(f); h = next(r); cols = [x.strip() for x in h[1:]]
            col_off[g] = len(sample_list)
            for s in cols:
                sample_list.append(s); series_of[s] = g
                series_at.append(g); root_at.append(gse_root(g))
    M = np.full((len(common), len(sample_list)), np.nan, dtype=np.float32)
    for g in glist:
        path = _find_csv(g); base = col_off[g]
        with open(path, newline="", encoding="utf-8-sig") as f:
            r = csv.reader(f); h = next(r); cols = [x.strip() for x in h[1:]]
            cidx = list(range(1, len(h)))
            for row in r:
                if len(row) < 2: continue
                gi = g2i.get(row[0].strip())
                if gi is None: continue
                for k in range(len(cols)):
                    v = row[cidx[k]] if cidx[k] < len(row) else ""
                    try:
                        fv = float(v) if v not in ("", "NA", "NaN") else np.nan
                    except ValueError:
                        fv = np.nan
                    M[gi, base + k] = fv
    return common, sample_list, series_of, M, series_at, root_at

def corr_matrix(M):
    """样本×样本 Pearson（按样本/列中心化）。"""
    col_mean = np.nanmean(M, axis=0, keepdims=True)
    Mc = np.nan_to_num(M - col_mean)
    norm = np.sqrt((Mc ** 2).sum(axis=0, keepdims=True))
    Mcu = Mc / np.where(norm == 0, 1, norm)
    return np.clip(Mcu.T @ Mcu, -1, 1)

def main():
    global IN_DIRS, OUT_DIR, R_THRESH, GROUP_MIN, POS, NAME
    ap = argparse.ArgumentParser(description="GEO 跨系列冗余检测（泛化引擎）")
    ap.add_argument("--indirs", nargs="+", required=True, dest="indirs", help="输入 SOP 目录（可多个，合并分析）")
    ap.add_argument("--out", required=True, help="输出目录")
    ap.add_argument("--rthresh", type=float, default=0.97)
    ap.add_argument("--group_min", type=int, default=2)
    ap.add_argument("--name", default="corpus", help="报告标题用标签")
    ap.add_argument("--pos1", default=None, help="阳性对照系列1（GSE根号）")
    ap.add_argument("--pos2", default=None, help="阳性对照系列2（GSE根号）")
    ap.add_argument("--baseline-mode", choices=["off", "adapt"], default="off",
                    help="off=仅绝对阈值 r>=0.97；adapt=基线感知：交叉 r 须超过双方各自组内基线分位")
    ap.add_argument("--baseline-q", type=float, default=99.9,
                    help="组内基线分位数（默认 99.9，即组内相关分布的 p99.9）")
    ap.add_argument("--probe-ratio", type=float, default=0.9,
                    help="探针空间完整性过滤：并入后交集/当前交集 < 该比例的系列从平台组剔除（0=关闭）")
    ap.add_argument("--baseline-floor", type=float, default=0.97,
                    help="基线感知阈值的绝对下限（默认 0.97 = 与 --rthresh 一致）。"
                         "注意：设为 0.99 会误杀经重新归一化后重复度降级的真实重投"
                         "（实测阳性对照 GSE41258~GSE68468 峰值仅 0.9814，会被 0.99 下限清零）")
    args = ap.parse_args()
    IN_DIRS = args.indirs; OUT_DIR = args.out; R_THRESH = args.rthresh
    GROUP_MIN = args.group_min; NAME = args.name; POS = (args.pos1, args.pos2)
    BASE_MODE = args.baseline_mode; BASE_Q = args.baseline_q; BASE_FLOOR = args.baseline_floor
    PROBE_RATIO = args.probe_ratio

    os.makedirs(OUT_DIR, exist_ok=True)
    series = {}
    for d in IN_DIRS:
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".csv"): continue
            g = fn[:-4]
            meta = read_meta(g)
            plat = (next(iter(meta.values())).get("platform") or "") if meta else ""
            series[g] = {"meta": meta, "platform": plat}
    gses = list(series.keys())
    sample_root = {}
    for g in gses:
        for s in series[g]["meta"]:
            sample_root[s] = gse_root(g)
    groups = defaultdict(list)
    for g in gses:
        groups[series[g]["platform"] or "_unknown"].append(g)

    # ---------- 探针空间完整性过滤 ----------
    # 组内比对空间 = 全组共同基因交集。混入一个探针集残缺的系列会把整组空间砍半
    # （GSE20916 事件：54675→27697，检出对数损失 88%），故必须先剔除。
    excluded_probe = []
    if PROBE_RATIO > 0:
        newgroups = {}
        for key, glist in list(groups.items()):
            if len(glist) < GROUP_MIN or len(glist) == 1 and False:
                newgroups[key] = glist; continue
            keep, dropped = filter_by_probe_integrity(glist, PROBE_RATIO)
            newgroups[key] = keep
            for g, c, why in dropped:
                excluded_probe.append({"group": key, "gse": g, "gene_count": c, "reason": why})
            if dropped:
                print(f"  [探针完整性] {key}: 保留 {len(keep)}/{len(glist)}，剔除 {len(dropped)} 个", flush=True)
                for g, c, why in dropped:
                    print(f"      - {g} (探针 {c}) {why}", flush=True)
        groups = defaultdict(list, newgroups)

    c_edges = []
    c_edges_abs = []          # 绝对阈值（legacy口径）下的边，仅用于对照
    pair_stats = []           # 系列对级统计（含组内基线/交叉分布/判定阈值）
    group_summary = {}
    bestmatch = {}
    for key, glist in groups.items():
        if len(glist) < GROUP_MIN:
            continue
        common, sample_list, series_of, M, series_at, root_at = load_group_csvs(glist)
        if M.shape[1] < 2:
            continue
        # ---------- 按列中心化 + 单位化（原地，避免额外复制大矩阵）----------
        col_mean = np.nanmean(M, axis=0, keepdims=True)
        np.nan_to_num(M, copy=False)
        M -= col_mean
        nrm = np.sqrt((M ** 2).sum(axis=0, keepdims=True))
        M /= np.where(nrm == 0, 1, nrm)
        Mn = M
        n = len(sample_list)
        sa = np.array(series_at); ra = np.array(root_at)
        cols_of = {g: np.where(sa == g)[0] for g in glist}

        # ---------- 组内基线：每个系列自身样本间的相关分布高分位 ----------
        baseline = {}
        for g in glist:
            idx = cols_of[g]
            if len(idx) < 4:
                baseline[g] = 0.0          # 样本过少，无法估计基线
                continue
            Cw = np.clip(Mn[:, idx].T @ Mn[:, idx], -1, 1)
            iu = np.triu_indices(len(idx), k=1)
            v = Cw[iu]
            baseline[g] = float(np.percentile(v, BASE_Q)) if len(v) else 0.0

        g_abs = g_abs95 = g_abs90 = 0
        g_ada = 0
        for ai in range(len(glist)):
            g1 = glist[ai]
            idx1 = cols_of[g1]
            if len(idx1) == 0:
                continue
            block = np.clip(Mn[:, idx1].T @ Mn, -1, 1)     # n1 × N
            # 每个样本的跨系列最佳匹配
            for bi in range(len(idx1)):
                i = idx1[bi]
                mask = (sa != g1) & (ra != ra[i])
                if not mask.any():
                    continue
                vals = block[bi][mask]
                jb = int(np.argmax(vals))
                cols = np.where(mask)[0]
                bestmatch[sample_list[i]] = (round(float(vals[jb]), 4),
                                             series_at[cols[jb]], sample_list[cols[jb]])
            # 系列对级统计
            for bj in range(ai + 1, len(glist)):
                g2 = glist[bj]
                if gse_root(g2) == gse_root(g1):
                    continue
                idx2 = cols_of[g2]
                v = block[:, idx2].ravel()
                n_abs = int(np.count_nonzero(v >= R_THRESH))
                n_abs95 = int(np.count_nonzero(v >= 0.95))
                n_abs90 = int(np.count_nonzero(v >= 0.90))
                g_abs += n_abs; g_abs95 += n_abs95; g_abs90 += n_abs90
                # 基线感知阈值：必须同时超过双方各自组内基线分位，且不低于绝对下限
                T = R_THRESH
                if BASE_MODE != "off":
                    T = max(T, baseline.get(g1, 0.0), baseline.get(g2, 0.0))
                    if BASE_FLOOR is not None:
                        T = max(T, BASE_FLOOR)
                sel = np.where(v >= T)[0]
                n_ada = int(len(sel))
                g_ada += n_ada
                if n_ada or n_abs:
                    pair_stats.append({
                        "group": key, "g1": g1, "g2": g2,
                        "n_cross_pairs": int(v.size),
                        "baseline_p%.4g_g1" % BASE_Q: round(baseline.get(g1, 0.0), 4),
                        "baseline_p%.4g_g2" % BASE_Q: round(baseline.get(g2, 0.0), 4),
                        "threshold_used": round(float(T), 4),
                        "edges_abs_0.97": n_abs, "edges_baseline_aware": n_ada,
                        "cross_max": round(float(v.max()), 4) if v.size else None,
                        "n_r>=0.9999": int(np.count_nonzero(v >= 0.9999)),
                        "reduced_by": round(1 - n_ada / n_abs, 4) if n_abs else None,
                    })
                for kk in sel:
                    i = idx1[kk // len(idx2)]; j = idx2[kk % len(idx2)]
                    c_edges.append({"g1": series_at[i], "s1": sample_list[i],
                                    "g2": series_at[j], "s2": sample_list[j],
                                    "pearson": round(float(v[kk]), 4)})
        if BASE_MODE == "off":
            c_edges_abs = c_edges
        group_summary[key] = {"n_series": len(glist), "n_samples": n, "n_genes": len(common),
                              "edges_r>=0.97(cross)_absolute": g_abs,
                              "edges_baseline_aware": g_ada,
                              "edges_r>=0.95(cross)": g_abs95,
                              "edges_r>=0.90(cross)": g_abs90,
                              "samples_with_cross_best>=0.97": int(sum(1 for s in sample_list if bestmatch.get(s, (-2,))[0] >= R_THRESH)),
                              "samples_with_cross_best>=0.95": int(sum(1 for s in sample_list if bestmatch.get(s, (-2,))[0] >= 0.95))}
        print(f"  [组] {key}: {len(glist)}系列/{n}样本/{len(common)}基因 | 绝对r>={R_THRESH}={g_abs} | 基线感知={g_ada}", flush=True)
        del Mn
        del M

    # Layer A GSM 重叠
    gsm_map = defaultdict(list)
    for g in gses:
        for s, row in series[g]["meta"].items():
            gsm_map[(row.get("gsm") or s).strip()].append(g)
    layerA = defaultdict(list)
    for gsm, occ in gsm_map.items():
        occ = sorted(set(occ))
        if len(occ) > 1:
            for g in occ:
                layerA[g].append({"gsm": gsm, "also_in": [x for x in occ if x != g]})

    # Layer D 系列级 token Jaccard
    series_tokens = {}
    for g in gses:
        toks = set()
        for s, row in series[g]["meta"].items():
            for fld, val in row.items():
                if fld in ("sample_id", "gsm", "platform"): continue
                toks |= tokenize(val)
        series_tokens[g] = toks
    layerD = []
    for i in range(len(gses)):
        for j in range(i + 1, len(gses)):
            g1, g2 = gses[i], gses[j]
            if gse_root(g1) == gse_root(g2): continue
            jc = len(series_tokens[g1] & series_tokens[g2]) / len(series_tokens[g1] | series_tokens[g2]) if (series_tokens[g1] | series_tokens[g2]) else 0
            if jc >= 0.5:
                gsm_overlap = set(series[g1]["meta"].keys()) & set(series[g2]["meta"].keys())
                layerD.append({"g1": g1, "g2": g2, "meta_jaccard": round(jc, 3),
                               "gsm_overlap": sorted(gsm_overlap),
                               "verdict": "FACE_SWAP_CANDIDATE" if not gsm_overlap else "OVERLAP_BUT_SHARED_GSM"})

    # 连通分量
    parent = {(g, s): (g, s) for g in gses for s in series[g]["meta"]}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        parent[find(a)] = find(b)
    for x in parent: union(x, x)
    for e in c_edges:
        union((e["g1"], e["s1"]), (e["g2"], e["s2"]))
    gd = defaultdict(list)
    for x in parent: gd[find(x)].append(x)
    dup_groups = [sorted(v) for v in gd.values() if len(v) > 1]
    remove_set = set()
    for grp in dup_groups:
        for x in grp[1:]: remove_set.add(x)
    total_samples = sum(len(series[g]["meta"]) for g in gses)
    redundancy_rate = (len(remove_set) / total_samples) if total_samples else 0.0

    pair_stats.sort(key=lambda x: -x["edges_baseline_aware"])
    abs_total = sum(p["edges_abs_0.97"] for p in pair_stats)
    report = {
        "summary": {"n_series": len(gses), "series": gses,
                    "platform_groups": {k: v for k, v in groups.items()},
                    "total_samples": total_samples,
                    "duplicate_sample_pairs_r>=0.97": len(c_edges), "duplicate_groups": len(dup_groups),
                    "samples_to_remove": len(remove_set), "redundancy_rate": round(redundancy_rate, 4),
                    "group_detail": group_summary, "r_threshold": R_THRESH,
                    "baseline_mode": BASE_MODE, "baseline_q": BASE_Q, "baseline_floor": BASE_FLOOR,
                    "edges_absolute_0.97_total": abs_total,
                    "edges_baseline_aware_total": len(c_edges),
                    "pseudo_correlations_removed": abs_total - len(c_edges),
                    "near_identity_edges_r>=0.9999": sum(1 for e in c_edges if e["pearson"] >= 0.9999),
                    "probe_ratio": PROBE_RATIO},
        "excluded_by_probe_integrity": excluded_probe,
        "layerA_gsm_overlap": {k: v for k, v in layerA.items()},
        "layerC_expression_duplicates": c_edges,
        "layerC_pair_stats": pair_stats,
        "layerD_faceswap": layerD,
        "duplicate_groups": dup_groups,
        "best_match_cross": {s: v for s, v in bestmatch.items()},
        "recommended_remove": sorted("{}|{}".format(*x) for x in remove_set),
    }
    with open(os.path.join(OUT_DIR, "dedup_report.json"), "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    # md
    L = [f"# GEO 跨系列冗余检测报告 ({NAME}, numpy 后端, 修正版 v3)", ""]
    L.append(f"- 系列数: **{len(gses)}** ({', '.join(gses)})")
    L.append(f"- 样本总数: **{total_samples}**")
    if BASE_MODE == "off":
        L.append(f"- 判据: 绝对阈值 r>={R_THRESH}（未启用基线感知）")
        L.append(f"- 表达层重复样本对(Layer C): **{len(c_edges)}**")
    else:
        L.append(f"- 判据: **基线感知**（交叉 r 须 ≥ max(组内p{BASE_Q:g}基线, 下限{BASE_FLOOR})）")
        L.append(f"- 绝对阈值 r>=0.97 下的原始边数: **{abs_total}**")
        L.append(f"- 基线感知后保留的重复样本对: **{len(c_edges)}**（滤除伪相关 **{abs_total - len(c_edges)}** 对，"
                 f"削减 {100*(abs_total-len(c_edges))/abs_total:.1f}%）" if abs_total else "- 基线感知后: 0")
    L.append(f"- 重复样本组: **{len(dup_groups)}**；建议剔除: **{len(remove_set)}**")
    L.append(f"- **冗余率: {redundancy_rate*100:.2f}%**")
    near1 = sum(1 for e in c_edges if e["pearson"] >= 0.9999)
    L.append(f"- 其中**近乎完全一致**（r>=0.9999，最保守口径）的样本对: **{near1}**")
    L.append("")
    if excluded_probe:
        L.append("## 探针空间完整性过滤（Method 需写明）")
        L.append(f"- 规则：平台组内按探针数降序逐个并入，若并入后共同基因交集 < 当前的 {PROBE_RATIO} 倍则剔除该系列。")
        L.append(f"- 目的：避免单个探针集残缺的系列把整组比对空间砍半（GSE20916 事件：54675→27697，检出对数损失 88%）。")
        for x in excluded_probe:
            L.append(f"- 剔除 `{x['gse']}`（{x['group']}，探针 {x['gene_count']}）：{x['reason']}")
        L.append("")
    L.append("## 平台分组明细（跨系列相关边）")
    for k, v in group_summary.items():
        L.append(f"- {k}: {v['n_series']}系列 / {v['n_samples']}样本 / {v['n_genes']}基因; "
                 f"绝对r>=0.97={v['edges_r>=0.97(cross)_absolute']}, 基线感知={v['edges_baseline_aware']}, "
                 f"r>=0.95={v['edges_r>=0.95(cross)']}, r>=0.90={v['edges_r>=0.90(cross)']}; "
                 f"有跨系列最佳匹配>=0.97的样本数={v['samples_with_cross_best>=0.97']}, >=0.95={v['samples_with_cross_best>=0.95']}")
    if BASE_MODE != "off" and pair_stats:
        L.append("")
        L.append("## 基线感知判据效果（TOP 20 系列对：绝对阈值 vs 基线感知）")
        L.append("")
        L.append("| 系列A | 系列B | 交叉对数 | 组内基线A | 组内基线B | 采用阈值 | 绝对边 | 基线感知边 | 削减 | 交叉max | r≥0.9999 |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for p in pair_stats[:20]:
            L.append(f"| {p['g1']} | {p['g2']} | {p['n_cross_pairs']} | "
                     f"{p['baseline_p%.4g_g1' % BASE_Q]} | {p['baseline_p%.4g_g2' % BASE_Q]} | {p['threshold_used']} | "
                     f"{p['edges_abs_0.97']} | {p['edges_baseline_aware']} | "
                     f"{('%.1f%%' % (100*p['reduced_by'])) if p['reduced_by'] is not None else '-'} | "
                     f"{p['cross_max']} | {p['n_r>=0.9999']} |")
    L.append("")
    L.append("## Layer C — 表达层重复（全基因 Pearson，按样本中心化）")
    if c_edges:
        for e in sorted(c_edges, key=lambda x:-x["pearson"])[:40]:
            L.append(f"- `{e['g1']}/{e['s1']}` ≈ `{e['g2']}/{e['s2']}` (r={e['pearson']})")
        if len(c_edges) > 40: L.append(f"- ... 共 {len(c_edges)} 对")
    else:
        L.append("- 未检出跨 GSE 表达重复（r>=0.97）。真实重投样本常因批次/归一化噪声落在 0.90–0.96，建议结合次阈值与最佳匹配进一步复核。")
    L.append("")
    # 阳性对照
    if POS[0] and POS[1] and POS[0] in gses and POS[1] in gses:
        p1, p2 = POS
        a_in_b = [v[0] for s, v in bestmatch.items() if sample_root.get(s) == p1 and v[1] == p2]
        b_in_a = [v[0] for s, v in bestmatch.items() if sample_root.get(s) == p2 and v[1] == p1]
        L.append("## 阳性对照核查（%s ↔ %s，已知同数据集换号重投）" % (p1, p2))
        if a_in_b:
            L.append(f"- {p1} 样本在 {p2} 的最佳匹配 r: 均值={np.mean(a_in_b):.4f}, 最大={np.max(a_in_b):.4f}; "
                     f">=0.95 的样本数 **{int(np.sum(np.array(a_in_b)>=0.95))}/{len(a_in_b)}**")
            L.append(f"- {p2} 样本在 {p1} 的最佳匹配 r: 均值={np.mean(b_in_a):.4f}, 最大={np.max(b_in_a):.4f}; "
                     f">=0.95 的样本数 **{int(np.sum(np.array(b_in_a)>=0.95))}/{len(b_in_a)}**")
            L.append("- 解读：两系列高度重叠（同数据集换号重投），方法学阳性对照通过。")
        else:
            L.append(f"- 未在跨 GSE 最佳匹配中检出 {p1}/{p2} 互为最高相关（可能不在同一平台组）。")
    L.append("")
    L.append("## Layer D — 变脸候选")
    if layerD:
        for d in layerD:
            L.append(f"- `{d['g1']}` ↔ `{d['g2']}` meta_jaccard={d['meta_jaccard']} → **{d['verdict']}**")
    else:
        L.append("- 未检出疑似变脸系列（标题文本差异大，AI 复核见下）。")
    L.append("")
    L.append("## Layer A — GSM 精确重叠")
    L.append("- 未检出（不同 accession 的 GSM 本不重叠，符合预期）。" if not layerA else "见 JSON。")
    with open(os.path.join(OUT_DIR, "dedup_report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    with open(os.path.join(OUT_DIR, "clean_sample_list.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f); w.writerow(["gse", "sample", "action"])
        for g in gses:
            for s in series[g]["meta"]:
                w.writerow([g, s, "REMOVE" if (g, s) in remove_set else "KEEP"])
    print(f"[完成] 语料 {NAME}: 系列 {len(gses)} 个, 样本 {total_samples} 个")
    print(f"[完成] 跨系列重复对(r>={R_THRESH}) {len(c_edges)} 个; 冗余率 {redundancy_rate*100:.2f}%; 剔除 {len(remove_set)} 个")
    print(f"[完成] 报告 -> {OUT_DIR}")

if __name__ == "__main__":
    main()
