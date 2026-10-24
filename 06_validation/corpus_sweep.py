#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
corpus_sweep.py — 完整 319 系列语料上的检测器数值复现 + 阈值敏感性扫描（任务 B 真实数据）

严格复现 paper/05_scripts/analyze_tier.py 的数值内核：
  * 全基因 Pearson（按样本/列中心化 + 单位化，float32，NaN→0）
  * 平台组内共同基因交集；探针空间完整性过滤（ratio=0.9）
  * 组内基线分位数；基线感知门 = max(τ, bA_q, bB_q)
并在一次扫描中同时输出：
  * τ ∈ {0.95,0.96,0.97,0.98,0.99} × q ∈ {99.0,99.5,99.9} 的门限对数 / 冗余率
  * pair-level evidence CSV（r ≥ 0.95），供 sweep_sensitivity.py 独立复算
  * 与 dedup_report.json 的逐项对账（validation）

数据来源（319 系列全部本地化）：
  G:/YYG/data/GEO_Tire3_sop, G:/YYG/data/GEO_Tire4_sop,
  J:/CRC_GEO_2026T1_sop_v2,  J:/CRC_GEO_2026T2_sop_v2
"""
import csv
import json
import os
import re
import sys
import time

import numpy as np

T0 = time.time()


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


DIRS = [r"G:/YYG/data/GEO_Tire3_sop", r"G:/YYG/data/GEO_Tire4_sop",
        r"J:/CRC_GEO_2026T1_sop_v2", r"J:/CRC_GEO_2026T2_sop_v2"]
REPORT = r"G:/YYG/paper/02_final_reports/dedup_report.json"
OUTDIR = r"G:/YYG/supplement/sweep_output"

TAUS = [0.95, 0.96, 0.97, 0.98, 0.99]
QS = [99.0, 99.5, 99.9]
PROBE_RATIO = 0.9
NAIVE_TS = [0.90, 0.95, 0.97, 0.9999]
PAIR_STORE_FLOOR = 0.95      # pair-level evidence 只存 r >= 0.95（门限 >= 0.95，足够复算全部单元格）

os.makedirs(OUTDIR, exist_ok=True)


def gse_root(g):
    return re.sub(r"_GPL\d+$", "", g)


def sel_g2i(genes, g2i):
    return [g2i[x] for x in genes if x in g2i]


ONLY = set(a for a in sys.argv[1:] if not a.startswith("-"))


# ---------------------------------------------------------------------------
# 1. 语料装配
# ---------------------------------------------------------------------------
def locate(series):
    loc = {}
    for d in DIRS:
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if fn.endswith(".csv"):
                loc.setdefault(fn[:-4], os.path.join(d, fn))
    missing = [g for g in series if g not in loc]
    if missing:
        raise SystemExit("缺少矩阵: %s" % missing)
    return loc


def fast_line_count(path):
    n = 0
    with open(path, "rb") as f:
        while True:
            b = f.read(1 << 22)
            if not b:
                break
            n += b.count(b"\n")
    return max(0, n - 1)


def gene_id_set(path):
    """与 analyze_tier.gene_id_set 一致：首列 split + strip 引号。"""
    s = set()
    with open(path, encoding="utf-8-sig") as f:
        next(f, None)
        for line in f:
            if not line.strip():
                continue
            s.add(line.split(",", 1)[0].strip().strip('"'))
    return s


def probe_filter(glist, counts, gsets):
    """复刻 filter_by_probe_integrity。"""
    dropped = []
    order = sorted([g for g in glist if counts[g] > 0], key=lambda g: -counts[g])
    dropped += [(g, counts[g], "0基因(无表达矩阵)") for g in glist if counts[g] == 0]
    keep, running = [], None
    for g in order:
        s = gsets[g]
        if running is None:
            running = s
            keep.append(g)
            continue
        inter = running & s
        frac = len(inter) / max(1, len(running))
        if frac >= PROBE_RATIO:
            keep.append(g)
            running = inter
        else:
            dropped.append((g, counts[g], "并入后交集 %d/%d=%.2f < %.2f" % (len(inter), len(running), frac, PROBE_RATIO)))
    return keep, dropped


def read_matrix(path):
    """-> (genes[去引号], samples[保留文件内引号风格], M float32 基因×样本)。

    样本名保留 csv.reader 的原始输出（含引号），与 dedup_report.json 的命名一致。
    """
    import pandas as pd
    with open(path, newline="", encoding="utf-8-sig") as f:
        h = next(csv.reader(f))
    names = ["__idx__"] + h[1:]
    dt = {c: np.float32 for c in names[1:]}
    try:
        df = pd.read_csv(path, index_col=0, skiprows=1, names=names, dtype=dt,
                         na_values=["", "NA", "NaN"], keep_default_na=False)
    except Exception:
        return read_matrix_slow(path)
    df = df[~df.index.duplicated(keep="last")]
    genes = [str(x).strip().strip('"') for x in df.index]
    samples = [str(c) for c in df.columns]
    M = np.ascontiguousarray(df.to_numpy(dtype=np.float32))
    return genes, samples, M


def read_matrix_slow(path):
    """兜底：与 analyze_tier.load_group_csvs 相同的逐行解析。"""
    genes, samples, rows = [], None, []
    with open(path, newline="", encoding="utf-8-sig") as f:
        r = csv.reader(f)
        h = next(r)
        samples = [x.strip() for x in h[1:]]
        for row in r:
            if len(row) < 2:
                continue
            genes.append(row[0].strip().strip('"'))
            rows.append([float(x) if x not in ("", "NA", "NaN") else np.nan for x in row[1:]])
    M = np.array(rows, dtype=np.float32)
    if M.shape[1] < len(samples):
        M = np.hstack([M, np.full((M.shape[0], len(samples) - M.shape[1]), np.nan, dtype=np.float32)])
    return genes, samples, M


# ---------------------------------------------------------------------------
# 2. 并查集（每单元格一套）
# ---------------------------------------------------------------------------
class UF:
    __slots__ = ("p",)

    def __init__(self, n):
        self.p = list(range(n))

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


# ---------------------------------------------------------------------------
# 3. 主流程
# ---------------------------------------------------------------------------
def main():
    rep = json.load(open(REPORT, encoding="utf-8"))
    series = list(rep["summary"]["series"])            # 原始 gses 顺序
    plat_of = {}
    for plat, gs in rep["summary"]["platform_groups"].items():
        for g in gs:
            plat_of[g] = plat
    for e in rep["excluded_by_probe_integrity"]:       # 过滤前平台归属
        plat_of[e["gse"]] = e["group"]

    loc = locate(series)
    log("语料装配完成: %d 系列" % len(series))

    # 预过滤平台组
    groups = {}
    for g in series:
        groups.setdefault(plat_of.get(g) or "_unknown", []).append(g)

    # 全局样本编号
    gidx = {}
    all_samples = []
    for g in series:
        with open(loc[g], newline="", encoding="utf-8-sig") as f:
            h = next(csv.reader(f))
        for s in h[1:]:
            s = s.strip()
            gidx[(g, s)] = len(all_samples)
            all_samples.append((g, s))
    N = len(all_samples)
    log("全局样本数 %d（报告 %d）" % (N, rep["summary"]["total_samples"]))

    # 15 个单元格的并查集 + 计数
    cells = [(t, q) for t in TAUS for q in QS]
    ufs = {c: UF(N) for c in cells}
    cell_pairs = {c: 0 for c in cells}          # 门限后对数
    cell_naive = {t: 0 for t in NAIVE_TS}       # 绝对阈值对数
    pair_evidence = []                          # r >= 0.95 的对
    bestmatch = {}
    group_results = {}
    excl_all = []

    for key, glist in groups.items():
        if len(glist) < 2:
            continue
        if ONLY and key not in ONLY:
            continue
        t0g = time.time()
        counts = {g: fast_line_count(loc[g]) for g in glist}
        gsets = {g: gene_id_set(loc[g]) for g in glist}
        log("[组 %s] %d 系列；基因集扫描 %.0fs" % (key, len(glist), time.time() - t0g))
        keep, dropped = probe_filter(glist, counts, gsets)
        for g, c, why in dropped:
            excl_all.append({"group": key, "gse": g, "gene_count": c, "reason": why})
        if not keep:
            continue
        log("[组 %s] %d 系列 -> 过滤后 %d（剔除 %d）" % (key, len(glist), len(keep), len(dropped)))
        for g, c, why in dropped:
            log("    - %s (探针 %d) %s" % (g, c, why))

        common = sorted(set.intersection(*[gsets[g] for g in keep]))
        g2i = {gg: i for i, gg in enumerate(common)}
        sample_list, series_at, root_at, col_off, nsamp = [], [], [], {}, {}
        for g in keep:
            with open(loc[g], newline="", encoding="utf-8-sig") as f:
                h = next(csv.reader(f))
            nsamp[g] = len(h) - 1
            col_off[g] = len(sample_list)
            for s in h[1:]:
                s = s.strip()
                sample_list.append(s)
                series_at.append(g)
                root_at.append(gse_root(g))
        n = len(sample_list)
        log("  共同基因 %d；样本 %d；加载矩阵..." % (len(common), n))

        M = np.full((len(common), n), np.nan, dtype=np.float32)
        for g in keep:
            genes, samples, Mg = read_matrix(loc[g])
            assert len(samples) == nsamp[g], "样本数不一致: %s (%d vs %d)" % (g, len(samples), nsamp[g])
            sel = np.array([i for i, x in enumerate(genes) if x in g2i], dtype=np.int64)
            M[sel_g2i(genes, g2i), col_off[g]:col_off[g] + len(samples)] = Mg[sel, :]
            del Mg
        log("  矩阵就绪 %.2f GB；中心化+单位化..." % (M.nbytes / 1e9))

        col_mean = np.nanmean(M, axis=0, keepdims=True)
        Mc = M - col_mean
        del M
        np.nan_to_num(Mc, copy=False)
        nrm = np.sqrt((Mc ** 2).sum(axis=0, keepdims=True))
        Mn = Mc
        Mn /= np.where(nrm == 0, 1, nrm)
        del nrm

        sa = np.array(series_at)
        ra = np.array(root_at)
        cols_of = {g: np.where(sa == g)[0] for g in keep}

        # 组内基线（三个分位）
        log("  计算组内基线分位...")
        base = {q: {} for q in QS}
        for g in keep:
            idx = cols_of[g]
            if len(idx) < 4:
                for q in QS:
                    base[q][g] = 0.0
                continue
            Cw = np.clip(Mn[:, idx].T @ Mn[:, idx], -1, 1)
            iu = np.triu_indices(len(idx), k=1)
            v = Cw[iu]
            for q in QS:
                base[q][g] = float(np.percentile(v, q)) if len(v) else 0.0

        g_naive = {t: 0 for t in NAIVE_TS}
        g_cell = {c: 0 for c in cells}
        log("  跨系列扫描（%d 系列对）..." % (len(keep) * (len(keep) - 1) // 2))
        for ai in range(len(keep)):
            g1 = keep[ai]
            idx1 = cols_of[g1]
            if len(idx1) == 0:
                continue
            block = np.clip(Mn[:, idx1].T @ Mn, -1, 1)     # n1 × N
            # 每样本跨系列最佳匹配
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
            for bj in range(ai + 1, len(keep)):
                g2 = keep[bj]
                if gse_root(g2) == gse_root(g1):
                    continue
                idx2 = cols_of[g2]
                v = block[:, idx2].ravel()
                for t in NAIVE_TS:
                    c = int(np.count_nonzero(v >= t))
                    cell_naive[t] += c
                    g_naive[t] += c
                for (t, q) in cells:
                    T = max(t, base[q][g1], base[q][g2])
                    sel = np.nonzero(v >= T)[0]
                    cell_pairs[(t, q)] += len(sel)
                    g_cell[(t, q)] += len(sel)
                    if len(sel):
                        i_idx = idx1[sel // len(idx2)]
                        j_idx = idx2[sel % len(idx2)]
                        uf = ufs[(t, q)]
                        for a, b in zip(i_idx.tolist(), j_idx.tolist()):
                            uf.union(gidx[(g1, sample_list[a])], gidx[(g2, sample_list[b])])
                # pair-level evidence (r>=0.95)
                hi = np.nonzero(v >= PAIR_STORE_FLOOR)[0]
                if len(hi):
                    rv = v[hi]
                    ii = idx1[hi // len(idx2)]
                    jj = idx2[hi % len(idx2)]
                    bA = base[99.9][g1]
                    bB = base[99.9][g2]
                    for a, b, r in zip(ii.tolist(), jj.tolist(), rv.tolist()):
                        pair_evidence.append((g1, sample_list[a], g2, sample_list[b],
                                              round(r, 4), round(bA, 4), round(bB, 4)))
            del block
            if ai % 40 == 0:
                log("    ...%d/%d 系列 | 绝对r>=0.97 累计 %d | 门限(0.97,99.9) 累计 %d"
                    % (ai + 1, len(keep), cell_naive[0.97], cell_pairs[(0.97, 99.9)]))
        del Mn

        group_results[key] = {"n_series": len(keep), "n_samples": n, "n_genes": len(common),
                              "naive": {str(t): g_naive[t] for t in NAIVE_TS},
                              "cells": {"%s|%s" % c: g_cell[c] for c in cells}}
        log("  [组 %s 完成] %s" % (key, json.dumps(group_results[key]["naive"])))
        with open(os.path.join(OUTDIR, "_checkpoint.json"), "w", encoding="utf-8") as f:
            json.dump({"groups": group_results, "excluded": excl_all,
                       "cell_pairs": {"%s|%s" % c: v for c, v in cell_pairs.items()},
                       "cell_naive": {str(t): v for t, v in cell_naive.items()},
                       "n_bestmatch": len(bestmatch), "n_pair_evidence": len(pair_evidence)},
                      f, ensure_ascii=False)

    # ---------------- 汇总 ----------------
    def cell_stats(cell):
        uf = ufs[cell]
        seen = {}
        for x in range(N):
            r = uf.find(x)
            seen.setdefault(r, 0)
            seen[r] += 1
        groups_n = sum(1 for v in seen.values() if v > 1)
        removed = sum(v - 1 for v in seen.values() if v > 1)
        return groups_n, removed

    log("汇总各单元格重复组...")
    cell_final = {}
    for c in cells:
        gn, rm = cell_stats(c)
        cell_final[c] = {"pairs": cell_pairs[c], "groups": gn, "removed": rm,
                         "rate_pct": round(100.0 * rm / N, 3)}
        log("  tau=%.2f q=%.1f: pairs=%d groups=%d removed=%d rate=%.3f%%"
            % (c[0], c[1], cell_pairs[c], gn, rm, 100.0 * rm / N))

    # bestmatch 统计
    bm_97 = sum(1 for v in bestmatch.values() if v[0] >= 0.97)
    bm_95 = sum(1 for v in bestmatch.values() if v[0] >= 0.95)

    # ---------------- 对账 ----------------
    v = []
    v.append("=" * 72)
    v.append("语料阈值扫描 — 与 dedup_report.json 对账（tau=0.97, q=99.9）")
    v.append("=" * 72)
    gd = rep["summary"]["group_detail"]
    for key in ["GPL570", "GPL96", "GPL10558", "GPL4133"]:
        if key not in group_results:
            continue
        mine_n, mine_b = group_results[key]["naive"]["0.97"], group_results[key]["cells"]["0.97|99.9"]
        rep_n = gd[key]["edges_r>=0.97(cross)_absolute"]
        rep_b = gd[key]["edges_baseline_aware"]
        ok = "OK " if (mine_n == rep_n and mine_b == rep_b) else "!!!"
        v.append("%s %-9s 绝对r>=0.97: %d vs %d | 基线感知: %d vs %d" % (ok, key, mine_n, rep_n, mine_b, rep_b))
    tot_n = cell_naive[0.97]
    tot_b = cell_pairs[(0.97, 99.9)]
    v.append("%s 全体     绝对r>=0.97: %d vs %d | 基线感知: %d vs %d"
             % ("OK " if (tot_n == rep["summary"]["edges_absolute_0.97_total"]
                          and tot_b == rep["summary"]["edges_baseline_aware_total"]) else "!!!",
                tot_n, rep["summary"]["edges_absolute_0.97_total"],
                tot_b, rep["summary"]["edges_baseline_aware_total"]))
    for key in ["GPL570", "GPL96", "GPL10558", "GPL4133"]:
        if key not in group_results:
            continue
        v.append("   %-9s 绝对r>=0.95: %d vs %d | r>=0.90: %d vs %d"
                 % (key, group_results[key]["naive"]["0.95"], gd[key]["edges_r>=0.95(cross)"],
                    group_results[key]["naive"]["0.9"], gd[key]["edges_r>=0.90(cross)"]))
    cf97 = cell_final[(0.97, 99.9)]
    v.append("%s 重复组 %d vs %d | 建议剔除 %d vs %d | 冗余率 %.4f%% vs %.2f%%"
             % ("OK " if (cf97["groups"] == rep["summary"]["duplicate_groups"]
                          and cf97["removed"] == rep["summary"]["samples_to_remove"]) else "!!!",
                cf97["groups"], rep["summary"]["duplicate_groups"],
                cf97["removed"], rep["summary"]["samples_to_remove"],
                cf97["rate_pct"], 100 * rep["summary"]["redundancy_rate"]))
    v.append("%s 近乎一致对(r>=0.9999) %d vs %d"
             % ("OK " if cell_naive[0.9999] == rep["summary"]["near_identity_edges_r>=0.9999"] else "!!!",
                cell_naive[0.9999], rep["summary"]["near_identity_edges_r>=0.9999"]))
    v.append("%s bestmatch 条目 %d vs %d | 跨最佳>=0.97 %d vs %d | >=0.95 %d vs %d"
             % ("OK " if len(bestmatch) == len(rep["best_match_cross"]) else "!!!",
                len(bestmatch), len(rep["best_match_cross"]), bm_97,
                sum(1 for x in rep["best_match_cross"].values() if x[0] >= 0.97),
                bm_95, sum(1 for x in rep["best_match_cross"].values() if x[0] >= 0.95)))
    v.append("%s 探针过滤剔除 %d vs %d" % ("OK " if len(excl_all) == len(rep["excluded_by_probe_integrity"]) else "!!!",
                                        len(excl_all), len(rep["excluded_by_probe_integrity"])))
    txt = "\n".join(v)
    print(txt)
    with open(os.path.join(OUTDIR, "validation.txt"), "w", encoding="utf-8") as f:
        f.write(txt + "\n")

    # ---------------- 输出 ----------------
    with open(os.path.join(OUTDIR, "sweep_corpus_table.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["correlation_threshold", "gated_pairs", "retained_samples", "redundancy_rate_pct"])
        for t in TAUS:
            c = cell_final[(t, 99.9)]
            w.writerow([t, c["pairs"], c["removed"], "%.2f" % c["rate_pct"]])

    with open(os.path.join(OUTDIR, "sweep_matrix.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["tau"] + ["p%g" % q for q in QS])
        for t in TAUS:
            w.writerow([t] + ["%.3f" % cell_final[(t, q)]["rate_pct"] for q in QS])

    with open(os.path.join(OUTDIR, "sweep_matrix_pairs.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["tau"] + ["p%g" % q for q in QS])
        for t in TAUS:
            w.writerow([t] + [cell_final[(t, q)]["pairs"] for q in QS])

    with open(os.path.join(OUTDIR, "pair_level_evidence.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["r", "p99_9_A", "p99_9_B", "sample_a", "sample_b", "jaccard", "series_a", "series_b"])
        for row in pair_evidence:
            w.writerow([row[4], row[5], row[6], row[1], row[3], 0.0, row[0], row[2]])

    with open(os.path.join(OUTDIR, "corpus_sweep_report.json"), "w", encoding="utf-8") as f:
        json.dump({"cell_final": {"%s|%s" % c: cell_final[c] for c in cells},
                   "cell_naive": {str(t): v2 for t, v2 in cell_naive.items()},
                   "group_results": group_results,
                   "excluded": excl_all,
                   "n_bestmatch": len(bestmatch),
                   "n_pair_evidence": len(pair_evidence),
                   "n_samples": N,
                   "validation": txt}, f, ensure_ascii=False, indent=2)
    log("输出完成 -> %s" % OUTDIR)


if __name__ == "__main__":
    main()
