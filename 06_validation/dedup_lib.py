#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dedup_lib.py — 与稿件检测器【完全一致】的验证内核。

重要：稿件真实检测器 (analyze_tier.py) 的方法是
  * 全基因 Pearson 相关；
  * 按列(样本)中心化 + 单位化；
  * 基线感知门：交叉 r 须 ≥ max(组内 p99.9 基线, 绝对下限 floor)。
本文件 1:1 复刻该逻辑，并增加从原始 GEO `series_matrix.txt.gz` 直读表达矩阵的
加载器，使验证脚本能直接在 D:/AIwork/GEO 的真实数据上运行，无需任何中间产物。

>>> 不要在这里塞"极端基因指纹"之类的近似 <<<
    之前版本用 top-variance 基因做指纹省内存，但那不是稿件方法，验证结论会失真。
    真实方法在常见对 (如 54675 基因 × 两三百样本) 上内存仅数十 MB，16 GB 机器足够。
"""
import os
import gzip
import re
import numpy as np
from collections import namedtuple

Series = namedtuple("Series", ["genes", "samples", "M", "path"])

# ---------------------------------------------------------------------------
# 1. 原始 GEO 表达矩阵加载（直读 series_matrix.txt.gz，镜像 geo_matrix_to_sop_stream）
# ---------------------------------------------------------------------------
def load_series_matrix(gz_path):
    """解析 GEO series_matrix.txt.gz。

    返回 Series(genes, samples, M, path)
      M: 基因(行) × 样本(列) 的 float32 矩阵（缺失值=NaN，与 analyze_tier 一致）。
    首列 ID_REF 为基因名；其余列为样本(GSM)。严格按 tab 切分，避免 expression
    文本里的逗号破坏列对齐。
    """
    genes = []
    samples = None
    rows = []
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        in_table = False
        header = None
        for line in f:
            line = line.rstrip("\n")
            if in_table:
                if line.startswith("!series_matrix_table_end"):
                    in_table = False
                    continue
                parts = line.split("\t")
                if header is None:
                    header = parts
                    samples = [p.strip() for p in parts[1:]]
                    continue
                if len(parts) > len(header):
                    parts = parts[:len(header)]
                elif len(parts) < len(header):
                    parts = parts + [""] * (len(header) - len(parts))
                genes.append(parts[0].strip().strip('"'))
                rows.append([float(x) if (x not in ("", "NA", "NaN") and x is not None)
                             else np.nan for x in parts[1:]])
                continue
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                header = None
    if not rows:
        raise ValueError("未找到表达表: %s" % gz_path)
    M = np.array(rows, dtype=np.float32)          # 基因 × 样本（注意：不要转置）
    return Series(genes, samples, M, gz_path)


def find_series_matrix(root, gse):
    """在 root 递归查找 `GSExxxx[_GPLyyy]_series_matrix.txt.gz`，返回路径或 None。"""
    gse = re.sub(r"_GPL\d+$", "", gse)             # GSE35982_GPL14767 -> GSE35982
    hit = None
    for cur, _dirs, fs in os.walk(root):
        for fn in fs:
            if fn.endswith("series_matrix.txt.gz") and gse in fn:
                p = os.path.join(cur, fn)
                # 同 GSE 多平台时优先体积更大者（主矩阵）
                if hit is None or os.path.getsize(p) > os.path.getsize(hit):
                    hit = p
    return hit


def build_corpus_index(root):
    """扫描 root 下所有 series_matrix，返回 {GSE: path}（GSE 根号去平台后缀）。"""
    idx = {}
    for cur, _dirs, fs in os.walk(root):
        for fn in fs:
            if fn.endswith("series_matrix.txt.gz"):
                m = re.search(r"(GSE\d+)", fn)
                if m:
                    idx.setdefault(m.group(1), os.path.join(cur, fn))
    return idx


# ---------------------------------------------------------------------------
# 2. 与 analyze_tier.py 完全一致的相关系数
# ---------------------------------------------------------------------------
def column_centered_pearson(M):
    """样本×样本 Pearson，按列(样本)中心化 + 单位化。

    与 analyze_tier.corr_matrix 一致：先 nan_to_num(NaN->0) 再中心化/归一化。
    M: 基因×样本。返回 样本×样本 相关矩阵（截断到 [-1,1]）。
    """
    col_mean = np.nanmean(M, axis=0, keepdims=True)
    Mc = np.nan_to_num(M - col_mean)
    norm = np.sqrt((Mc ** 2).sum(axis=0, keepdims=True))
    Mcu = Mc / np.where(norm == 0, 1, norm)
    return np.clip(Mcu.T @ Mcu, -1, 1)


def align_two(serA, serB):
    """对两系列取共同基因并拼接为 基因×(nA+nB) 矩阵。"""
    setA, setB = set(serA.genes), set(serB.genes)
    common = sorted(setA & setB)
    if len(common) < 50:
        raise ValueError("共同基因过少(%d)，无法比对" % len(common))
    ia = {g: i for i, g in enumerate(serA.genes)}
    ib = {g: i for i, g in enumerate(serB.genes)}
    idxa = [ia[g] for g in common]
    idxb = [ib[g] for g in common]
    Mc = np.hstack([serA.M[idxa, :], serB.M[idxb, :]])
    return Mc, len(serA.samples), len(serB.samples), len(common)


def intra_baseline(C_block, q=99.9):
    """组内相关分布的上分位（作为基线感知阈值）。"""
    n = C_block.shape[0]
    if n < 4:
        return 0.0
    iu = np.triu_indices(n, k=1)
    v = C_block[iu]
    return float(np.percentile(v, q)) if len(v) else 0.0


def effective_threshold(floor, bA, bB, adapt=True):
    """基线感知阈值 = max(floor, bA, bB)（adapt）；否则仅 floor。"""
    if not adapt:
        return float(floor)
    return float(max(floor, bA, bB))


# ---------------------------------------------------------------------------
# 3. 验证指标
# ---------------------------------------------------------------------------
def pair_recall(serA, serB, floor=0.97, q=99.9, adapt=True, mode="A_to_B"):
    """在两个【已知重复】系列上算召回率（回答“检测器能找回已知重复吗”）。

    mode:
      'A_to_B' : A 的每个样本最佳跨系列匹配是否落在 B 且 >= T
      'B_to_A' : 反向
      'both'   : 双向并集
    返回 dict: nA,nB,common,max_cross,mean_best,cross_max,threshold,recall,...
    """
    Mc, nA, nB, common = align_two(serA, serB)
    C = column_centered_pearson(Mc)
    bA = intra_baseline(C[:nA, :nA], q)
    bB = intra_baseline(C[nA:, nA:], q)
    T = effective_threshold(floor, bA, bB, adapt)
    cross = C[:nA, nA:]
    bestA = np.max(cross, axis=1)                 # A->B 每个样本最佳 r
    bestB = np.max(cross, axis=0)                 # B->A
    if mode == "A_to_B":
        rec = float(np.mean(bestA >= T))
    elif mode == "B_to_A":
        rec = float(np.mean(bestB >= T))
    else:
        rec = float(np.mean(np.maximum(bestA, bestB) >= T))
    return {
        "nA": nA, "nB": nB, "common_genes": common,
        "baseline_A": round(bA, 4), "baseline_B": round(bB, 4),
        "threshold": round(T, 4),
        "cross_max": round(float(cross.max()), 4),
        "mean_best_A2B": round(float(bestA.mean()), 4),
        "mean_best_B2A": round(float(bestB.mean()), 4),
        "recall": round(rec, 4),
        "n_pass_A": int(np.sum(bestA >= T)),
        "n_pass_B": int(np.sum(bestB >= T)),
    }


def twin_recall(M_base, M_twin, floor=0.97, q=99.9, adapt=True, seed=0):
    """spike-in 基准：base 样本 i 的孪生在 nA+i，检测能否被精确召回。

    M_base, M_twin: 基因×样本，列一一对应（twin_i 是 base_i 的扰动副本）。
    返回 recall + 明细（每档扰动用）。
    """
    nA = M_base.shape[1]
    nB = M_twin.shape[1]
    if nA != nB:
        raise ValueError("base 与 twin 列数必须相等")
    Mc = np.hstack([M_base, M_twin])
    C = column_centered_pearson(Mc)
    bA = intra_baseline(C[:nA, :nA], q)
    bB = intra_baseline(C[nA:, nA:], q)
    T = effective_threshold(floor, bA, bB, adapt)
    correct = 0
    rs = []
    for i in range(nA):
        # base_i 在所有 twin 列中的最佳匹配
        vals = C[i, nA:]
        jb = int(np.argmax(vals))
        r_ij = float(vals[jb])
        ok = (jb == i) and (r_ij >= T)
        if ok:
            correct += 1
        rs.append(r_ij)
    return {
        "n": nA, "threshold": round(T, 4),
        "baseline_A": round(bA, 4), "baseline_B": round(bB, 4),
        "recall": round(correct / nA, 4),
        "mean_best_r": round(float(np.mean(rs)), 4),
        "min_best_r": round(float(np.min(rs)), 4),
    }


def sweep_threshold(serA, serB, floors=(0.95, 0.96, 0.97, 0.98, 0.99),
                    qs=(99.0, 99.9, 99.99), adapt=True):
    """阈值鲁棒性扫描：返回每组合下的 recall 与 max_cross。"""
    Mc, nA, nB, common = align_two(serA, serB)
    C = column_centered_pearson(Mc)
    bA = intra_baseline(C[:nA, :nA], 99.9)
    bB = intra_baseline(C[nA:, nA:], 99.9)
    cross = C[:nA, nA:]
    bestA = np.max(cross, axis=1)
    bestB = np.max(cross, axis=0)
    stacked = np.concatenate([bestA, bestB])      # A∪B 每个样本的最佳跨系列匹配
    out = []
    for fl in floors:
        for q in qs:
            T = effective_threshold(fl, bA, bB, adapt)
            rec = float(np.mean(stacked >= T))
            out.append({"floor": fl, "q": q, "threshold": round(T, 4),
                        "recall": round(rec, 4),
                        "cross_max": round(float(cross.max()), 4)})
    return out


# ---------------------------------------------------------------------------
# 3b. 链路复现 / 门限特异性 指标（用于稿件锚点的诚实验证）
# ---------------------------------------------------------------------------
def pair_link_report(serA, serB, floor=0.97, q=99.9, adapt=True):
    """对【已知重复 / 假阳性塌缩】类稿件锚点，给出诚实可辩护的指标。

    注意：稿件所举“已知重复”对多为“峰值 r 极高、仅部分样本真重复”的关联 accession，
    并非全样本完美复制；故不用“全员 recall=1.0”这种错误指标。

    返回：
      cross_max_peak_r      : 跨系列相关矩阵的最大值（复现稿件报告的峰值 r）
      naive_pairs_above_floor : cross≥floor 的(样本×样本)对总数（朴素匹配数，量级对应稿件 2825）
      gated_pairs           : 经基线感知门限后、双向互最佳且 ≥T 的留存对（量级对应稿件 66）
      recall_at_T           : A 中最佳跨系列匹配 ≥T 的样本比例（参考量）
    """
    Mc, nA, nB, common = align_two(serA, serB)
    C = column_centered_pearson(Mc)
    bA = intra_baseline(C[:nA, :nA], q)
    bB = intra_baseline(C[nA:, nA:], q)
    T = effective_threshold(floor, bA, bB, adapt)
    cross = C[:nA, nA:]   # 修复: 非方阵时 C[:nA,nB:] 会错位成子块
    cross_max = float(cross.max())
    naive = int(np.sum(cross >= floor))
    am = np.argmax(cross, axis=1)
    bm = np.argmax(cross, axis=0)
    gated = 0
    for i in range(nA):
        j = int(am[i])
        if int(bm[j]) == i and cross[i, j] >= T:
            gated += 1
    bestA = cross.max(axis=1)
    recall_at_T = float(np.mean(bestA >= T))
    return {
        "nA": nA, "nB": nB, "common_genes": common,
        "baseline_A": round(bA, 4), "baseline_B": round(bB, 4),
        "threshold": round(T, 4), "cross_max": round(cross_max, 4),
        "naive_pairs_above_floor": naive, "gated_pairs": gated,
        "recall_at_T": round(recall_at_T, 4),
    }


def pair_fpr(serA, serB, floor=0.97, q=99.9, adapt=True):
    """阴性对照：返回与 pair_link_report 相同的结构，主用 recall_at_T 作误报率。"""
    return pair_link_report(serA, serB, floor, q, adapt)


# ---------------------------------------------------------------------------
# 4. spike-in 扰动发生器（作用于 基因×样本 矩阵）
# ---------------------------------------------------------------------------
def make_twin(M, kind, level, seed=1):
    """对表达矩阵 M(基因×样本) 制造孪生副本，分三种真实扰动：
      noise    : 加性高斯噪声，level = SNR(dB)。SNR 越低越难。
      dropout  : 随机令 level 比例的基因缺失(NaN)，模拟平台/处理差异。
      scramble : 每个样本独立打乱 level 比例基因的取值顺序，破坏基因谱模式。
    返回扰动后的 M_twin（同形状）。
    """
    rng = np.random.default_rng(seed)
    Y = M.copy()
    if kind == "noise":
        sig = np.nanstd(M, axis=0, keepdims=True) + 1e-9
        # level = SNR(dB)。60dB≈无噪；0dB=噪声能量=信号能量（重扰动）。
        # 注意：不要把 0 当作“无扰动”——0dB 是重度噪声。
        noise_std = sig / (10 ** (level / 20.0))
        Y = M + rng.normal(0, noise_std, size=M.shape)
    elif kind == "dropout":
        p = float(level)
        if p <= 0:
            return Y
        mask = rng.random(M.shape) < p
        Y = M.copy()
        Y[mask] = np.nan
    elif kind == "scramble":
        p = float(level)
        if p <= 0:
            return Y
        n_genes = M.shape[0]
        k = max(1, int(round(p * n_genes)))
        for j in range(M.shape[1]):
            idx = rng.choice(n_genes, size=k, replace=False)
            Y[idx, j] = M[rng.permutation(idx), j]
    elif kind == "renorm":
        # 良性重处理 / 跨平台重提交：真实重复常经不同归一化或不同芯片探针平台。
        #  * 每样本随机仿射(a·x+b)：列中心化+单位化会自动消除，r 不受影响 -> 演示不变性。
        #  * 基因级【乘性】批次偏置(跨样本一致)：模拟不同探针/平台的基因特异性增益差异，
        #    该效应不会被列中心化消除，会真实降低 r，形成梯度，展示检测器能容忍多大
        #    程度的“良性跨平台分歧”。
        #  -> 该模式证明检测器对“良性预处理差异”鲁棒，而非只对字节一致副本敏感。
        rng2 = np.random.default_rng(seed + 7)
        n_genes, n_samp = M.shape
        for j in range(n_samp):
            a = rng2.uniform(0.8, 1.2)
            b = rng2.normal(0, 1.0)
            Y[:, j] = a * M[:, j] + b
        if float(level) > 0:
            # 基因特异性增益（对数正态，level=增益偏离 1.0 的 log-sd）
            gene_gain = np.exp(rng2.normal(0, float(level), size=n_genes))
            Y = Y * gene_gain[:, None]
        # 极轻技术噪声（SNR 40 dB）模拟重测
        sig = np.nanstd(M, axis=0, keepdims=True) + 1e-9
        Y = Y + rng2.normal(0, sig / (10 ** (40 / 20.0)), size=M.shape)
    else:
        raise ValueError("unknown kind: %s" % kind)
    return Y
