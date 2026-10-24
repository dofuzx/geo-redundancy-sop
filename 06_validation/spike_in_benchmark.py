#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
工作流 2 — 可控 spike-in 扰动基准（最硬证据：纯本地、零外部依赖、非过拟合）

思路：从 D:/AIwork/GEO 取一个真实系列（默认自动挑样本数最多的 GPL570 系列），
以其中 N 个真实样本作底本，对其人为制造“孪生”副本并施加三种真实扰动，
检测“检测器（真实方法）能否把已知孪生对精确召回”。

扰动梯度：
  noise      : 加性高斯噪声，level = SNR(dB)（60→0，越小越难）——技术重测噪声（真实重复区间）
  renorm     : 每样本随机仿射 + 基因级批次偏置(level=偏差幅度)——不同归一化/跨平台重提交（良性分歧）
  dropout    : 随机令 level 比例基因缺失（0→0.5）——对抗性破坏（真实重复不会出现）
  scramble   : 每样本打乱 level 比例基因的取值（0→0.8）——对抗性破坏（真实重复不会出现）

输出：spike_in_results.csv（每档 recall/阈值/max_r）+ spike_in_recall.png
      （召回率 vs 扰动强度曲线，三色对应三种扰动）。

>>> 为什么这能回应审稿人：孪生标签是人为构造的“完美真值”，且底本是真实表达值；
    检测器在零扰动时应 recall≈1，随扰动增强平滑下降——证明它确实在“找重复”，
    而不是只在自己语料上调参。
"""
import argparse
import csv
import gzip
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dedup_lib import (load_series_matrix, find_series_matrix, build_corpus_index,
                       make_twin, twin_recall)

ROOT = "D:/AIwork/GEO"
FLOOR = 0.97
Q = 99.9
GRID = {
    "noise_SNRdB": [60, 40, 30, 20, 10, 0],
    "renorm_bias": [0.0, 0.1, 0.3, 0.5, 0.8],
    "dropout_frac": [0.0, 0.05, 0.15, 0.30, 0.50],
    "scramble_frac": [0.0, 0.1, 0.3, 0.5, 0.8],
}


def _platform_of(gz_path):
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!Series_platform_id"):
                v = line.split("\t", 1)[1].strip() if "\t" in line else ""
                return v.strip('"')                       # 去掉 GEO 引号 "GPL570"
            if line.startswith("!series_matrix_table_begin"):
                break
    return ""


def _n_samples_of(gz_path):
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                header = f.readline().rstrip("\n").split("\t")
                return max(0, len(header) - 1)
    return 0


def pick_large_series(root, min_samples=80, platform="GPL570"):
    """自动挑样本数最多且平台匹配的本地系列作为底本。"""
    idx = build_corpus_index(root)
    best = None
    for gse, p in idx.items():
        try:
            plat = _platform_of(p).strip('"')
        except Exception:
            plat = ""
        if platform and platform not in plat:
            continue
        try:
            n = _n_samples_of(p)
        except Exception:
            n = 0
        if n >= min_samples and (best is None or n > best[1]):
            best = (gse, n, p)
    return best


def run_real(root, gse, n_samples, seed):
    p = find_series_matrix(root, gse) if gse else None
    if not p:
        pick = pick_large_series(root)
        if not pick:
            raise SystemExit("[错误] 未在 %s 找到合适的 GPL570 系列" % root)
        gse, n_all, p = pick
        print("[自动选底本] %s (%d 样本)" % (gse, n_all))
    ser = load_series_matrix(p)
    n_all = ser.M.shape[1]
    if n_samples > n_all:
        n_samples = n_all
    rng = np.random.default_rng(seed)
    cols = rng.choice(n_all, size=n_samples, replace=False)
    M_base = ser.M[:, cols]                     # 真实表达值底本（基因×n_samples）
    print("[底本] %s：取 %d 个真实样本，基因数=%d" % (gse, n_samples, ser.M.shape[0]))

    rows = []
    for kind, levels in GRID.items():
        for lv in levels:
            M_twin = make_twin(M_base, kind.split("_")[0], lv, seed=seed + 1)
            o = twin_recall(M_base, M_twin, floor=FLOOR, q=Q)
            rows.append({"perturbation": kind, "level": lv,
                         "recall": o["recall"], "threshold": o["threshold"],
                         "mean_best_r": o["mean_best_r"], "min_best_r": o["min_best_r"]})
            print("  %-16s level=%-5s recall=%.3f mean_r=%.3f thr=%.4f"
                  % (kind, lv, o["recall"], o["mean_best_r"], o["threshold"]))
    return gse, rows


def run_demo():
    """自检：合成数据验证脚本逻辑（不依赖真实数据）。"""
    rng = np.random.default_rng(0)
    n, g = 120, 8000
    base = rng.normal(0, 1, size=(n, g))
    coeff = rng.normal(0, 1, size=(n, 5))
    module = rng.normal(0, 1, size=(5, g))
    M_base = (base + coeff @ module).T            # 基因×样本
    rows = []
    for kind, levels in GRID.items():
        for lv in levels:
            M_twin = make_twin(M_base, kind.split("_")[0], lv, seed=1)
            o = twin_recall(M_base, M_twin, floor=FLOOR, q=Q)
            rows.append({"perturbation": kind, "level": lv,
                         "recall": o["recall"], "mean_best_r": o["mean_best_r"]})
            print("  %-16s level=%-5s recall=%.3f" % (kind, lv, o["recall"]))
    return "DEMO", rows


def plot(gse, rows, out_png):
    order = ["noise_SNRdB", "renorm_bias", "dropout_frac", "scramble_frac"]
    colors = {"noise_SNRdB": "#c0392b", "renorm_bias": "#8e44ad",
              "dropout_frac": "#2980b9", "scramble_frac": "#27ae60"}
    # 英文标签，避免 matplotlib 缺中文字形
    xmap = {"noise_SNRdB": "SNR (dB)  ← harder", "renorm_bias": "batch bias (×std)",
            "dropout_frac": "dropout fraction", "scramble_frac": "scramble fraction"}
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.4))
    for ax, kind in zip(axes.ravel(), order):
        sub = [r for r in rows if r["perturbation"] == kind]
        xs = [r["level"] for r in sub]
        ys_r = [r["recall"] for r in sub]
        ys_m = [r.get("mean_best_r", float("nan")) for r in sub]
        ax.plot(xs, ys_m, "s--", color=colors[kind], lw=1.5, ms=5,
                label="mean best r")
        ax.plot(xs, ys_r, "o-", color=colors[kind], lw=2, ms=7, label="recall@thr")
        ax.axhline(1.0, ls=":", color="gray", lw=1)
        ax.set_title(kind, fontsize=11)
        ax.set_xlabel(xmap[kind], fontsize=9)
        ax.set_ylabel("recall / mean best r")
        ax.set_ylim(-0.05, 1.08)
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(alpha=0.3)
    fig.suptitle("Spike-in perturbation benchmark — real detector sensitivity\n"
                 "base: %s (real GEO expression matrix)" % gse, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=300)
    print("[图] 已保存 %s" % out_png)


def main():
    ap = argparse.ArgumentParser(description="spike-in 扰动基准（真实检测器）")
    ap.add_argument("--mode", choices=["real", "demo"], default="real")
    ap.add_argument("--root", default=ROOT, help="GEO 数据根目录（含 Tire3/Tire4/...）")
    ap.add_argument("--gse", default=None, help="指定底本系列（默认自动挑最大 GPL570）")
    ap.add_argument("--nsamples", type=int, default=100, help="取多少真实样本作底本")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="spike_in_results.csv")
    ap.add_argument("--png", default="spike_in_recall.png")
    args = ap.parse_args()

    if args.mode == "demo":
        gse, rows = run_demo()
    else:
        gse, rows = run_real(args.root, args.gse, args.nsamples, args.seed)

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["perturbation", "level", "recall",
                                          "threshold", "mean_best_r", "min_best_r"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("[CSV] 已保存 %s" % args.out)
    plot(gse, rows, args.png)


if __name__ == "__main__":
    main()
