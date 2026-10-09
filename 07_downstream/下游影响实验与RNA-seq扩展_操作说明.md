# 加分实验操作说明：下游影响（分类器泄漏 AUC 膨胀）与 RNA-seq 扩展

> **【执行状态 · 2026-10-09】实验 A 已执行**，脚本 `downstream_leakage.py`，结果 `downstream_leakage_results.json`，图 `Fig7_leakage_auc.png`，并已在稿件 §3.9 / Table 4 / Fig7 报告。核心结果：GPL570 子集（9,892 样本 / 179 系列，2,102 正常 / 7,790 肿瘤）上，naive 随机 CV AUC≈0.945，study-isolated（按 GSE 分组）CV 跌至≈0.668–0.708，泄漏致 AUC 膨胀≈0.237；1,824 个冗余样本泄漏率 92.7%，但去除后 AUC 变动<0.01（说明膨胀主因是 study mixing 而非样本重复）。实验 B（RNA-seq 扩展）尚未执行，见下文 §B。

> 本文件是**操作说明**，供你决定是否纳入论文。两步实验均为「可选加分项」：
> - 实验 A 直接呼应稿件 §2.8 / §4.3 局限⑤，可发展为独立的「下游影响」Results 小节；
> - 实验 B 呼应稿件 §2.1 末段与 §4.3 局限④，把结论从 microarray 推广到 RNA-seq。
> 以下给出目标、数据准备、方法、伪代码与论文写法，**不自动执行**，由你定夺。

---

## 实验 A：分类器泄漏导致 AUC 膨胀（下游影响量化）

### A.1 科学问题
冗余样本若同时出现在训练集与测试集（数据泄漏，data leakage），会让机器学习分类器的性能（AUC）被**人为夸大**。本实验要定量证明：不做去冗余（KEEP 全量）得到的「高 AUC」有一部分来自泄漏，而非真实泛化能力。

### A.2 数据准备
- 输入：本包 `02_final_reports/clean_sample_list.csv`（25,460 行，含 `sample, gse, action` 三列，`action∈{KEEP, REMOVE}`）。
- 表达矩阵：由各 `GSE.csv`（`01_manuscript/manuscript_figs` 之外，原始 SOP 产物）按 sample 拼成 `features × samples` 矩阵 `X`，标签 `y` 来自每个系列的表型（肿瘤/正常，或分子亚型）。
- 构建两个语料：
  - **Corpus-RAW**：含全部 25,460 样本（不去冗余）。
  - **Corpus-CLEAN**：仅保留 `action == KEEP` 的样本（剔除 2,169 冗余）。

### A.3 方法
1. **任务**：二分类（如 CRC 肿瘤 vs 正常；或 CMS 亚型判别）。
2. **两种划分**：
   - *随机划分*（默认 `train_test_split`，`stratify=y`）：会自然把同一患者的冗余样本分到两侧 → 引入泄漏。
   - *按系列划分*（group = GSE，用 `GroupShuffleSplit`）：同一系列不跨训练/测试 → 消除系列级泄漏（但仍可能有跨系列冗余，可用 Corpus-CLEAN 进一步消除）。
3. **模型**：LogisticRegression / SVM / RandomForest（用 scikit-learn，标准化特征）。
4. **评估**：5×10 折（或 10 折重复 5 次）CV，报告 **AUC（均值 ± 95% CI）**。
5. **核心对比**（论文主表）：

   | 语料 | 划分方式 | AUC (95% CI) | 说明 |
   |---|---|---|---|
   | RAW | 随机 | AUC_raw_rand | 含全部泄漏 ↑ |
   | RAW | 按系列 | AUC_raw_group | 去掉系列级泄漏 |
   | CLEAN | 随机 | AUC_clean_rand | 去掉样本级冗余 |
   | CLEAN | 按系列 | AUC_clean_group | 最干净基线 |

6. **量化指标**：
   - **泄漏膨胀量** `ΔAUC = AUC_raw_rand − AUC_clean_group`（估计被冗余虚抬的性能）。
   - **泄漏比例** `leak_frac =`（跨 train/test 的冗余样本对数）/（总测试样本），由 `clean_sample_list.csv` 的 REMOVE 组在随机划分下实际落入两侧的比例估计。

### A.4 伪代码（可直接落为脚本 `downstream_leakage.py`）
```python
import pandas as pd, numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold, GroupShuffleSplit
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

csv = pd.read_csv("02_final_reports/clean_sample_list.csv")
keep = set(csv.loc[csv.action=="KEEP","sample"])
# X: DataFrame features×samples ; y: Series indexed by sample ; gse: Series indexed by sample
def evaluate(X, y, groups, cv):
    pipe = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
    return cross_val_score(pipe, X.T, y, groups=groups, cv=cv, scoring="roc_auc")

# RAW = all samples ; CLEAN = keep only KEEP
for name, mask in [("RAW", X.columns), ("CLEAN", [s for s in X.columns if s in keep])]:
    Xs, ys, gs = X[mask], y[mask], gse[mask]
    auc_rand  = evaluate(Xs, ys, groups=None, cv=StratifiedKFold(10, shuffle=True))
    auc_group = evaluate(Xs, ys, groups=gs,  cv=GroupShuffleSplit(n_splits=10))
    print(name, "random", auc_rand.mean(), "group", auc_group.mean())
```

### A.5 论文写法建议
- 新增 **§3.8（或 §3.9）Downstream impact: classifier leakage inflates AUC**，上表 + 一段解释。
- Discussion 可强调：「8.52% 冗余若不处理，在 CV 中通过随机划分会把一部分性能『算』进模型」，呼应 Kenn 等关于 over-fitting 的警告 [7]。
- 这是本稿独有的「影响证据」，多数同类方法论文只报检测不报后果，属明显加分。

---

## 实验 B：RNA-seq 扩展（把结论推广到计数数据）

### B.1 科学问题
本稿语料为 microarray。RNA-seq 是否也存在同样规模的跨研究冗余？NCBI 对 RNA-seq 提供**统一流程计算的 count 矩阵**（稿件引文 [2] Clough 2024），可直接跨研究比较，使本 SOP 可平行迁移。

### B.2 数据获取
- 来源：GEO 中 `GPLxxx` 对应 RNA-seq 系列，或 SRA。优先用 NCBI-computed quantification（如 GEO 的 `*.tsv.gz` count 表，或 recount3/trimmed counts）。
- 枚举：复用 `05_scripts/gen_tier4_list.py` 的思路，把平台过滤改为 RNA-seq 关键词（Homo sapiens, RNA-seq, colorectal），用 NCBI E-utilities 拉 CRC RNA-seq 系列清单。
- 表达矩阵：count 矩阵（features=genes × samples），无需 probe 完整性（改用「基因检出率完整性」，见 B.4）。

### B.3 方法适配（复用 `analyze_tier.py` 四层框架）
- **Layer A（GSM 重叠）**：完全复用。
- **Layer B（元数据 Jaccard）**：复用。
- **Layer C（表达重复）**：microarray 用 Pearson；RNA-seq count 建议改用 **Spearman 相关** 或 **对数化后的 Pearson**，并对零值/低表达做下采样或过滤；指纹法（极端基因集合）可保留。
  - 阈值可保持 r ≥ 0.97（Spearman）或做平台级校准。
- **基线感知判据（§2.4）**：完全复用——跨样本相关须超过各自系列内 p99.9 基线、下限 0.97。RNA-seq 的组内基线通常更高，正好验证判据的「自适应」价值。
- **探针完整性 → 基因完整性过滤（§2.5 类比）**：改为「每样本基因检出数 ≥ 全组中位数的 90%」，剔除极低质量/稀疏计数系列，避免比对空间塌陷。

### B.4 伪代码（落为 `analyze_tier_rnaseq.py`，在现有脚本加 `--mode rnaseq` 即可）
```python
# 载入 count 矩阵 X (genes × samples)；按样本基因检出率过滤
detected = (X > 0).sum(axis=0)
keep_samples = detected[detected >= 0.9*detected.median()].index
# Layer C 用 Spearman
from scipy.stats import spearmanr
def pair_r(a,b): return spearmanr(a,b).correlation
# 其余与 analyze_tier.py 同：指纹预筛 → Spearman 确认 → 基线感知门 → 语义复核
```

### B.5 运行与报告
- 输出同结构：`dedup_report_rnaseq.{json,md}`、`clean_sample_list_rnaseq.csv`。
- 报告格式对齐 microarray 结果：RNA-seq 系列数 / 样本数 / 冗余率 / 近一致对数 / 阳性对照（如有已知重投对）/ 基线感知削减比例。

### B.6 论文写法建议
- 作为 **§3.9 RNA-seq extension（或独立小节）**，与 microarray 结果并列，强调「方法学平台无关、结论可迁移」。
- Discussion 局限④可直接改写为「RNA-seq 扩展已执行，结论一致」，显著提升普适性说服力。
- 需要额外计算资源：RNA-seq 样本数可能更大，但基因维度（~2 万）远小于探针（5 万），内存反而更省。

---

## 决策建议（供你判断）

| 实验 | 工作量 | 风险 | 加分价值 | 建议 |
|---|---|---|---|---|
| A 分类器泄漏 | 低（半天，复用现有 CSV） | 低 | 高（独有「影响」证据） | **强烈建议加入** |
| B RNA-seq 扩展 | 中（1–2 天，需枚举+重跑） | 中（需新数据下载） | 高（普适性） | 若时间允许建议加入；可作为 Supplementary |

> 两实验都不影响已完成的检测与冗余率主结论，仅作为「影响证据」与「泛化证据」补充。
> 决定后告诉我，我可据此生成对应脚本并把结果并入稿件（`downstream_leakage.py` / `analyze_tier_rnaseq.py` 及新 Results 小节）。
