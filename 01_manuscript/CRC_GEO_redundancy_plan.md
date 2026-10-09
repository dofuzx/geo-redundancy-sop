# 结直肠癌转录组公共数据库冗余问题投稿论文 · 研究方案与大纲
# Cross-GSE Sample Redundancy in Colorectal Cancer Transcriptomics: A Detection & Impact Study — Research Plan & Outline

> 版本：v1.0　拟定日期：2026-10-07　负责人：付泽娴（河北工程大学医学院）
> 定位：以结直肠癌（Colorectal Cancer, CRC）为疾病语料，系统量化 GEO 中跨 GSE 样本重复，并评估其对荟萃分析与标志物发现的偏移影响。

---

## 1. 研究背景与立题依据（Background & Rationale）

### 1.1 已确证的事实锚点（带引文）
- **共性且可量化**：Bgee 团队整合 GEO + ArrayExpress 时发现约 **14% 数据含重复内容**（完整/部分重复实验、芯片复用）[1]。
- **不同 accession 含相同样本是已知机制**：DupChecker 指出"不同登录号的数据集可能含与多个样本号关联的重复样本"，未去重导致**假阳性、误导性聚类、模型过拟合**[2]。
- **换壳隐蔽性强**：Kenn 等证实不同 GSE 中"相等的 CEL 文件常被赋新 GSM 编号掩盖"，简单比较无法识别[3]。
- **数据库自身踩坑**：Gemma 策展超 10,000 研究，发现同数据二次提交产生不同 GSM，"仅靠 GSM 比对抓不到，需人工策展"[4]。
- **真实案例**：社区报道 GSE41258（2012）与 GSE68468（2015）"数据集变脸术"——编号/联系人均不同但样本构成高度一致。

### 1.2 本研究的空白（Gap）
已有工作集中于**微阵列时代、跨物种通用库**。尚缺：①**特定疾病（CRC）专项、可行动的冗余率**；②同时覆盖**微阵列 + RNA-seq**（NCBI 统一量化 counts 使 RNA-seq 跨 GSE 比对首次可行，Bgee 2013 明确当时未覆盖 RNA-seq 重复）；③**量化对下游任务的真实偏移**（而非仅"标记重复"）；④可复现的**检测流程 + 去冗余精选语料库**。

### 1.3 核心假设（Hypothesis）
在 GEO 的 CRC 转录组系列中，存在**不可忽视比例的跨 GSE 样本重复/复用**；未检测去除会**虚增效应量、污染差异基因集、使分类器表观性能过拟合**；一套多层检测流程可有效量化并缓解该问题。

---

## 2. 研究目标（Objectives）
1. 构建 CRC 转录组 GEO 语料（微阵列 + RNA-seq），量化跨 GSE 冗余率。
2. 建立可复现的"多层冗余检测流程"，并分类冗余类型。
3. 量化冗余对下游分析（差异基因、标志物分类器）的偏移幅度。
4. 发布去冗余精选 CRC 语料库 + 检测工具 + 报告清单（MEMO 式）。

---

## 3. 方法（Methods）

### 3.1 语料构建（Corpus）
- 检索式：GEO 中 `colorectal cancer` / `CRC` / `colon adenocarcinoma` + `Homo sapiens` + `expression profiling`，收集 **目标 100–300 个系列**。
- 分层：微阵列（如 GPL570/HG-U133Plus2、GPL96 等）与 RNA-seq 分别统计。
- 记录每系列的 GSE、平台、样本数、提交实验室、发表年份。

### 3.2 数据获取（Acquisition）
- 微阵列：`GEOquery`（R）或 `GEOparse`（Python）拉取表达矩阵。
- RNA-seq：优先使用 **NCBI 统一计算的量化 counts**（raw count，跨研究可比）[5]，避免提交者异质处理。

### 3.3 多层冗余检测流程（Detection Pipeline，核心创新）
- **Layer A — GSM 精确重叠**：跨系列比对 Sample  accession，直接命中。
- **Layer B — 元数据模糊匹配**：样本标题、characteristics、患者 ID、提交实验室相似度（Jaccard / 编辑距离）。
- **Layer C — 表达层重复**：跨系列样本两两相关性（Pearson r）或距离，阈值（如 r > 0.99）标记疑似重复；原始 CEL 用 DupChecker 的 MD5 指纹[2]。
- **Layer D — "变脸"检测**：两 GSE 样本元数据向量集合的 Jaccard 相似度异常高（编号/联系人不同但构成一致）→ 疑似重新上号。

### 3.4 冗余类型分类（Typology）
完整重复（full duplicate）／部分重叠（partial overlap）／对照组复用（control reuse）／SuperSeries-SubSeries 嵌套／换号重投（re-accessioning）。

### 3.5 偏移量化（Impact Quantification，关键证据）
- **差异基因稳定性**：limma（微阵列）/ DESeq2（RNA-seq）在"含重复"vs"去重复"下比较 DE 基因集重叠度（Jaccard）与效应量变化。
- **分类器过拟合**：训练 tumor vs normal 分类器，演示重复样本被随机分入训练/测试时 AUC 虚高；去重复后性能回落。
- **报告指标**：总体冗余率、各平台冗余率、最常见类型、偏移幅度（如 AUC 膨胀百分点、DE 基因集变化比例）。

### 3.6 交付物（Deliverables）
- 开源 R/Python 检测流程（GitHub）。
- 去冗余精选 CRC 转录组语料库（含来源 provenance）。
- 公共数据库元分析质控报告清单（MEMO-style checklist）。

---

## 4. 论文结构大纲（Manuscript Outline）
1. **Introduction** — 公共库挖掘价值、重复问题已知证据、CRC 专项空白、假设。
2. **Methods** — 语料、获取、多层检测流程、偏移实验设计。
3. **Results**
   - 3.1 CRC 语料概览与冗余总率。
   - 3.2 冗余类型分布（含真实"变脸"案例复现）。
   - 3.3 微阵列 vs RNA-seq 冗余差异。
   - 3.4 对差异基因集与分类器性能的偏移量化。
   - 3.5 检测流程性能（对人工策展金标准的 precision/recall）。
4. **Discussion** — 成因、对精准医学标志物的影响、局限（阈值模糊、合法生物学重复、元数据缺失）、对 GEO 提交规范的启示。
5. **Conclusion & Recommendations** — 强制去重步骤、策展资源优先、报告清单。

---

## 5. 目标期刊（Tiered）
| 梯队 | 期刊 | 适配点 |
|---|---|---|
| 首选 | *Briefings in Bioinformatics* | 方法+资源+影响，IF≈9 |
| 数据_note | *Scientific Data* / *GigaScience* | 若强调"去冗余精选语料库"为数据集 |
| 方法 | *BMC Bioinformatics* / *Database* (Oxford) | 检测流程方法学 |

---

## 6. 时间线（Proposed Timeline）
- **M1**：语料构建 + 流程开发 + 小样本试点（15–20 个 GSE）。
- **M2**：全语料扫描 + 偏移量化分析。
- **M3**：语料库策展 + 初稿撰写。
- **M4**：修改 + 投稿。

---

## 7. 风险与对策（Risks）
- **"重复"阈值模糊** → 设敏感/严格双阈值，报告区间；抽手工验证金标准。
- **合法生物学重复误判** → 结合元数据区分技术重复 vs 独立样本。
- **元数据缺失限制 Layer B/D** → 以表达层（Layer C）为主、元数据为辅。

---

## 8. 参考文献（References）
[1] ROSIKIEWICZ M, COMTE A, NIKNEJAD A, et al. Uncovering hidden duplicated content in public transcriptomics data[J]. Database (Oxford), 2013, 2013: bat010. DOI:10.1093/database/bat010.
[2] SHENG Q, SHYR Y, CHEN X. DupChecker: a bioconductor package for checking high-throughput genomic data redundancy in meta-analysis[J]. BMC Bioinformatics, 2014, 15: 323. DOI:10.1186/1471-2105-15-323.
[3] KENN M, CACSIRE CASTILLO-TONG D, SINGER C F, et al. Microarray Normalization Revisited for Reproducible Breast Cancer Biomarkers[J]. BioMed Research International, 2020, 2020: 1363827. DOI:10.1155/2020/1363827.
[4] LIM S, FREUDENBERG J, AN C, et al. Curation of over 10,000 transcriptomic studies to enable data reuse[R]. bioRxiv, 2020. DOI:10.1101/2020.07.13.201442.
[5] NCBI GEO. RNA-seq quantifications from GEO (GEOquery vignette)[EB/OL]. https://bioconductor.posit.co/packages/3.24/bioc/vignettes/GEOquery/inst/doc/rnaseq.html.
