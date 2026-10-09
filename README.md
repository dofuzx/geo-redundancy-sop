# CRC GEO 跨系列冗余论文 · 投稿资料包
# Submission Package — Cross-GSE Redundancy in CRC Transcriptomics

> 打包日期：2026-10-09　负责人：付泽娴（河北工程大学医学院）
> 目标期刊：*Briefings in Bioinformatics*（首选）
> 稿件版本：v4（领域级 319 系列 / 25,460 样本 / 冗余率 8.52%；含实验A分类器泄漏 §3.9 与实验B RNA-seq 扩展 §3.10；8 图 5 表、16 篇参考文献）

## 目录结构

```
paper/
├── 研究主线_Research_Main_Line.md   ← 论文研究主线复盘（一句话主线/故事线/证据链/与旧版对照/待办）
├── 01_manuscript/                   ← 稿件与图表
│   ├── Manuscript_BIB_CRC_GEO_redundancy.docx   （v4，已嵌入 8 图，可直接投稿编辑）
│   ├── Manuscript_BIB_CRC_GEO_redundancy.md     （源文件，改后可由 build_docx.py 重建）
│   ├── CRC_GEO_redundancy_plan.md               （研究方案 v1.0，立题与大纲出处）
│   └── manuscript_figs/Fig1–Fig8 *.png          （300 dpi 出版级图表；Fig7=泄漏AUC，Fig8=RNA-seq扩展）
├── 02_final_reports/                ← 最终检测结果（权威数据源）
│   ├── dedup_report.md                          （人类可读报告：基线感知 TOP20、阳性对照、变脸候选）
│   ├── dedup_report.json                        （pair 级全量证据 2.9MB）
│   └── clean_sample_list.csv                    （25,460 样本 gse/sample/action 去留清单）
├── 03_case_evidence/                ← 案例证据与分型
│   ├── SOP_step4_语义复核.md                    （五层判定标准 + Juntendo 群 + GSE29621 伪相关论证）
│   └── series_metadata.json                     （系列级 title/summary/submitter/PMID/relation）
├── 04_download_lists/               ← 语料下载与清单
│   ├── Tier4_下载操作指南.md
│   ├── crc_gse_tier4_full.txt（247 条）/ crc_gse_tier4_highconf.txt（109 条）
│   ├── crc_gse_tier4_detail.tsv / crc_gse_download_list_tier3_4.txt
├── 05_scripts/                      ← 全流程可复现脚本（Python，仅依赖 numpy/python-docx）
│   ├── download_tier3_4.py          （断点续传下载器，进度条+完整性校验）
│   ├── check_gz_integrity.py        （gz 完整性检查）
│   ├── geo_matrix_to_sop_stream.py  （流式转换：series_matrix.gz → csv/meta.tsv）
│   ├── extract_series_meta.py       （系列元数据抽取）
│   ├── gen_tier4_list.py            （E-utilities 领域级系列枚举）
│   ├── analyze_tier.py              （核心检测器：四层+基线感知+探针完整性过滤）
│   ├── make_figs.py                 （本文 Fig1–6 生成脚本）
│   ├── make_fig8.py                 （Fig8 RNA-seq 扩展图脚本）
│   ├── enum_rnaseq.py / download_rnaseq.py / convert_rnaseq.py / analyze_tier_rnaseq.py
│   │                                （实验B 全流程：枚举→下载→鲁棒转换→命名空间感知检测）
│   └── build_docx.py                （md→docx 稿件构建，支持插图）
├── 06_github_upload/                ← GitHub 仓库上传操作指南（三种方式 + 认证 + 后处理）
├── 07_downstream/                   ← 加分实验操作说明（A 分类器泄漏 / B RNA-seq）+ 实验A结果与脚本
└── 08_rnaseq_extension/             ← 实验B 完整产出（检测报告 md/json、去留清单、Fig8、枚举清单）
```

## 核心数字（最终口径，勿与他版混用）

| 指标 | 值 |
|---|---|
| 语料 | 319 系列 / 25,460 样本（Tier1–4 四目录合并，无重名冲突） |
| 探针完整性剔除 | 16 系列（8 个无表达矩阵 + 8 个探针比 <0.9） |
| 绝对 r≥0.97 标记对 | 34,497 → 基线感知保留 **2,956**（削减 91.4%） |
| 近完全一致对 r≥0.9999 | 1,238 |
| 重复组 / 建议剔除样本 | 1,527 / 2,169 |
| **冗余率** | **8.52%** |
| 阳性对照 GSE41258↔GSE68468 | 154/381 与 141/147 ≥0.95，峰值 r=0.981（判据全通过） |
| Layer D 变脸候选 | 52 对 |

## 实验 B · RNA-seq 扩展核心数字（§3.10）

| 指标 | 值 |
|---|---|
| 枚举 / 有矩阵 / 解析成功 | 1,719 / 545 / **464 系列（10,321 样本）** |
| 命名空间分组 | ENSG 278 / SYMBOL 144 / OTHER 42；完整性过滤后可用 363 系列 |
| 进入检测 | 363 系列 / 5,916 样本（低深度 QC 剔除 2,289 样本） |
| 绝对 ρ≥0.97 标记对 | 42 → 基线感知保留 **6**（削减 85.7%） |
| 重复组 / 建议剔除 | 1 组（GSE236688/236691/278411 器官样培养群）/ 4 样本 |
| **RNA-seq 冗余率** | **0.07%**（microarray 为 8.52%，低两个数量级） |
| Layer D 变脸候选 / Layer A 重名 | 28 对 / 109 处 |

## 投稿前待办（状态）

1. ✅ **合作者署名与 Funding 已补齐**：稿件首页作者顺序 Bai Jing（共一作，‡）、Yu Yonggang（共一作，†）、Liu Fei、Chen Shengqiang、Liu Guoqiang、Zexian Fu（通讯，*，†）；5 个机构标注；Funding 段已写入河北省医学科研 project No. 20260723 与河北省教育厅科研 project No. ZD2022039。docx 已同步重建。
2. ⏳ `geo-redundancy-sop` GitHub 仓库：本地仓库已 `init` 并提交（30 文件）；按 `06_github_upload/GitHub_上传操作指南.md` 建空仓库并 `git push` 后，把稿件两处占位 `https://github.com/[your-GitHub-username]/geo-redundancy-sop` 替换为真实地址再推一次。
3. ✅ **下游影响实验 A 已执行并写入稿件**：tumour/normal 分类器在 GPL570 子集（9,892 样本 / 179 系列）上，naive 随机 CV AUC≈0.945，而 study-isolated（按 GSE 分组）CV 跌至≈0.668–0.708，泄漏致 AUC 膨胀≈0.24；1,824 个冗余样本泄漏率 92.7% 但去除后 AUC 变动<0.01。已入稿件 §3.9、Table 4、Fig7。DE 稳定性分支未执行。
4. ✅ **实验 B（RNA-seq 扩展）已执行并写入稿件**：1,719 枚举 → 464 系列解析成功（10,321 样本）→ 命名空间感知检测（ENSG/SYMBOL/OTHER 分组 + 多数基因空间 + 覆盖率过滤 + 日志表头/注释列防御）→ 363 系列/5,916 样本，42→6 重复边，冗余率 0.07%。已入稿件 §2.9、§3.10、Table 5、Fig8；产出在 `08_rnaseq_extension/`，脚本在 `05_scripts/`。
