# CRC GEO 跨系列冗余论文 · 投稿资料包
# Submission Package — Cross-GSE Redundancy in CRC Transcriptomics

> 打包日期：2026-10-09　负责人：付泽娴（河北工程大学医学院）
> 目标期刊：*Briefings in Bioinformatics*（首选）
> 稿件版本：v3（领域级 319 系列 / 25,460 样本 / 冗余率 8.52%，含 6 图 3 表、16 篇参考文献）

## 目录结构

```
paper/
├── 研究主线_Research_Main_Line.md   ← 论文研究主线复盘（一句话主线/故事线/证据链/与旧版对照/待办）
├── 01_manuscript/                   ← 稿件与图表
│   ├── Manuscript_BIB_CRC_GEO_redundancy.docx   （v3，已嵌入 6 图，可直接投稿编辑）
│   ├── Manuscript_BIB_CRC_GEO_redundancy.md     （源文件，改后可由 build_docx.py 重建）
│   ├── CRC_GEO_redundancy_plan.md               （研究方案 v1.0，立题与大纲出处）
│   └── manuscript_figs/Fig1–Fig6 *.png          （300 dpi 出版级图表）
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
└── 05_scripts/                      ← 全流程可复现脚本（Python，仅依赖 numpy/python-docx）
    ├── download_tier3_4.py          （断点续传下载器，进度条+完整性校验）
    ├── check_gz_integrity.py        （gz 完整性检查）
    ├── geo_matrix_to_sop_stream.py  （流式转换：series_matrix.gz → csv/meta.tsv）
    ├── extract_series_meta.py       （系列元数据抽取）
    ├── gen_tier4_list.py            （E-utilities 领域级系列枚举）
    ├── analyze_tier.py              （核心检测器：四层+基线感知+探针完整性过滤）
    ├── make_figs.py                 （本文 6 图生成脚本）
    └── build_docx.py                （md→docx 稿件构建，支持插图）
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

## 投稿前待办

1. 合作者/署名与 Funding 补齐（稿件中标注 [To be completed]）。
2. `geo-redundancy-sop` GitHub 仓库定稿，替换文中占位 URL。
3. （可选加分）下游影响实验：DE 稳定性 + 分类器泄漏 AUC 膨胀。
4. （可选）RNA-seq 扩展：NCBI 统一 counts 复跑同一管线。
