# Submission Package — Cross-GSE Redundancy in CRC Transcriptomics

> Packaged: 2026-10-09　PI: Zexian Fu (Hebei University of Engineering, School of Medicine)
> Manuscript version: v4 (domain-level 319 series / 25,460 samples / 8.52% redundancy; includes Experiment A classifier leakage §3.9 and Experiment B RNA-seq extension §3.10; 8 figures, 5 tables, 16 references)

## Directory Structure

```
paper/
├── Research_Main_Line.md   ← Research main line recap (one-sentence thread / story line / evidence chain / vs. prior version / to-do)
├── 01_manuscript/                   ← Manuscript and figures
│   ├── Manuscript_BIB_CRC_GEO_redundancy.docx   (v4, 8 figures embedded, ready for editorial submission)
│   ├── Manuscript_BIB_CRC_GEO_redundancy.md     (source; rebuildable via build_docx.py)
│   ├── CRC_GEO_redundancy_plan.md               (research plan v1.0, origin of hypothesis and outline)
│   └── manuscript_figs/Fig1–Fig8 *.png          (300 dpi publication-grade figures; Fig7 = leakage AUC, Fig8 = RNA-seq extension)
├── 02_final_reports/                ← Final detection results (authoritative data source)
│   ├── dedup_report.md                          (human-readable report: baseline-aware TOP20, positive control, faceswap candidates)
│   ├── dedup_report.json                        (pair-level full evidence, 2.9MB)
│   └── clean_sample_list.csv                    (25,460 samples gse/sample/action keep-drop list)
├── 03_case_evidence/                ← Case evidence and typing
│   ├── SOP_step4_语义复核.md                    (five-layer criteria + Juntendo cluster + GSE29621 pseudo-correlation argument)
│   └── series_metadata.json                     (series-level title/summary/submitter/PMID/relation)
├── 04_download_lists/               ← Corpus download and lists
│   ├── Tier4_下载操作指南.md
│   ├── crc_gse_tier4_full.txt (247 entries) / crc_gse_tier4_highconf.txt (109 entries)
│   ├── crc_gse_tier4_detail.tsv / crc_gse_download_list_tier3_4.txt
├── 05_scripts/                      ← Fully reproducible pipeline (Python, depends only on numpy/python-docx)
│   ├── download_tier3_4.py          (resumable downloader with progress bar + integrity check)
│   ├── check_gz_integrity.py        (gz integrity check)
│   ├── geo_matrix_to_sop_stream.py  (streaming conversion: series_matrix.gz → csv/meta.tsv)
│   ├── extract_series_meta.py       (series metadata extraction)
│   ├── gen_tier4_list.py            (E-utilities domain-level series enumeration)
│   ├── analyze_tier.py              (core detector: four layers + baseline-aware + probe-integrity filter)
│   ├── make_figs.py                 (Fig1–6 generation scripts)
│   ├── make_fig8.py                 (Fig8 RNA-seq extension figure script)
│   ├── enum_rnaseq.py / download_rnaseq.py / convert_rnaseq.py / analyze_tier_rnaseq.py
│   │                                (Experiment B full pipeline: enumerate → download → robust convert → namespace-aware detect)
│   └── build_docx.py                (md→docx manuscript builder, supports figures)
├── 06_github_upload/                ← GitHub repository upload guide (three methods + auth + post-processing)
├── 07_downstream/                   ← Bonus-experiment spec (A classifier leakage / B RNA-seq) + Experiment A results and scripts
└── 08_rnaseq_extension/             ← Experiment B full outputs (detection report md/json, keep-drop list, Fig8, enumeration lists)
```

## Core Numbers (final scope — do not mix with other versions)

| Metric                                             | Value                                                                          |
| -------------------------------------------------- | ------------------------------------------------------------------------------ |
| Corpus                                             | 319 series / 25,460 samples (Tier1–4 four-directory merge, no name collisions) |
| Probe-integrity removals                           | 16 series (8 with no expression matrix + 8 with probe ratio <0.9)              |
| Absolute r≥0.97 flagged pairs                      | 34,497 → baseline-aware retained **2,956** (91.4% reduction)                   |
| Near-identical pairs r≥0.9999                      | 1,238                                                                          |
| Duplicate groups / samples recommended for removal | 1,527 / 2,169                                                                  |
| **Redundancy rate**                                | **8.52%**                                                                      |
| Positive control GSE41258↔GSE68468                 | 154/381 and 141/147 ≥0.95, peak r=0.981 (all criteria passed)                  |
| Layer D faceswap candidates                        | 52 pairs                                                                       |

## Experiment B · RNA-seq Extension Core Numbers (§3.10)

| Metric                                                | Value                                                                         |
| ----------------------------------------------------- | ----------------------------------------------------------------------------- |
| Enumerated / with matrix / parsed successfully        | 1,719 / 545 / **464 series (10,321 samples)**                                 |
| Namespace grouping                                    | ENSG 278 / SYMBOL 144 / OTHER 42; 363 series usable after integrity filtering |
| Entered detection                                     | 363 series / 5,916 samples (2,289 samples dropped by low-depth QC)            |
| Absolute ρ≥0.97 flagged pairs                         | 42 → baseline-aware retained **6** (85.7% reduction)                          |
| Duplicate groups / samples recommended for removal    | 1 group (GSE236688/236691/278411 organoid-culture cohort) / 4 samples         |
| **RNA-seq redundancy rate**                           | **0.07%** (microarray is 8.52%, two orders of magnitude lower)                |
| Layer D faceswap candidates / Layer A name collisions | 28 pairs / 109 occurrences                                                    |

##


