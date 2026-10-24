# Submission Package — Cross-GSE Redundancy in CRC Transcriptomics

> Packaged: 2026-10-24 (v1.0.1 release)　PI: Zexian Fu (Hebei University of Engineering, School of Medicine)
> Manuscript version: BMC Bioinformatics submission (domain-level 319 series / 25,460 samples / 8.52% redundancy; includes the downstream classifier-leakage experiment and the RNA-seq extension; 8 figures + Figure S1, 8 tables, 23 references)
> License: MIT

This package is the released implementation accompanying the manuscript
*Detecting cross-GSE sample redundancy and "face-swapping" datasets in colorectal cancer
transcriptomics: a reproducible SOP and field-level quantification*. It contains the
detection pipeline, the standard operating procedure, the reporting checklist, the
complete corpus accession list, the final detection reports, and an independent
validation toolkit.

## Directory Structure

```
geo-redundancy-sop/
├── README.md                        ← this file
├── SOP.md                           ← standard operating procedure (5 steps, parameters, validation protocol)
├── REPORTING_CHECKLIST.md           ← MIAME-inspired reporting checklist for cross-GSE meta-analyses
├── Supplementary_Table_S1.csv       ← complete 319-series corpus accession list (platform, tier, status, samples)
├── Research_Main_Line.md            ← research main line recap (story line / evidence chain / vs. prior version)
├── 01_manuscript/                   ← manuscript figures
│   └── manuscript_figs/Fig1–Fig8 *.png          (300 dpi publication-grade figures)
├── 02_final_reports/                ← final detection results (authoritative data source)
│   ├── dedup_report.md                          (human-readable report: baseline-aware TOP20, positive control, faceswap candidates)
│   ├── dedup_report.json                        (pair-level full evidence, 2.9 MB)
│   └── clean_sample_list.csv                    (25,460 samples: gse/sample/action keep-drop list)
├── 03_case_evidence/                ← case evidence and typing
│   ├── SOP_step4_语义复核.md                    (semantic-review criteria + Juntendo cluster + GSE29621 pseudocorrelation argument)
│   └── series_metadata.json                     (series-level title/summary/submitter/PMID/relation)
├── 04_download_lists/               ← corpus enumeration and download lists
│   ├── Tier4_下载操作指南.md
│   ├── crc_gse_tier4_full.txt (247 entries) / crc_gse_tier4_highconf.txt (109 entries)
│   ├── crc_gse_tier4_detail.tsv / crc_gse_download_list_tier3_4.txt
├── 05_scripts/                      ← fully reproducible pipeline (Python; numpy only for the microarray arm)
│   ├── download_tier3_4.py          (resumable downloader with progress bar + integrity check)
│   ├── check_gz_integrity.py        (gz integrity check)
│   ├── geo_matrix_to_sop_stream.py  (streaming conversion: series_matrix.gz → csv/meta.tsv)
│   ├── extract_series_meta.py       (series metadata extraction)
│   ├── gen_tier4_list.py            (E-utilities field-level series enumeration)
│   ├── analyze_tier.py              (core detector: four layers + baseline-aware gate + probe-integrity filter)
│   ├── make_figs.py / make_fig8.py  (figure generation)
│   ├── enum_rnaseq.py / download_rnaseq.py / convert_rnaseq.py / analyze_tier_rnaseq.py
│   │                                (RNA-seq arm: enumerate → download → robust convert → namespace-aware detect)
│   └── build_docx.py                (md→docx manuscript builder, supports figures)
├── 06_validation/                   ← independent-validation toolkit (see 06_validation/README.md)
│   ├── dedup_lib.py                 (verification kernel identical to the manuscript detector)
│   ├── spike_in_benchmark.py        (constructed ground truth, perturbation families)
│   ├── threshold_robustness.py      (link recovery, gate specificity, negative controls)
│   ├── corpus_sweep.py              (319-series numerical reproduction + threshold-sensitivity scan)
│   ├── rerun_verified.py            (authoritative re-verification run)
│   ├── external_known_dups.py / gemma_positive_control.py / verify_positive_control.py
│   ├── make_figure_S1.py / make_rebuttal_figures.py
│   ├── *_verified.csv / sweep_*.csv / corpus_sweep_report.json   (results)
│   └── figures/                     (F01–F14 validation figures, T01–T06 tables)
├── 07_downstream/                   ← downstream impact experiment (classifier leakage) + results and scripts
└── 08_rnaseq_extension/             ← RNA-seq extension outputs (detection report md/json, keep-drop list, Fig8, lists)
```

## Core Numbers (final scope — do not mix with other versions)

| Metric                                             | Value                                                                          |
| -------------------------------------------------- | ------------------------------------------------------------------------------ |
| Corpus                                             | 319 series / 25,460 samples (Tier 1–4 enumeration: 13 + 9 + 52 + 247 = 321 units, minus 2 structural exclusions → 319) |
| Probe-integrity removals                           | 16 series (8 with no expression matrix + 8 with probe ratio < 0.9)              |
| Absolute r ≥ 0.97 flagged pairs                    | 34,497 → baseline-aware retained **2,956** (91.4% reduction)                   |
| Near-identical pairs r ≥ 0.9999                    | 1,238                                                                          |
| Duplicate groups / samples recommended for removal | 1,527 / 2,169                                                                  |
| **Redundancy rate**                                | **8.52%** (2,169 / 25,460)                                                     |
| Positive control GSE41258 ↔ GSE68468               | 137/390 and 177/390 samples with reciprocal best match ≥ 0.95, peak r = 0.982 (verified re-run) |
| Layer D face-swap candidates                       | 52 series-pairs (review queue, not auto-deleted)                               |
| RNA-seq extension                                  | 363 series / 5,916 samples analysed; 42 → 6 pairs retained; redundancy rate 0.07% |

## Corpus accounting (denominators)

- Enumerated: 321 series units across four tiers (pilot 13 + 9, disease-focused 52,
  field-level 247); no accession appeared in more than one tier.
- Structural exclusions before detection: GSE17538's GPL570 arm (SuperSeries of
  GSE17536/GSE17537, already present) and GSE20916 (27,697/54,675 probes) → **319**.
- Detection scope: cross-series correlation ran within the four cross-comparable
  platform groups (286 series / 22,548 samples); 17 single-platform series without a
  comparable partner were inventoried but not compared; 16 series were excluded by the
  probe-space integrity filter. The redundancy rate uses the full-corpus denominator
  (25,460 samples).
- The threshold sweep (`06_validation/corpus_sweep.py`) reproduces the operating point
  to within 0.1% (2,954 pairs / 2,168 samples at τ = 0.97, p99.9); the final corpus run
  in `02_final_reports/` is authoritative.

## Experiment B · RNA-seq Extension Core Numbers

| Metric                                                | Value                                                                         |
| ----------------------------------------------------- | ----------------------------------------------------------------------------- |
| Enumerated / with matrix / parsed successfully        | 1,719 / 545 / **464 series (10,321 samples)**                                 |
| Namespace grouping                                    | ENSG 278 / SYMBOL 144 / OTHER 42; 363 series usable after integrity filtering |
| Entered detection                                     | 363 series / 5,916 samples (2,289 samples dropped by low-depth QC)            |
| Absolute ρ ≥ 0.97 flagged pairs                       | 42 → baseline-aware retained **6** (85.7% reduction)                          |
| Duplicate groups / samples recommended for removal    | 1 group (GSE236688/236691/278411 organoid-culture cohort) / 4 samples         |
| **RNA-seq redundancy rate**                           | **0.07%** (microarray is 8.52%, two orders of magnitude lower)                |
| Layer D face-swap candidates / Layer A name collisions | 28 pairs / 109 occurrences                                                    |

## Citation

Bai J, Yu Y, Liu F, Chen S, Liu G, Fu Z. Detecting cross-GSE sample redundancy and
"face-swapping" datasets in colorectal cancer transcriptomics: a reproducible SOP and
field-level quantification. (manuscript in review) — and this package:
https://github.com/dofuzx/geo-redundancy-sop
