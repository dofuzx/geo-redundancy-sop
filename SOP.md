# geo-redundancy-sop — Standard Operating Procedure (SOP)

**Cross-GSE sample-redundancy and "face-swapping" detection in public transcriptomics repositories**

- Version: 1.0 (2026-10-24)
- Companion manuscript: *Detecting cross-GSE sample redundancy and "face-swapping" datasets in colorectal cancer transcriptomics: a reproducible SOP and field-level quantification* (BMC Bioinformatics submission)
- License: MIT (see LICENSE)
- Reference corpus results: 319 GEO series / 25,460 samples; 34,497 pairs flagged at absolute r ≥ 0.97 → 2,956 retained by the baseline-aware gate (91.4% reduction); 1,527 duplicate groups; 2,169 samples recommended for removal; corpus redundancy rate 8.52%.

---

## 0. Scope and intended use

This SOP is a mandatory **pre-processing step** for any analysis that aggregates expression
series from public repositories (GEO, ArrayExpress/BioStudies, SRA). It detects two
failure modes that identifier-based checks cannot see:

1. **Cross-series sample redundancy** — the same biological sample deposited under
   different series, usually with re-assigned GSM identifiers.
2. **"Face-swapping" (re-accessioning)** — an entire dataset re-submitted under a new
   GSE number, often with a different submitter or contact.

It operates on **processed expression matrices plus series metadata**; raw CEL/FASTQ
files are not required. Both a microarray arm and an RNA-seq arm are specified.

Run this SOP **before** any cross-series integration, batch correction, meta-analysis,
or classifier evaluation.

## 1. Definitions

| Term | Definition |
|---|---|
| Cross-series duplicate pair | Two samples in different series whose expression vectors are (near-)identical (per-sample-centred Pearson/Spearman r ≥ threshold) |
| Face-swap candidate | Series pair with high series-level metadata Jaccard (≥ 0.5) and zero GSM overlap |
| Platform-level pseudocorrelation | High cross-series r arising from same-tissue/same-platform cohort similarity, **without** shared biological material |
| Baseline-aware gate | Acceptance rule: cross-series r must exceed **both** series' own intra-series p99.9 baseline **and** an absolute floor (0.97) |
| Redundancy rate | Samples recommended for removal ÷ total corpus samples |

## 2. Prerequisites

- Python ≥ 3.9; `numpy` (microarray arm); `scipy` (RNA-seq arm, column ranking);
  `matplotlib` (figures only). See `06_validation/requirements.txt`.
- Disk: series matrices are streamed, so memory stays constant (a 285-sample ×
  54,675-probe series processes in < 1 GB RAM).
- Network: access to the NCBI GEO FTP repository and E-utilities.

## 3. Step 1 — Fetch and standardise series with integrity checks

1. **Enumerate the corpus.** Query NCBI E-utilities (`db=gds`) with disease terms
   (e.g. `colorectal cancer`, `colon adenocarcinoma`, `CRC`), organism `Homo sapiens`,
   entry type `GSE`, dataset type as required. Record the query, date and result count.
   Assemble in tiers (pilot → disease-focused → field-level) and record which tier each
   accession came from; **no accession may be enumerated twice**.
2. **Count series per platform arm.** A dual-platform accession contributes one unit
   per platform matrix; record the rule used.
3. **Apply structural exclusions before detection** and record them explicitly:
   SuperSeries whose member series are already present (e.g. GSE17538's GPL570 arm,
   whose members GSE17536/GSE17537 are in the corpus) and series with degenerate
   matrices (e.g. GSE20916, 27,697 of 54,675 probes).
4. **Download** each `<GSE>_series_matrix.txt.gz` from the GEO FTP with a resumable,
   progress-reporting downloader; reject truncated archives (gzip integrity check).
5. **Convert** each archive by streaming to `<GSE>.csv` (genes × samples) and
   `<GSE>_meta.tsv` (per-sample metadata) with constant memory.
6. **Assign platforms per sample** from `Sample_platform_id`, never from the
   series-level declaration — this correctly splits dual-platform series.
7. **RNA-seq arm:** download submitter-provided gene-level matrices from the same FTP
   supplementary directories; parse robustly (skip log/comment lines, drop annotation
   columns by name blacklist and low-cardinality test, strip Ensembl version suffixes,
   upper-case symbols); reject series with < 5,000 gene rows, < 4 sample columns, or
   non-unique identifiers.

## 4. Step 2 — Apply the feature-space completeness filter

A single incomplete series silently truncates the common comparison space for the
whole platform group and destroys sensitivity (motivating case GSE20916: 27,697 of
54,675 probes → 88% loss of detectable pairs).

**Microarray arm (probe-space filter):** within each platform group, merge series in
descending order of feature count; exclude any series whose inclusion would reduce the
common intersection below 0.9 × its current size; exclude zero-feature series
unconditionally. Record every exclusion with its probe ratio.

**RNA-seq arm (namespace-aware filter):** partition series by identifier namespace
(Ensembl / symbol / other); build a majority gene space per namespace (features present
in ≥ 50% of the group's series); exclude series covering < 60% of their group space;
mean-impute missing features; do not compare across namespaces (report as structurally
incomparable). Add count-specific QC: drop samples whose gene-detection count is below
0.9 × the group median; remove genes detected in < 5% of samples.

## 5. Step 3 — Run the four-layer detector with the baseline-aware gate

Evaluate every cross-series sample pair and every series pair through four layers:

- **Layer A — GSM exact overlap.** Compare sample accessions across series; any sample
  appearing in more than one series is reported directly. Catches only un-concealed
  duplication.
- **Layer B — metadata fuzzy similarity.** Tokenise free-text fields (title,
  characteristics, source, contact) per sample; Jaccard similarity surfaces candidate
  near-duplicates for review. Hint layer only.
- **Layer C — expression-level duplication (core).** Reduce each sample to an
  extreme-gene fingerprint (highest-/lowest-expressed genes); fingerprint-Jaccard
  pre-screens candidate pairs; confirm by per-sample-centred Pearson correlation
  (microarray) or per-sample-centred Spearman on ranks (RNA-seq) on the full feature
  vector. Per-sample centring removes the platform/batch mean.
- **Layer D — "face-swap" detection.** For every series pair, compare the union of
  metadata tokens by Jaccard; similarity ≥ 0.5 combined with absent GSM overlap flags
  a likely re-accessioned dataset.

**Baseline-aware acceptance gate (pseudocorrelation control).** Accept a cross-series
pair as a duplicate candidate only if

```
r_cross ≥ max( q_99.9(intra-series r of series A),
               q_99.9(intra-series r of series B),
               0.97 )
```

The floor (0.97) prevents degradation for homogeneous series; the percentile term
adapts to each platform's dynamic range. **Calibrate the floor against a verified
positive control**: in the reference corpus, raising the floor to 0.99 removed half of
the genuine duplicates of the verified re-accession pair, so 0.97 was retained.

**Default parameters** (all exposed for sensitivity analysis):

| Parameter | Default | Range tested |
|---|---|---|
| Fingerprint size (extreme genes, each tail) | 100 | — |
| Fingerprint Jaccard pre-screen | 0.6 | — |
| Correlation threshold τ | 0.97 | 0.95–0.99 |
| Baseline percentile | p99.9 | p99.0–p99.9 |
| Baseline floor | 0.97 | — |
| Metadata Jaccard (Layer D) | 0.5 | 0.4–0.6 |
| Probe-space merge ratio | 0.9 | — |

Report the threshold-sensitivity grid (τ × baseline percentile) of the corpus-level
redundancy rate; in the reference corpus the rate spans 8.02%–8.56% along the fixed
p99.9 line and 8.02%–12.35% across the full grid.

## 6. Step 4 — AI or expert semantic confirmation on residual candidates

Pass candidates from Layers B–D and the correlation stage to a semantic-review step
that reads both series' titles, abstracts, submitters, contacts, sample counts and
phenotypic composition, and adjudicates each borderline cluster as:

(a) the same experiment re-accessioned → **remove**;
(b) independent but similar cohorts → **keep**;
(c) legitimate biological replication → **keep**;
(d) platform-level pseudocorrelation → **keep** (gate removes the flags, not samples).

Record every adjudication (pair, decision, rationale). Semantic review finalises the
removal/retention decision and prevents over-deletion of genuinely independent samples
that moderate correlation alone would flag.

## 7. Step 5 — Analyse only the KEEP set, and report

Outputs (see `02_final_reports/` in this package):

- `dedup_report.json` — full pair-level evidence (every retained pair with r, both
  series' baselines, threshold used; Layer A/B/D evidence; duplicate groups).
- `dedup_report.md` — human-readable summary (gate effect on top series pairs,
  positive-control check, face-swap review queue).
- `clean_sample_list.csv` — per-sample KEEP/REMOVE with provenance.

Report, at minimum: enumerated vs analysed series/samples (state the denominator);
raw flags at the absolute threshold; pairs retained after the gate; near-identical
pairs (r ≥ 0.9999); duplicate groups; samples recommended for removal; corpus
redundancy rate; face-swap candidate queue size; the typology and disposition of
retained look-alikes. Use the companion `REPORTING_CHECKLIST.md` for the full list.

## 8. Validation protocol (recommended)

1. **Constructed spike-in ground truth.** Take real samples as substrate; generate
   exact twin copies under perturbation families (additive noise, renormalisation bias,
   gene dropout, gene-value scrambling) spanning the realistic re-submission envelope
   and an adversarial envelope. Expect recall = 1.000 inside the realistic envelope and
   a smooth, monotonic degradation outside it (`06_validation/spike_in_benchmark.py`).
2. **Verified re-accession positive control.** A community-documented re-accessioned
   pair must be recovered end-to-end from primary data (peak cross-series r reproduced;
   genuine duplicates retained by the gate). Use it to calibrate the floor.
3. **Specificity / false-positive rate.** Run known-independent series pairs; expect
   naive flags to collapse under the gate and recall at the operating point to be 0.
4. **Threshold sensitivity.** Reproduce the corpus-level rate across the τ × percentile
   grid (`06_validation/corpus_sweep.py`); confirm the operating point is not a knife-edge.

## 9. Package layout

```
geo-redundancy-sop/
├── README.md                        ← this package overview and core numbers
├── SOP.md                           ← this document
├── REPORTING_CHECKLIST.md           ← reporting checklist for cross-GSE meta-analyses
├── Supplementary_Table_S1.csv       ← complete 319-series corpus accession list
├── Research_Main_Line.md            ← study main line and evidence chain
├── 01_manuscript/                   ← manuscript figures
├── 02_final_reports/                ← dedup_report.json/.md, clean_sample_list.csv
├── 03_case_evidence/                ← semantic-review record, series metadata
├── 04_download_lists/               ← tiered enumeration and download lists
├── 05_scripts/                      ← acquisition, conversion, detector, RNA-seq arm
├── 06_validation/                   ← independent-validation toolkit (spike-in,
│                                      threshold robustness, corpus sweep) + results
├── 07_downstream/                   ← downstream impact experiments (classifier leakage)
└── 08_rnaseq_extension/             ← RNA-seq extension outputs
```

## 10. Citation

If you use this SOP, cite the companion manuscript and this package
(https://github.com/dofuzx/geo-redundancy-sop).
