# Reporting Checklist for Cross-GSE Meta-Analyses

**geo-redundancy-sop companion checklist** — a MIAME-inspired minimum-information list
for any study that aggregates expression series from public repositories.
Version 1.0 (2026-10-24). License: MIT.

Rationale: public archives are populated by self-submission with no cross-study
sample-level de-duplication, so the independence of aggregated samples cannot be
assumed. Roughly one in twelve colorectal-cancer transcriptomic samples in the
reference corpus of the companion manuscript was materially duplicated across series,
and 91.4% of naive correlation flags were platform-level pseudocorrelation rather than
duplication. Both facts are invisible to identifier-level checks and to
correlation-only workflows.

Report each item below (mark "not applicable" with a reason where appropriate).

## A. Corpus definition

- [ ] **A1.** Database(s) queried and exact search terms, filters (organism, disease,
      entry/dataset type), and query date.
- [ ] **A2.** Inclusion and exclusion criteria, including minimum sample and feature
      counts and handling of SuperSeries and dual-platform accessions.
- [ ] **A3.** Series-counting rule (per accession vs per platform arm).
- [ ] **A4.** Complete accession list of the final corpus supplied as a supplementary
      table (accession, platform, tier, per-series sample count).
- [ ] **A5.** Enumerated → final corpus accounting (structural exclusions listed with
      reasons; denominator stated explicitly).

## B. Data acquisition and standardisation

- [ ] **B1.** Data source and file version (e.g. GEO FTP `_series_matrix.txt.gz`),
      access date(s), and integrity checks applied to downloads.
- [ ] **B2.** Matrix condition used (submitter-processed gene-level matrices vs
      uniformly computed counts vs raw files) and why.
- [ ] **B3.** Normalisation state of the matrices (raw counts / TPM / FPKM / log) and
      identifier namespace(s) (Ensembl / HGNC symbol / probe ID), including version
      handling.
- [ ] **B4.** Platform assignment method (per-sample annotation preferred over
      series-level declaration) and treatment of dual-platform series.

## C. Feature-space integrity

- [ ] **C1.** Common comparison space defined per platform group (or per identifier
      namespace for RNA-seq).
- [ ] **C2.** Completeness criterion applied before any cross-series correlation
      (e.g. exclude series whose inclusion would shrink the common space below 0.9×).
- [ ] **C3.** Every series excluded by the completeness filter listed with its
      coverage ratio and reason.

## D. Redundancy detection

- [ ] **D1.** Detection method described at the level of layers/mechanisms
      (identifier overlap, metadata similarity, expression-level comparison,
      series-level face-swap detection).
- [ ] **D2.** All thresholds reported: correlation threshold, fingerprint pre-screen,
      baseline percentile and floor, metadata Jaccard.
- [ ] **D3.** Whether a baseline-aware acceptance rule was used to separate true
      duplication from platform-level pseudocorrelation; if not, why not.
- [ ] **D4.** Floor/threshold calibration against a verified positive control
      (state what tightening or relaxing the threshold does to known duplicates).
- [ ] **D5.** Semantic or expert review record for borderline candidates
      (pair, decision, rationale), or a statement that none remained.

## E. Redundancy results

- [ ] **E1.** Raw flag count at the absolute threshold and the retained count after
      the baseline-aware gate (with the reduction percentage).
- [ ] **E2.** Near-identical pair count at a conservative threshold (e.g. r ≥ 0.9999).
- [ ] **E3.** Duplicate groups; samples recommended for removal; corpus redundancy
      rate with its denominator.
- [ ] **E4.** Face-swap candidate series-pairs and their disposition (reported, not
      auto-deleted).
- [ ] **E5.** Typology of detected redundancy and the disposition of retained
      look-alikes (independent same-disease cohorts kept, pseudocorrelation kept).

## F. Downstream analysis integrity

- [ ] **F1.** Cross-validation scheme stated; for multi-study data, series-grouped
      (study-isolated) splits used or the leakage risk explicitly acknowledged.
- [ ] **F2.** Effect of de-duplication on the reported endpoint quantified
      (e.g. ΔAUC under each split scheme), not merely asserted.
- [ ] **F3.** Analysis performed on the KEEP set only; removed samples and the
      deduplicated lists made available.

## G. Availability

- [ ] **G1.** Code repository (URL, version/commit) with the full pipeline, including
      the acquisition, conversion and detection tools.
- [ ] **G2.** Machine-readable pair-level evidence and per-sample KEEP/REMOVE lists.
- [ ] **G3.** SOP document and parameter defaults, sufficient for independent re-running.
- [ ] **G4.** Independent-validation evidence (constructed ground truth, verified
      positive control, specificity/false-positive rate, threshold sensitivity).

## H. Limitations

- [ ] **H1.** Scope limits stated (e.g. cross-platform and cross-namespace pairs out
      of scope; heuristic thresholds; biological replication vs redundancy).
- [ ] **H2.** Technology transfer limits stated where an extension arm exists
      (e.g. RNA-seq estimate based on submitter-processed matrices).
