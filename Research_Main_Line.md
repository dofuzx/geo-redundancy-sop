# Research Main Line — Cross-GSE Redundancy in CRC Transcriptomics

> Version: v2.0 (review edition, 2026-10-09)　PI: Zexian Fu
> Corresponding manuscript: `Manuscript_BIB_CRC_GEO_redundancy.docx` (v3, 6 figures embedded, 16 references)
> 

---

## One-Sentence Main Line

**Public GEO CRC transcriptomic data harbours non-negligible cross-series sample redundancy (8.52%), including duplicates thoroughly hidden as "face-swaps" (re-submission under new accessions); this study provides a reproducible SOP of "four-layer detection + dual safety gates + AI semantic review", delivers the first disease-specific, domain-level quantification (319 series / 25,460 samples), and demonstrates that without a baseline-aware criterion a naive correlation threshold misclassifies 91.4% of platform-level pseudo-correlations as duplicates.**

---

## 1. Story Line (Problem → Method → Evidence → Conclusion)

### 1.1 Problem (Why)

- Self-submission to public repositories with no cross-repository deduplication means the same biological sample can appear under multiple accessions [Rosikiewicz 2013: ~14% affected; DupChecker 2014: 64–231 duplicate CELs across just three colorectal cancer series].
- The most insidious form = **face-swap (accession re-submission)**: an entire series gets a new GSE number and a new contact, and its samples receive new GSM numbers → every "compare-by-identifier" check fails completely [Kenn 2020; Gemma curation experience].
- Gap: prior work consisted of **microarray-era, cross-species general-repository** estimates; **a CRC-specific, actionable redundancy rate had never been measured**.

### 1.2 Method (How) — the four-step methodological main line

1. **Four-layer detection framework** (no dependence on raw CEL files):
   - A: exact GSM overlap (catches unhidden duplicates)
   - B: fuzzy metadata Jaccard (clue layer)
   - C: extreme gene-fingerprint prescreen + per-sample centered Pearson r (core layer, catches face-swaps)
   - D: series-level metadata Jaccard + zero GSM overlap (specifically catches face-swaps)
2. **Safety gate 1: probe-space integrity filter** — series with incomplete probe sets halve the common comparison space for the whole group (GSE20916: 54,675→27,697, 88% loss of detected pairs); rule = greedy union, drop if intersection <0.9×.
3. **Safety gate 2: baseline-aware criterion (the study's largest methodological innovation)** — cross-series r must be ≥ max(intra-series p99.9 baseline of both sides, 0.97 floor). The 0.97 floor was calibrated against the positive control (0.99 would wrongly kill true duplicates).
4. **AI/LLM semantic review (SOP step-4)** — adjudicates residual candidates as "same-experiment re-submission / independent similar cohort / legitimate duplicate / pseudo-correlation" to prevent over-deletion.

### 1.3 Evidence Chain (What we found)

| Evidence                          | Content                                                                                                                                                                                                | Figure         |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------- |
| Positive control                  | GSE41258↔GSE68468 (community-documented face-swap) recovered end-to-end: 154/381 and 141/147 samples mutually best-match ≥0.95, peak r=0.981                                                           | Fig 5          |
| Three-generation face-swap chain  | Juntendo cluster: GSE18105(2009)→GSE22598(2010)→GSE32323(2011), same 34 arrays published three times at r=1.000; the third step gave all-new GSMs → **Layer A fully blind, only Layer C can catch it** | Fig 4          |
| Exact copies                      | GSE12251↔GSE23597, GSE14580↔GSE16879 (the latter partly with identical GSM numbers) show batches at r=1.000                                                                                            | §3.6 text      |
| Pseudo-correlation counterexample | GSE29621↔GSE17536: absolute threshold flagged 2,825 pairs → baseline-aware retained 66 (97.7% reduction), r≥0.9999 pairs = 0 → KEEP                                                                    | Fig 3, Table 3 |
| Domain-level totals               | 34,497 absolute flagged pairs → 2,956 retained (91.4% reduction); final redundancy rate **8.52%** (2,169/25,460), 1,527 duplicate groups, 1,238 pairs at r≥0.9999, 52 Layer D faceswap candidates      | Fig 2, §3.7    |

### 1.4 Conclusion (So what)

- **Roughly 1 in 12 CRC transcriptomic samples is substantially duplicated with a sample from another series** — a systematic sample-size inflation and "study" covariate distortion for a typical ten-series meta-analysis.
- **A naive correlation-threshold auto-deduplication would delete far more than it should** (91.4% is pseudo-correlation) — this is an important warning to the field, and the core contribution that distinguishes this paper from "yet another deduplication script".
- Call to action: enforce a deduplication SOP before any cross-series integration (a standard step at the level of batch correction [ComBat/SVA]); repositories should link samples across accessions at submission time.

---

## 2. Evolution vs. the Earlier Manuscript (21 series / 1.97%)

| Dimension       | Old v2 manuscript                                                            | This v3 manuscript (post-review)                                                                                                                                          |
| --------------- | ---------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Corpus          | 21 series / 3,045 samples (Tier1+2 pilot)                                    | **319 series / 25,460 samples** (four-layer full corpus, domain-level)                                                                                                    |
| Redundancy rate | 1.97% (60/3,045)                                                             | **8.52%** (2,169/25,460)                                                                                                                                                  |
| Methodology     | four layers + AI semantics                                                   | four layers + AI semantics + **baseline-aware criterion** + **probe-integrity filter**                                                                                    |
| Key lesson      | pseudo-correlation not recognised (GSE29621 once leaned toward misjudgement) | pseudo-correlation systematically quantified (91.4%) and eliminated by the criterion                                                                                      |
| Case evidence   | only GSE41258/68468 pair                                                     | added **Juntendo three-generation face-swap chain**, two exact-copy chains, L5 typing table                                                                               |
| Figures/tables  | 0 figures, 2 tables                                                          | **6 figures, 3 tables** (workflow, corpus composition, criterion calibration, evidence chain, positive control, probe filter)                                             |
| References      | 10 (incl. 2 weak sources: community blog + vignette)                         | **16 all peer-reviewed** (added GEO 2023 update NAR 2024, Gemma formally published 2012, CRC meta-dataset Sci Data 2021, ComBat 2007, SVA 2007, Nat Rev Genet 2010, etc.) |

---

## 3. Five-Type Redundancy Taxonomy (logic behind paper Table 2)

1. **L1 Exact copy**: batches at r=1.000, shared/shifted identifiers → delete.
2. **L2 Identity masking (face-swap)**: same data under new GSMs → delete; **only Layer C can catch it**, the best demonstration of "Layer A+C complementarity necessity" (the Juntendo chain).
3. **L3 Same-group reuse**: cohort overlap from the same submitter/institution → merge/delete.
4. **L4 Heterogeneous similarity**: different institutions, no exact repeat, criterion independent → KEEP (AI review gates it).
5. **L5 Pseudo-correlation**: surface similarity caused by elevated platform/disease baseline → KEEP (eliminated by the baseline-aware criterion at the flag level, not by deleting samples).

---

## 4. Three "Must-Write-Into-the-Paper" Methodological Lessons Established in the Review

1. **Shared-GSM attribution overwrite bug**: early attribution used a `{sample_name: series}` dictionary, where a later-loaded series overwrote an earlier one → silently lost duplicate pairs (e.g. GSE14580↔GSE16879, 30 pairs at r=1.0 once missed). Fix = attribute by column index via `series_at/root_at`. **Lesson: never use a name as a unique key in multi-sample shared-identifier scenarios.**
2. **The baseline floor must not be chosen arbitrarily**: a 0.99 floor wrongly killed the positive control (51→0); the 0.97 floor was retained after calibration. **Lesson: every threshold needs a control-calibration report.**
3. **Probe integrity must precede correlation computation**: the GSE20916 event proved incomplete series collapse the group's "common space" and lose 88% of detections. **Lesson: run a feature-space health check before any cross-series comparison.**

---


