# 06_validation — Independent-Validation Toolkit

Independent verification suite for the cross-GSE redundancy detector released with the
companion manuscript. Nothing in this folder is used to tune the detector; it exists to
**falsify** three failure modes: (i) recall only on self-tuned data, (ii) thresholds
fitted to the corpus, (iii) unexamined false positives on genuinely independent series.

## Scripts

| Script | Purpose |
|---|---|
| `dedup_lib.py` | Verification kernel that 1:1 reproduces the manuscript detector (full-gene Pearson, per-column centring + unitisation, baseline-aware gate) and reads expression matrices directly from raw GEO `series_matrix.txt.gz` archives — no intermediate artefacts required. |
| `spike_in_benchmark.py` | Constructed ground truth: exact twin copies of real samples under four perturbation families (additive noise, renormalisation bias, gene dropout, gene-value scrambling) across the realistic re-submission envelope and an adversarial envelope. |
| `threshold_robustness.py` | Link recovery on anchor pairs, gate specificity (naive → gated collapse), and negative-control false-positive rate. |
| `corpus_sweep.py` | Full 319-series corpus reproduction of the detector's numerical core plus the two-dimensional threshold-sensitivity scan (τ ∈ {0.95…0.99} × baseline percentile ∈ {p99.0…p99.9}). |
| `rerun_verified.py` | Re-runs the workflows on local CSV corpora with a corrected kernel and emits the **authoritative verified numbers** (an earlier index bug in `threshold_robustness.py`'s pair-link report mis-sliced non-square matrices; see script docstring). |
| `external_known_dups.py` | External positive-control discovery via GEO/Gemma sample-identifier overlap (Gemma's free-text curation notes do not expose cross-accession pair links — see §"External curation coverage" of the manuscript). |
| `gemma_positive_control.py` | Two-stage Gemma query for externally documented duplicate samples. |
| `verify_positive_control.py` | Positive-control spot check. |
| `make_figure_S1.py` | Re-draws the Figure S1 threshold-sensitivity heatmap from `sweep_matrix.csv`. |
| `make_rebuttal_figures.py` | Regenerates the 14 validation figures (F01–F14) and 6 result tables (T01–T06) in `figures/`. |

## Result files (authoritative = `*_verified.csv` / `rerun_verified.py` output)

| File | Content |
|---|---|
| `spike_in_results_verified.csv` | Recall by perturbation family/level at the calibrated gate (T ≈ 0.9886): recall 1.000 for noise ≥ 20 dB, renorm ≤ 0.1; 0.000 for dropout ≥ 0.05, scramble ≥ 0.1, renorm ≥ 0.3. |
| `threshold_robustness_verified.csv` | Anchor link recovery (GSE41258↔GSE68468: peak r 0.982, 26/51 naive pairs retained; GSE12251↔GSE23597: peak r 1.000, 23/258 retained), specificity (GSE29621↔GSE17536 collapse), and 4 negative controls (FPR = 0). |
| `spike_in_results.csv`, `threshold_robustness.csv` | Pre-correction runs (superseded; retained for provenance). |
| `sweep_corpus_table.csv` | Corpus redundancy rate vs correlation threshold at fixed p99.9 (0.95→8.56%, 0.97→8.52%, 0.99→8.02%). |
| `sweep_matrix.csv`, `sweep_matrix_pairs.csv` | Full two-dimensional τ × percentile grid (rates and pair counts). |
| `corpus_sweep_report.json` | Machine-readable sweep output including a reconciliation block against `02_final_reports/dedup_report.json`. |

## Known reconciliation note

The sweep and the final corpus run differ immaterially at the operating point
(τ = 0.97, p99.9): sweep 2,954 gated pairs / 2,168 samples / 1,528 groups vs final run
2,956 / 2,169 / 1,527 (< 0.1%, from candidate pre-screening). **The final corpus run
(`02_final_reports/`) is authoritative** and is what the manuscript reports; the sweep
is used for sensitivity analysis only.

## Requirements

`numpy`, `scipy`, `matplotlib` (see `requirements.txt`). Inputs: raw GEO
`series_matrix.txt.gz` archives (microarray arm) or the converted `<GSE>.csv` layout.
