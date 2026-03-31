# PROJECT_CONTEXT

## Purpose
This repository contains MATLAB code for mmWave on-wafer calibration workflows:
- legacy 2-port SOLR calibration/validation (`SOLR Cal`)
- active Phase I standard extraction from fabricated-die measurements (`Phase I Extraction`)

The immediate project focus is **Phase I extraction**:
- derive de-embedded `Short/Open/Load` and reciprocal-thru standards from probe-tip-calibrated measurements
- feed these extracted standards into later 4-port SOLR work

---

## Repository Layout

### Top-level folders
- `Phase I Extraction/`  
  Active code for fabricated-die standard extraction.
- `SOLR Cal/`  
  Older 2-port SOLR implementation and validation scripts from pre-fab flow.
- `Measured Data/PhaseI Raw Standards/`  
  Current measured Phase I input data (`.s2p`) for all bands and standards.
- `Backup/Phase I Extraction/`  
  Backup snapshot of prior Phase I scripts.

### Top-level context files
- `PNA_X_mmWave_SOLR_Context.txt`  
  Measurement/procedure context for PNA-X based workflow.
- `Table X Measurement.txt` and `Table X Measurement Revised.txt`  
  Measurement plan and revised implementation mapping.
- `AMPC_arxiv.pdf`, `AMPC_arxiv.txt`  
  Background paper/context.

---

## MATLAB Entry Points

### Active Phase I run
- `Phase I Extraction/main_phase1_extraction.m`  
  Clears workspace/console/figures, adds path, runs extraction.
- `Phase I Extraction/run_phase1_extraction.m`  
  Main driver called by `main_phase1_extraction.m`.

### Legacy 2-port SOLR run
- `SOLR Cal/mainCalAndValid.m`  
  Legacy simulation/validation script for 2-port SOLR.

---

## Phase I Architecture

## 1) Discovery and config
- `phase1_config.m` defines:
  - band list (`0-67`, `67-115`, `110-170`)
  - line model params and search grids
  - refinement weights and options
  - edge/overlap weighting (recently added)
- `discover_phase1_measurements.m`, `parse_phase1_filename.m`
  discover + map raw files by band and standard key.

## 2) Line model fitting
- `extract_phase1_band.m` calls `fit_line_alpha_for_band(...)`
- Coarse-to-fine search over:
  - `alpha`
  - `eps_eff`
  - `Zc`
- Objective terms include:
  - SOL passivity/target/load penalties
  - straight-thru line match / return-loss
  - regularization terms
  - embedding passivity penalty

## 3) Embedding extraction
- `extract_embedding_from_thru.m`
  extracts baseline embedding `E(f)` from `P1P2` straight thru:
  - `THRU = E * L * E`
  - matrix half-factor via `expm(0.5*logm(...))`
  - sign continuity across frequency

## 4) Embedding refinement (current)
- `refine_embedding_with_sol.m`
- Current mode is **weighted Stage A**:
  - refine `S11/S22` only
  - keep `S12/S21` fixed to baseline
  - include P1/P2 SOL constraints
  - include P12 thru reconstruction + P34 verification penalties
  - include local passivity penalty
  - newly added frequency weighting near band edges and overlap windows

## 5) De-embedding and export
- One-port: `deembed_oneport_standard.m`
- Two-port: `deembed_twport_standard.m`
- Touchstone readers/writers:
  - `read_touchstone_file.m`
  - `write_touchstone_1port.m`
  - `write_touchstone_2port.m`
- Outputs saved per band to:
  - `Phase I Extraction/outputs/touchstone/`
  - `Phase I Extraction/outputs/mat/`

## 6) Diagnostics and plots
- `plot_phase1_band_results.m`
  (per-band SOL/thru + fit/embedding/refinement diagnostics)
- `plot_phase1_stitched_results.m`
  stitched full-band curves using overlap alignment/blending
- Notes written to `outputs/notes/`:
  - `PhaseI_ObjectiveBreakdown_<band>.txt`
  - `PhaseI_EmbeddingComparison_<band>.txt`
  - `PhaseI_RefinementSummary_<band>.txt`

---

## Data Flow (Phase I)

1. Read raw `.s2p` standards from `Measured Data/PhaseI Raw Standards/`
   (with backward-compatible fallback to the older folder names).
2. Fit line model parameters per band (`alpha/eps_eff/Zc`).
3. Extract baseline embedding `E(f)` from `THRU_P1P2_STRAIGHT`.
4. Optionally refine `E(f)` using measured P1/P2 SOL and thru constraints.
5. De-embed:
   - `SHORT/OPEN/LOAD` (P1/P2) -> extracted one-port standards
   - reciprocal thrus -> extracted two-port standards
6. Write outputs (`.s1p/.s2p` + `.mat`) and diagnostics.

---

## Output Artifacts

### Primary extracted artifacts
- `outputs/touchstone/EXTRACTED_*.s1p`
- `outputs/touchstone/EXTRACTED_*.s2p`
- `outputs/mat/EXTRACTED_*.mat`
- `outputs/mat/EXTRACTED_EMBEDDING_<band>.mat`
- `outputs/mat/PHASE1_RESULTS_ALL.mat`

### Diagnostic artifacts
- `outputs/figures/PhaseI_*.jpg` and `.fig`
- `outputs/notes/PhaseI_*.txt`

---

## Current Status Snapshot (important for resume)

The repo is mid-iteration on embedding-level refinement.

Recent path:
- Stage B (allowing `S12/S21` movement) was tried and rolled back due instability/discontinuities.
- Stage C (global control-point refinement) was tried and rolled back due poor behavior.
- Current code has been moved to a **weighted Stage A** approach (S11/S22-only refinement), adding edge/overlap emphasis.

Practical implications:
- existing files in `outputs/*` may reflect runs from earlier variants
- after moving to new PC, run a fresh extraction before judging current behavior

### Outstanding technical issue
- De-embedded `open/short` physicality remains challenging, especially at high frequencies.
- Transmission-only approximation diagnostic (`Gamma_meas/(S12*S21)`) remains far from full model in many cases, indicating reflective terms are still dominant.

---

## Resume Checklist on New PC

1. Restore repository to same path (or update MATLAB working directory).
2. Ensure MATLAB toolboxes are available:
   - RF Toolbox (or compatible reader path; Phase I uses custom reader but legacy SOLR uses RF toolbox helpers)
   - Parallel Computing Toolbox (optional; code falls back without parallel)
3. Confirm raw files exist under:
   - `Measured Data/PhaseI Raw Standards/`
4. Run:
   - `Phase I Extraction/main_phase1_extraction.m`
5. Inspect:
   - `outputs/notes/PhaseI_RefinementSummary_*.txt`
   - `outputs/figures/PhaseI_SOL_Magnitude_*.jpg`
   - `outputs/figures/PhaseI_THRU_S21_Magnitude_*.jpg`
6. If continuing optimization work, start from current `refine_embedding_with_sol.m` and `phase1_config.m`.

---

## Key Files for Future Codex Session

If Codex is reopened on a new machine, provide these first:
- `PROJECT_CONTEXT.md` (this file)
- `Phase I Extraction/phase1_config.m`
- `Phase I Extraction/refine_embedding_with_sol.m`
- `Phase I Extraction/extract_phase1_band.m`
- `Phase I Extraction/plot_phase1_band_results.m`
- latest `outputs/notes/PhaseI_RefinementSummary_*.txt`

And optionally for historical rationale:
- `Phase I Extraction/Implementation Notes.txt`
- `Phase I Extraction/Corrections and Adjustments Notes.txt`

---

## Notes on Backward Compatibility

- `Backup/Phase I Extraction/` preserves prior script versions for diff/recovery.
- `SOLR Cal/` remains functional as legacy reference and has not been integrated yet with current extracted Phase I outputs.
