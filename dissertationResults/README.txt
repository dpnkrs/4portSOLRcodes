dissertationResults

This module converts the saved Phase I extraction outputs into dissertation-facing
metrics, tables, and section-ready notes. It does not rerun the extraction itself.

Entry point:
    main_dissertation_results

Inputs:
    ..\Phase I Extraction\outputs\mat\PHASE1_RESULTS_ALL.mat
    ..\Phase I Extraction\outputs\mat\FINALFIT_OPEN_STITCHED_FINAL.mat
    ..\Phase I Extraction\outputs\mat\FINALFIT_SHORT_STITCHED_FINAL.mat
    ..\Phase I Extraction\outputs\mat\FINALFIT_LOAD_STITCHED_FINAL.mat

Outputs:
    outputs\mat
    outputs\tables
    outputs\notes

Primary outputs:
    DissertationResultsSummary.mat
    Table5_1_Metrics.csv
    Table5_2_AdoptedStandards.csv
    BandwiseStandardMetrics.csv
    notes\Section5_1_2_EvaluationFramework.txt
    notes\Section5_1_3_OpenShort.txt
    notes\Section5_1_4_Load.txt
    notes\Section5_1_5_AdoptedStandards.txt
    notes\WriteupSuggestions.txt

Interpretation choices:
    - accepted bandwise standards remain the primary calibration standards
    - stitched final fits are treated as reporting/downstream smooth references
    - direct model closure residual is retained as an algebraic sanity check only
    - held-out redundant residual is computed as the more meaningful validation metric
