Phase I Extraction

This folder contains the MATLAB implementation for Phase I calibration
standard extraction from probe-tip-calibrated measurements.

Entry point:
    run_phase1_extraction

Default raw-data location:
    ..\Measured Data\PhaseI Raw Standards

Outputs:
    outputs\touchstone
    outputs\mat
    outputs\figures

Current assumptions:
    - one shared embedding network per band
    - straight thru model uses:
          length = 25.44 um
          Zc fit per band
          alpha and eps_eff fit per band from a soft objective
    - embedding E(f) is extracted from P1-P2 straight thru
    - P3-P4 straight thru is used as a verification constraint
    - reflective terms of E(f) can be softly refined using measured
      Short/Open/Load, while keeping transmission behavior anchored by the
      thru-based extraction
    - optimization uses a coarse-to-fine search
    - parfor is used when Parallel Computing Toolbox is available
    - extracted reciprocal standards are saved as full 2-port S-parameters
    - raw SOL files are .s2p with the driven-port reflection selected as:
          P1/P3 -> S11
          P2/P4 -> S22

Supported filename styles:
    Short_P1_0-67.s2p
    Open_P1_67-115.s2p
    Load_P1_110-170.s2p
    THRU_P1P2_STRAIGHT_0-67.s2p
    Thru_Straight_P1P2_0-67.s2p

Backward-compatible raw-data fallbacks:
    ..\Measured Data\Raw Standards
    ..\Measured Data\RawStandards
