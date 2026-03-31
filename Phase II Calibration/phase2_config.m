function config = phase2_config()
%PHASE2_CONFIG Central configuration for downstream Phase II calibration.

scriptDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(scriptDir);
measuredRoot = fullfile(rootDir, 'Measured Data');

config = struct();
config.script_dir = scriptDir;
config.root_dir = rootDir;
config.measured_root = measuredRoot;

config.phase1_output_mat_dir = fullfile(rootDir, 'Phase I Extraction', 'outputs', 'mat');
config.phase2_raw_dir = fullfile(measuredRoot, 'PhaseII Raw Standards');
config.phase2_switch_dir = fullfile(measuredRoot, 'PhaseII Switch Terms');
config.phase2_dut_dir = fullfile(measuredRoot, 'PhaseII DUTs');

config.output_dir = fullfile(scriptDir, 'outputs');
config.touchstone_dir = fullfile(config.output_dir, 'touchstone');
config.mat_dir = fullfile(config.output_dir, 'mat');
config.figure_dir = fullfile(config.output_dir, 'figures');
config.note_dir = fullfile(config.output_dir, 'notes');
config.interim_dir = fullfile(scriptDir, 'interimOutputs');
config.interim_mat_dir = fullfile(config.interim_dir, 'mat');
config.interim_figure_dir = fullfile(config.interim_dir, 'figures');
config.interim_note_dir = fullfile(config.interim_dir, 'notes');

config.bands = {'0-67', '67-115', '110-170'};
config.band_keys = matlab.lang.makeValidName(strrep(config.bands, '-', '_'));

config.port_labels = {'P1', 'P2', 'P3', 'P4'};
config.z0 = 50;
config.switch_den_floor = 1e-6;
config.cal_den_floor = 1e-9;
config.switch_interp_method = 'linear';
config.stitch_interp_method = 'pchip';
config.figure_visible = 'on';
config.switch_ratio_mode = 'a_over_b';
config.switch_debug_pairs = struct( ...
    'pair', {'P1P2', 'P3P4'}, ...
    'geometry', {'STRAIGHT', 'STRAIGHT'});

config.reference_pairs = struct( ...
    'pair', {'P1P2', 'P3P4', 'P1P3', 'P1P4', 'P2P3', 'P2P4'}, ...
    'geometry', {'STRAIGHT', 'STRAIGHT', 'ARC', 'ARC', 'ARC', 'ARC'});

config.holdout_pairs = struct( ...
    'pair', {'P1P3', 'P1P4', 'P2P3', 'P2P4'}, ...
    'geometry', {'DIAG', 'DIAG', 'DIAG', 'DIAG'});

config.standard_names = {'Open', 'Short', 'Load'};
config.reference_standard_files = struct( ...
    'Open', fullfile(config.phase1_output_mat_dir, 'FINALFIT_OPEN_STITCHED_FINAL.mat'), ...
    'Short', fullfile(config.phase1_output_mat_dir, 'FINALFIT_SHORT_STITCHED_FINAL.mat'), ...
    'Load', fullfile(config.phase1_output_mat_dir, 'FINALFIT_LOAD_STITCHED_FINAL.mat'));

config.diagonal_weight_straight = 2;
config.diagonal_weight_orthogonal = 1;
config.walkthrough_filename = fullfile(config.script_dir, 'PhaseII Walkthrough.txt');
end
