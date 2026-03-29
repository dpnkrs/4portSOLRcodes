function config = phase1_config(rawDataDir, touchstoneDir, matDir)
%PHASE1_CONFIG Shared configuration for Phase I extraction.

config = struct();
config.raw_data_dir = rawDataDir;
config.touchstone_output_dir = touchstoneDir;
config.mat_output_dir = matDir;

config.band_names = {'0-67', '67-115', '110-170'};
config.band_keys = {'BAND_0_67', 'BAND_67_115', 'BAND_110_170'};

config.line_length_m = 25.44e-6;
config.eps_eff = 3.7245;
config.z0 = 50;
config.zc = 50;
config.alpha_np_per_m = 0;
config.alpha_search_np_per_m = linspace(0, 1500, 11);
config.eps_eff_search = linspace(3.2, 4.2, 11);
config.zc_search = linspace(40, 65, 11);
config.alpha_refine_halfspan_np_per_m = 150;
config.eps_eff_refine_halfspan = 0.12;
config.zc_refine_halfspan = 3;
config.alpha_refine_points = 9;
config.eps_eff_refine_points = 9;
config.zc_refine_points = 9;
config.use_parallel = true;
config.alpha_fit_weights = struct( ...
    'open_short_passivity', 100, ...
    'open_short_target', 20, ...
    'load_magnitude', 2, ...
    'straight_line_match', 20, ...
    'straight_return_loss', 10, ...
    'embedding_passivity', 30, ...
    'alpha_regularization', 5, ...
    'eps_regularization', 2, ...
    'zc_regularization', 2);
config.alpha_regularization_center_np_per_m = 250;
config.alpha_regularization_scale_np_per_m = 250;
config.eps_eff_regularization_center = 3.7245;
config.eps_eff_regularization_scale = 0.2;
config.zc_regularization_center = 50;
config.zc_regularization_scale = 5;
config.embedding_refine_with_sol = true;
config.embedding_refine_weights = struct( ...
    'open_short_passivity_p1', 180, ...
    'open_short_target_p1', 45, ...
    'load_magnitude_p1', 30, ...
    'open_short_passivity_p2', 180, ...
    'open_short_target_p2', 45, ...
    'load_magnitude_p2', 30, ...
    'thru_reconstruct', 140, ...
    'p34_line_match', 50, ...
    'p34_return_loss', 25, ...
    'baseline_s11', 90, ...
    'baseline_s12', 220, ...
    'baseline_s21', 220, ...
    'baseline_s22', 120, ...
    'smoothness_s11', 35, ...
    'smoothness_s12', 120, ...
    'smoothness_s21', 120, ...
    'smoothness_s22', 45, ...
    'reciprocity', 80, ...
    'embedding_passivity_mean', 90, ...
    'embedding_passivity_local', 180, ...
    'denominator_conditioning', 40, ...
    'trust_region_s12s21', 120, ...
    'open_highband_boost', 1.15, ...
    'simulation_ref', 8, ...
    'simulation_raw', 0);
config.simulation_prior_use_load = false;
config.simulation_prior_alignment = true;
config.simulation_prior_alignment_robust = true;
config.embedding_refine_max_iter = 80;
config.embedding_refine_max_fun_evals = 250;
config.embedding_refine_edge_span_hz = 5e9;
config.embedding_refine_edge_boost = 0.9;
config.embedding_refine_overlap_windows_hz = [63e9, 70e9; 112e9, 118e9];
config.embedding_refine_overlap_boost = 0.8;
config.embedding_refine_edge_freeze_hz = 1.5e9;
config.embedding_refine_edge_taper_hz = 3e9;
config.embedding_refine_highband_start_hz = 160e9;
config.embedding_refine_highband_damp_start_hz = 155e9;
config.embedding_refine_highband_baseline_s22_scale = 1.8;
config.embedding_refine_highband_strong_s22_start_hz = 160e9;
config.embedding_refine_highband_strong_s22_scale = 2.4;
config.embedding_refine_highband_open_short_p2_scale = 0.75;
config.embedding_refine_highband_load_p2_scale = 0.85;
config.embedding_refine_denominator_floor = 0.012;
config.embedding_refine_s12s21_trust_radius = 0.04;
config.embedding_refine_transition_windows_hz = [23e9, 27e9; 38e9, 42e9];
config.embedding_refine_transition_p34_scale = 0.25;
config.embedding_refine_transition_trust_scale = 2.0;
config.embedding_refine_transition_baseline_s11_scale = 1.8;
config.embedding_refine_transition_baseline_s22_scale = 1.8;
config.embedding_refine_transition_soft_hz = 1.0e9;
config.embedding_refine_0_67_seconddiff_window_hz = [20e9, 45e9];
config.embedding_refine_0_67_seconddiff_weight_s11 = 120;
config.embedding_refine_0_67_seconddiff_weight_s22 = 140;
config.embedding_refine_0_67_s22_focus_window_hz = [21.5e9, 24.5e9];
config.embedding_refine_0_67_s22_focus_weight = 260;
config.embedding_refine_0_67_disable_sim_window_hz = [20e9, 45e9];
config.embedding_root_jump_threshold = 0.20;
config.embedding_root_blend_floor = 0.25;
config.embedding_root_enable_jump_blend = false;
config.simulation_priors = struct('enabled', false);
config.sol_passivity_enforcement = false;
config.sol_passivity_target = 0.999;
config.sol_passivity_smoothing_window = 7;

config.supported_standard_types = {'short', 'open', 'load'};
config.supported_geometries = {'straight', 'arc', 'diagonal', 'diag'};
end





