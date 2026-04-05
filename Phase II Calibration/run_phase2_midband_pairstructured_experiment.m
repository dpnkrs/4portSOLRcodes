function run_phase2_midband_pairstructured_experiment()
%RUN_PHASE2_MIDBAND_PAIRSTRUCTURED_EXPERIMENT
% Controlled experiment: use downstream straight references as the absolute midband anchor
% and preserve each orthogonal pair's own excess-loss offset relative to that straight baseline.

baseConfig = phase2_config();
baseResultsPath = fullfile(baseConfig.mat_dir, 'PHASE2_RESULTS_ALL.mat');
if ~exist(baseResultsPath, 'file')
    error('Baseline Phase II results not found: %s', baseResultsPath);
end
loaded = load(baseResultsPath, 'results');
baselineResults = loaded.results;

[pairScales, straightModel] = derive_pairstructured_midband_scales(baselineResults);

expRoot = fullfile(baseConfig.interim_dir, 'pairStructuredExperiment');
expConfig = baseConfig;
expConfig.output_dir = fullfile(expRoot, 'outputs');
expConfig.touchstone_dir = fullfile(expConfig.output_dir, 'touchstone');
expConfig.mat_dir = fullfile(expConfig.output_dir, 'mat');
expConfig.figure_dir = fullfile(expConfig.output_dir, 'figures');
expConfig.note_dir = fullfile(expConfig.output_dir, 'notes');
expConfig.interim_dir = fullfile(expRoot, 'interim');
expConfig.interim_mat_dir = fullfile(expConfig.interim_dir, 'mat');
expConfig.interim_figure_dir = fullfile(expConfig.interim_dir, 'figures');
expConfig.interim_note_dir = fullfile(expConfig.interim_dir, 'notes');
ensure_output_dirs_local(expConfig);
apply_plot_defaults_local();

experimentResults = baselineResults;
bandIdx = find(strcmpi(cellfun(@(b) b.band, baselineResults.bands, 'UniformOutput', false), '67-115'), 1, 'first');
if isempty(bandIdx)
    error('67-115 band not found in baseline results.');
end

baseBand = baselineResults.bands{bandIdx};
modBand = baseBand;
for idxPair = 1:numel(modBand.reference_pair_results)
    pr = modBand.reference_pair_results(idxPair);
    key = sprintf('%s_%s', upper(pr.pair), upper(pr.geometry));
    if isfield(pairScales, key)
        scaleInfo = pairScales.(key);
        modBand.reference_pair_results(idxPair).error_terms.t21 = modBand.reference_pair_results(idxPair).error_terms.t21 .* scaleInfo.scale_mag;
        modBand.reference_pair_results(idxPair).error_terms.t12 = modBand.reference_pair_results(idxPair).error_terms.t12 .* scaleInfo.scale_mag;
        modBand.reference_pair_results(idxPair).pairstructured_experiment = scaleInfo;
    else
        modBand.reference_pair_results(idxPair).pairstructured_experiment = struct('applied', false);
    end
    modBand.reference_pair_results(idxPair).reference_corrected = correct_network_from_error_terms_local( ...
        modBand.reference_pair_results(idxPair).reference_switch_corrected, ...
        modBand.reference_pair_results(idxPair).error_terms, ...
        expConfig.cal_den_floor);
    modBand.reference_pair_results(idxPair).reference_metrics = summarize_pair_metrics_local( ...
        modBand.reference_pair_results(idxPair).reference_corrected);
end

for idxHold = 1:numel(modBand.holdout_results)
    holdCfg = expConfig.holdout_pairs(idxHold);
    referenceIdx = find(strcmpi({modBand.reference_pair_results.pair}, holdCfg.pair), 1, 'first');
    holdNet = read_touchstone_nport(modBand.holdout_results(idxHold).file);
    modBand.holdout_results(idxHold) = correct_pair_network_local( ...
        expConfig, holdNet, modBand.reference_pair_results(referenceIdx), holdCfg.pair, upper(holdCfg.geometry), true);
end

pairBlocks = cell(1, numel(modBand.reference_pair_results));
for idxPair = 1:numel(modBand.reference_pair_results)
    pairBlocks{idxPair} = correct_pair_network_local( ...
        expConfig, modBand.dut.raw, modBand.reference_pair_results(idxPair), ...
        modBand.reference_pair_results(idxPair).pair, modBand.reference_pair_results(idxPair).geometry, false);
end
pairBlocks = [pairBlocks{:}];
modBand.dut.pair_corrected_blocks = pairBlocks;
modBand.dut.corrected = assemble_corrected_4port_from_pairs(expConfig, pairBlocks);
modBand.dut.summary = summarize_dut_metrics_local(modBand.dut.corrected.S);

experimentResults.bands{bandIdx} = modBand;
experimentResults.stitched = build_stitched_outputs_local(experimentResults.bands, expConfig);

save_experiment_outputs(expRoot, expConfig, baselineResults, experimentResults, pairScales, straightModel);

bundle = struct();
bundle.baseline = baselineResults;
bundle.experiment = experimentResults;
bundle.pair_scales = pairScales;
bundle.straight_model = straightModel;
save(fullfile(expRoot, 'PHASE2_MIDBAND_PAIRSTRUCTURED_EXPERIMENT.mat'), 'bundle');
end

function [pairScales, straightModel] = derive_pairstructured_midband_scales(results)
bandIdx = find(strcmpi(cellfun(@(b) b.band, results.bands, 'UniformOutput', false), '67-115'), 1, 'first');
bandResult = results.bands{bandIdx};
pairScales = struct();

straightMask = strcmpi({bandResult.reference_pair_results.geometry}, 'STRAIGHT');
straightPairs = bandResult.reference_pair_results(straightMask);
if numel(straightPairs) ~= 2
    error('Expected exactly two straight reference pairs in 67-115.');
end

straightCorr21 = [];
straightCorr12 = [];
straightSw21 = [];
straightSw12 = [];
for idx = 1:numel(straightPairs)
    straightCorr21(:, idx) = mag_db(squeeze(straightPairs(idx).reference_corrected(2,1,:))); %#ok<AGROW>
    straightCorr12(:, idx) = mag_db(squeeze(straightPairs(idx).reference_corrected(1,2,:))); %#ok<AGROW>
    straightSw21(:, idx) = mag_db(squeeze(straightPairs(idx).reference_switch_corrected(2,1,:))); %#ok<AGROW>
    straightSw12(:, idx) = mag_db(squeeze(straightPairs(idx).reference_switch_corrected(1,2,:))); %#ok<AGROW>
end

straightModel = struct();
straightModel.target21_db = mean(straightCorr21, 2, 'omitnan');
straightModel.target12_db = mean(straightCorr12, 2, 'omitnan');
straightModel.switch21_db = mean(straightSw21, 2, 'omitnan');
straightModel.switch12_db = mean(straightSw12, 2, 'omitnan');
straightModel.freq = bandResult.freq(:);

for idx = 1:numel(bandResult.reference_pair_results)
    pr = bandResult.reference_pair_results(idx);
    if strcmpi(pr.geometry, 'STRAIGHT')
        continue;
    end
    key = sprintf('%s_%s', upper(pr.pair), upper(pr.geometry));
    arcSw21Db = mag_db(squeeze(pr.reference_switch_corrected(2,1,:)));
    arcSw12Db = mag_db(squeeze(pr.reference_switch_corrected(1,2,:)));
    arcCorr21Db = mag_db(squeeze(pr.reference_corrected(2,1,:)));
    arcCorr12Db = mag_db(squeeze(pr.reference_corrected(1,2,:)));

    desired21Db = straightModel.target21_db + (arcSw21Db - straightModel.switch21_db);
    desired12Db = straightModel.target12_db + (arcSw12Db - straightModel.switch12_db);

    delta21Db = arcCorr21Db - desired21Db;
    delta12Db = arcCorr12Db - desired12Db;
    deltaDb = median([delta21Db; delta12Db], 'omitnan');
    scaleMag = 10^(deltaDb / 20);

    pairScales.(key) = struct( ...
        'applied', true, ...
        'pair', pr.pair, ...
        'geometry', pr.geometry, ...
        'scale_db', deltaDb, ...
        'scale_mag', scaleMag, ...
        'baseline_corrected_s21_db', median(arcCorr21Db, 'omitnan'), ...
        'baseline_corrected_s12_db', median(arcCorr12Db, 'omitnan'), ...
        'desired_s21_db', median(desired21Db, 'omitnan'), ...
        'desired_s12_db', median(desired12Db, 'omitnan'));
    end
end

function corrected = correct_pair_network_local(config, network, pairModel, pairKey, geometry, isHoldout)
ports = phase2_pair_key_to_ports(pairKey);
Sraw = phase2_extract_pair_submatrix(network.S, pairKey);
switchTerms = pairModel.switch_terms;

if numel(network.freq) ~= numel(pairModel.freq) || any(abs(network.freq(:) - pairModel.freq(:)) > 1)
    switchTerms = struct();
    switchTerms.gamma_forward = interp1(pairModel.freq, pairModel.switch_terms.gamma_forward, network.freq, 'linear', 'extrap');
    switchTerms.gamma_reverse = interp1(pairModel.freq, pairModel.switch_terms.gamma_reverse, network.freq, 'linear', 'extrap');
    errorTerms = interpolate_error_terms_local(pairModel.error_terms, pairModel.freq, network.freq);
else
    errorTerms = pairModel.error_terms;
end

Sswitch = apply_switch_correction_2port(Sraw, switchTerms.gamma_forward, switchTerms.gamma_reverse, config.switch_den_floor);
Scorrected = correct_network_from_error_terms_local(Sswitch, errorTerms, config.cal_den_floor);

metrics = struct();
metrics.reciprocity_mismatch = abs(squeeze(Scorrected(2, 1, :)) - squeeze(Scorrected(1, 2, :)));
metrics.return_loss = min(-20 * log10(max(abs(squeeze(Scorrected(1, 1, :))), 1e-12)), ...
    -20 * log10(max(abs(squeeze(Scorrected(2, 2, :))), 1e-12)));
metrics.insertion_loss_db = -20 * log10(max(abs(squeeze(Scorrected(2, 1, :))), 1e-12));
metrics.phase_s21_deg = unwrap(angle(squeeze(Scorrected(2, 1, :)))) * 180 / pi;
metrics.phase_s12_deg = unwrap(angle(squeeze(Scorrected(1, 2, :)))) * 180 / pi;

corrected = struct();
corrected.file = network.file;
corrected.pair = pairKey;
corrected.geometry = geometry;
corrected.ports = ports;
corrected.freq = network.freq(:);
corrected.S_raw = Sraw;
corrected.S_switch = Sswitch;
corrected.S_corrected = Scorrected;
corrected.corrected = Scorrected;
corrected.metrics = metrics;
corrected.is_holdout = isHoldout;
end

function termsInterp = interpolate_error_terms_local(errorTerms, freqIn, freqOut)
fields = fieldnames(errorTerms);
termsInterp = struct();
for idx = 1:numel(fields)
    values = errorTerms.(fields{idx});
    realPart = interp1(freqIn, real(values), freqOut, 'linear', 'extrap');
    imagPart = interp1(freqIn, imag(values), freqOut, 'linear', 'extrap');
    termsInterp.(fields{idx}) = complex(realPart, imagPart);
end
end

function corrected = correct_network_from_error_terms_local(Sraw, errorTerms, denominatorFloor)
corrected = nan(size(Sraw));
for idxFreq = 1:size(Sraw, 3)
    corrected(:, :, idxFreq) = apply_error_correction_8term_local(Sraw(:, :, idxFreq), errorTerms, idxFreq, denominatorFloor);
end
end

function metrics = summarize_pair_metrics_local(S)
s21 = squeeze(S(2, 1, :));
s12 = squeeze(S(1, 2, :));
s11 = squeeze(S(1, 1, :));
s22 = squeeze(S(2, 2, :));
metrics = struct();
metrics.mean_recip_mismatch = mean(abs(s21 - s12), 'omitnan');
metrics.max_recip_mismatch = max(abs(s21 - s12), [], 'omitnan');
metrics.median_insertion_loss_db = median(-20 * log10(max(abs(s21), 1e-12)), 'omitnan');
metrics.max_return_loss_db = max(-20 * log10(max([abs(s11); abs(s22)], 1e-12)), [], 'omitnan');
end

function metrics = summarize_dut_metrics_local(S)
nFreq = size(S, 3);
offDiag = [];
recip = [];
for row = 1:4
    for col = 1:4
        if row == col
            continue;
        end
        offDiag = [offDiag; squeeze(S(row, col, :)).']; %#ok<AGROW>
        if row < col
            recip = [recip; abs(squeeze(S(row, col, :)).' - squeeze(S(col, row, :)).')]; %#ok<AGROW>
        end
    end
end

diagVals = zeros(4, nFreq);
for idxPort = 1:4
    diagVals(idxPort, :) = squeeze(S(idxPort, idxPort, :)).';
end

port1 = squeeze(S(2:4, 1, :));
medianPort1 = median(-20 * log10(max(abs(port1), 1e-12)), 2, 'omitnan');
sortedMedian = sort(medianPort1);

metrics = struct();
metrics.median_return_loss_db = median(-20 * log10(max(abs(diagVals(:)), 1e-12)), 'omitnan');
metrics.median_reciprocity_mismatch = median(recip(:), 'omitnan');
metrics.max_reciprocity_mismatch = max(recip(:), [], 'omitnan');
metrics.isolation_db = max(median(-20 * log10(max(abs(offDiag), 1e-12)), 2, 'omitnan'));
if numel(sortedMedian) >= 2
    metrics.coupling_balance_db = abs(sortedMedian(1) - sortedMedian(2));
else
    metrics.coupling_balance_db = NaN;
end
magDb = 20 * log10(max(abs(offDiag), 1e-12));
metrics.smoothness_rms_db = rms(diff(magDb, 2, 2), 'all');
end

function stitched = build_stitched_outputs_local(bandResults, config)
stitched = struct();
dutNetworks = cell(1, numel(bandResults));
for idxBand = 1:numel(bandResults)
    dutNetworks{idxBand} = bandResults{idxBand}.dut.corrected;
end
stitched.dut = stitch_phase2_networks(dutNetworks, config.stitch_interp_method);

for idxHold = 1:numel(config.holdout_pairs)
    pairKey = config.holdout_pairs(idxHold).pair;
    geomKey = upper(config.holdout_pairs(idxHold).geometry);
    holdNetworks = cell(1, numel(bandResults));
    for idxBand = 1:numel(bandResults)
        hit = find(strcmpi({bandResults{idxBand}.holdout_results.pair}, pairKey), 1, 'first');
        holdNetworks{idxBand} = struct( ...
            'freq', bandResults{idxBand}.holdout_results(hit).freq, ...
            'S', bandResults{idxBand}.holdout_results(hit).S_corrected, ...
            'z0', config.z0);
    end
    stitched.holdout.(pairKey).(geomKey) = stitch_phase2_networks(holdNetworks, config.stitch_interp_method);
end
end

function save_experiment_outputs(expRoot, config, baselineResults, experimentResults, pairScales, straightModel)
write_experiment_summary(expRoot, baselineResults, experimentResults, pairScales, straightModel);
plot_experiment_comparison(expRoot, baselineResults, experimentResults);

save(fullfile(expRoot, 'outputs', 'mat', 'PHASE2_PAIRSTRUCTURED_RESULTS.mat'), 'experimentResults', 'pairScales', 'straightModel');
write_touchstone_nport(fullfile(expRoot, 'outputs', 'touchstone', 'CORRECTED_LangeCoupler_STITCHED.s4p'), ...
    experimentResults.stitched.dut.freq, experimentResults.stitched.dut.S, config.z0);
end

function write_experiment_summary(expRoot, baselineResults, experimentResults, pairScales, straightModel)
lines = {};
lines{end+1} = 'Phase II 67-115 pair-structured experiment';
lines{end+1} = '=========================================';
lines{end+1} = '';
lines{end+1} = sprintf('Straight-model median corrected S21 / S12 (dB): %.3f / %.3f', ...
    median(straightModel.target21_db, 'omitnan'), median(straightModel.target12_db, 'omitnan'));
lines{end+1} = sprintf('Straight-model median switch S21 / S12 (dB): %.3f / %.3f', ...
    median(straightModel.switch21_db, 'omitnan'), median(straightModel.switch12_db, 'omitnan'));
lines{end+1} = '';

scaleKeys = fieldnames(pairScales);
lines{end+1} = 'Applied ARC pair scales in 67-115';
for idx = 1:numel(scaleKeys)
    s = pairScales.(scaleKeys{idx});
    lines{end+1} = sprintf('  %s %s scale: %.3f dB (mag %.4f) | baseline S21 %.3f dB -> desired %.3f dB', ...
        s.pair, s.geometry, s.scale_db, s.scale_mag, s.baseline_corrected_s21_db, s.desired_s21_db);
end
lines{end+1} = '';

bandIdx = find(strcmpi(cellfun(@(b) b.band, baselineResults.bands, 'UniformOutput', false), '67-115'), 1, 'first');
baseBand = baselineResults.bands{bandIdx};
expBand = experimentResults.bands{bandIdx};
lines{end+1} = '67-115 reference pairs';
for idxPair = 1:numel(baseBand.reference_pair_results)
    bpr = baseBand.reference_pair_results(idxPair);
    epr = expBand.reference_pair_results(idxPair);
    lines{end+1} = sprintf('  %s %s baseline / experiment corrected |S21|: %.3f / %.3f dB', ...
        bpr.pair, bpr.geometry, median(mag_db(squeeze(bpr.reference_corrected(2,1,:))), 'omitnan'), ...
        median(mag_db(squeeze(epr.reference_corrected(2,1,:))), 'omitnan'));
end
lines{end+1} = '';
lines{end+1} = '67-115 holdouts';
for idxHold = 1:numel(baseBand.holdout_results)
    bhr = baseBand.holdout_results(idxHold);
    ehr = expBand.holdout_results(idxHold);
    lines{end+1} = sprintf('  %s %s baseline / experiment corrected |S21|: %.3f / %.3f dB', ...
        bhr.pair, bhr.geometry, median(mag_db(squeeze(bhr.S_corrected(2,1,:))), 'omitnan'), ...
        median(mag_db(squeeze(ehr.S_corrected(2,1,:))), 'omitnan'));
end
lines{end+1} = '';
lines{end+1} = '67-115 DUT';
lines{end+1} = sprintf('  Baseline / experiment median S21: %.3f / %.3f dB', ...
    median(mag_db(squeeze(baseBand.dut.corrected.S(2,1,:))), 'omitnan'), ...
    median(mag_db(squeeze(expBand.dut.corrected.S(2,1,:))), 'omitnan'));
lines{end+1} = sprintf('  Baseline / experiment median S31: %.3f / %.3f dB', ...
    median(mag_db(squeeze(baseBand.dut.corrected.S(3,1,:))), 'omitnan'), ...
    median(mag_db(squeeze(expBand.dut.corrected.S(3,1,:))), 'omitnan'));
lines{end+1} = sprintf('  Baseline / experiment median S41: %.3f / %.3f dB', ...
    median(mag_db(squeeze(baseBand.dut.corrected.S(4,1,:))), 'omitnan'), ...
    median(mag_db(squeeze(expBand.dut.corrected.S(4,1,:))), 'omitnan'));

write_text_file_local(fullfile(expRoot, 'PhaseII_PairStructuredExperiment_Summary.txt'), lines);
end

function plot_experiment_comparison(expRoot, baselineResults, experimentResults)
fig = figure('Visible', 'on', 'Color', 'w');
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
title(ax1, 'Stitched Lange coupler transmission: baseline vs pair-structured experiment');
xlabel(ax1, 'Frequency (GHz)'); ylabel(ax1, '|S_{x1}| (dB)');

ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on');
title(ax2, 'Stitched DIAG hold-outs: baseline vs pair-structured experiment');
xlabel(ax2, 'Frequency (GHz)'); ylabel(ax2, '|S_{21}| (dB)');

freqBase = baselineResults.stitched.dut.freq / 1e9;
freqExp = experimentResults.stitched.dut.freq / 1e9;
colors = lines(6);
labels = {'S21', 'S31', 'S41'};
for idxPort = 2:4
    plot(ax1, freqBase, mag_db(squeeze(baselineResults.stitched.dut.S(idxPort,1,:))), '--', ...
        'Color', colors(idxPort - 1, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax1, freqExp, mag_db(squeeze(experimentResults.stitched.dut.S(idxPort,1,:))), '-', ...
        'Color', colors(idxPort - 1, :), 'LineWidth', 1.4, 'DisplayName', labels{idxPort - 1});
end
legend(ax1, 'Location', 'best');

holdNames = fieldnames(experimentResults.stitched.holdout);
for idxHold = 1:numel(holdNames)
    pairKey = holdNames{idxHold};
    baseNet = baselineResults.stitched.holdout.(pairKey).DIAG;
    expNet = experimentResults.stitched.holdout.(pairKey).DIAG;
    c = colors(idxHold + 2, :);
    plot(ax2, baseNet.freq / 1e9, mag_db(squeeze(baseNet.S(2,1,:))), '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax2, expNet.freq / 1e9, mag_db(squeeze(expNet.S(2,1,:))), '-', 'Color', c, 'LineWidth', 1.4, 'DisplayName', pairKey);
end
legend(ax2, 'Location', 'bestoutside');

save_debug_figure_local(fig, fullfile(expRoot, 'outputs', 'figures', 'PhaseII_PairStructuredExperiment_Comparison'));
end

function ensure_output_dirs_local(config)
dirs = {config.output_dir, config.touchstone_dir, config.mat_dir, config.figure_dir, config.note_dir, ...
    config.interim_dir, config.interim_mat_dir, config.interim_figure_dir, config.interim_note_dir};
for idx = 1:numel(dirs)
    if ~exist(dirs{idx}, 'dir')
        mkdir(dirs{idx});
    end
end
end

function apply_plot_defaults_local()
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultTextColor', 'k');
set(groot, 'defaultAxesColorOrder', lines(12));
set(groot, 'defaultAxesLineStyleOrder', '-');
set(groot, 'defaultLineLineWidth', 1.2);
end

function values = mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function save_debug_figure_local(fig, basePath)
folder = fileparts(basePath);
if ~exist(folder, 'dir')
    mkdir(folder);
end
set(fig, 'Color', 'w');
ax = findall(fig, 'type', 'axes');
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
drawnow;
exportgraphics(fig, [basePath, '.jpg'], 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, [basePath, '.pdf'], 'BackgroundColor', 'white', 'ContentType', 'vector');
savefig(fig, [basePath, '.fig']);
close(fig);
end

function write_text_file_local(filePath, lines)
folder = fileparts(filePath);
if ~exist(folder, 'dir')
    mkdir(folder);
end
fid = fopen(filePath, 'w');
cleanup = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end
