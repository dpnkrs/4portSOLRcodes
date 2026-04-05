function run_phase2_midband_normalization_experiment()
%RUN_PHASE2_MIDBAND_NORMALIZATION_EXPERIMENT Test 67-115 reciprocal-target normalization.

baseConfig = phase2_config();
baseResultsPath = fullfile(baseConfig.mat_dir, 'PHASE2_RESULTS_ALL.mat');
if ~exist(baseResultsPath, 'file')
    error('Baseline Phase II results not found: %s', baseResultsPath);
end
loaded = load(baseResultsPath, 'results');
baselineResults = loaded.results;

ratio = derive_midband_normalization_ratio(baselineResults);

expRoot = fullfile(baseConfig.interim_dir, 'midbandNormalizationExperiment');
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
expConfig.transmission_anchor.enabled = true;
expConfig.transmission_anchor.bands = {'67-115'};
expConfig.reciprocal_target_normalization.enabled = true;
expConfig.reciprocal_target_normalization.band = '67-115';
expConfig.reciprocal_target_normalization.ratio = ratio;

experimentResults = run_phase2_solr(expConfig);
write_experiment_summary(expRoot, baselineResults, experimentResults, ratio);
plot_experiment_comparison(expRoot, baselineResults, experimentResults);

bundle = struct();
bundle.baseline = baselineResults;
bundle.experiment = experimentResults;
bundle.normalization_ratio = ratio;
save(fullfile(expRoot, 'PHASE2_MIDBAND_NORMALIZATION_EXPERIMENT.mat'), 'bundle');
end

function ratio = derive_midband_normalization_ratio(results)
bandIdx = find(strcmpi(cellfun(@(b) b.band, results.bands, 'UniformOutput', false), '67-115'), 1, 'first');
if isempty(bandIdx)
    error('67-115 band not found in baseline results');
end
bandResult = results.bands{bandIdx};
straightMask = strcmpi({bandResult.reference_pair_results.geometry}, 'STRAIGHT');
straightPairs = bandResult.reference_pair_results(straightMask);
allRatios = [];
for idx = 1:numel(straightPairs)
    pr = straightPairs(idx);
    sw21 = squeeze(pr.reference_switch_corrected(2,1,:));
    sw12 = squeeze(pr.reference_switch_corrected(1,2,:));
    tgt21 = squeeze(pr.reference_phase1_target.S(2,1,:));
    tgt12 = squeeze(pr.reference_phase1_target.S(1,2,:));
    allRatios = [allRatios; sw21 ./ tgt21; sw12 ./ tgt12]; %#ok<AGROW>
end
allRatios = allRatios(isfinite(real(allRatios)) & isfinite(imag(allRatios)) & abs(allRatios) > 1e-12);
if isempty(allRatios)
    error('No valid straight-pair ratios found for normalization');
end
magDb = 20 * log10(abs(allRatios));
qLo = quantile(magDb, 0.15);
qHi = quantile(magDb, 0.85);
keep = magDb >= qLo & magDb <= qHi;
allRatios = allRatios(keep);
magDb = magDb(keep);
ratioDb = median(magDb, 'omitnan');
ratioPhase = angle(mean(exp(1j * angle(allRatios)), 'omitnan'));
ratio = 10^(ratioDb/20) * exp(1j * ratioPhase);
end

function write_experiment_summary(expRoot, baselineResults, experimentResults, ratio)
lines = {};
lines{end+1} = 'Phase II 67-115 reciprocal-target normalization experiment';
lines{end+1} = '====================================================';
lines{end+1} = '';
lines{end+1} = sprintf('Applied normalization ratio (67-115 only): %.4f dB, %.3f deg', 20*log10(abs(ratio)), angle(ratio)*180/pi);
lines{end+1} = '';

for idxBand = 1:numel(baselineResults.bands)
    baseBand = baselineResults.bands{idxBand};
    expBand = experimentResults.bands{idxBand};
    lines{end+1} = sprintf('Band %s', baseBand.band);
    lines{end+1} = sprintf('  DUT return loss baseline / experiment (dB): %.3f / %.3f', baseBand.dut.summary.median_return_loss_db, expBand.dut.summary.median_return_loss_db);
    lines{end+1} = sprintf('  DUT reciprocity baseline / experiment: %.5f / %.5f', baseBand.dut.summary.median_reciprocity_mismatch, expBand.dut.summary.median_reciprocity_mismatch);
    lines{end+1} = sprintf('  DUT isolation surrogate baseline / experiment (dB): %.3f / %.3f', baseBand.dut.summary.isolation_db, expBand.dut.summary.isolation_db);
    for idxPair = 1:numel(baseBand.reference_pair_results)
        bpr = baseBand.reference_pair_results(idxPair);
        epr = expBand.reference_pair_results(idxPair);
        baseS21 = median(mag_db(squeeze(bpr.reference_corrected(2,1,:)) - 0 .* squeeze(bpr.reference_corrected(2,1,:))), 'omitnan'); %#ok<NASGU>
        baseDelta21 = median(mag_db(squeeze(bpr.reference_corrected(2,1,:))) - mag_db(squeeze(bpr.reference_phase1_target.S(2,1,:))), 'omitnan');
        expDelta21 = median(mag_db(squeeze(epr.reference_corrected(2,1,:))) - mag_db(squeeze(epr.reference_phase1_target.S(2,1,:))), 'omitnan');
        baseDelta12 = median(mag_db(squeeze(bpr.reference_corrected(1,2,:))) - mag_db(squeeze(bpr.reference_phase1_target.S(1,2,:))), 'omitnan');
        expDelta12 = median(mag_db(squeeze(epr.reference_corrected(1,2,:))) - mag_db(squeeze(epr.reference_phase1_target.S(1,2,:))), 'omitnan');
        lines{end+1} = sprintf('  %s %s baseline / experiment mismatch S21: %.3f / %.3f dB', bpr.pair, bpr.geometry, baseDelta21, expDelta21);
        lines{end+1} = sprintf('  %s %s baseline / experiment mismatch S12: %.3f / %.3f dB', bpr.pair, bpr.geometry, baseDelta12, expDelta12);
    end
    lines{end+1} = '';
end
write_text_file(fullfile(expRoot, 'PhaseII_MidbandNormalizationExperiment_Summary.txt'), lines);
end

function plot_experiment_comparison(expRoot, baselineResults, experimentResults)
fig = figure('Visible', 'on', 'Color', 'w');
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

bandOrder = {'0-67', '67-115', '110-170'};
baseBands = containers.Map();
expBands = containers.Map();
for idx = 1:numel(baselineResults.bands)
    baseBands(baselineResults.bands{idx}.band) = baselineResults.bands{idx};
    expBands(experimentResults.bands{idx}.band) = experimentResults.bands{idx};
end

ax1 = nexttile; hold(ax1, 'on');
title(ax1, 'Straight-reference corrected mismatch vs Phase I target');
xlabel(ax1, 'Frequency (GHz)'); ylabel(ax1, 'Corrected - target |S_{21}| (dB)'); grid(ax1, 'on'); yline(ax1, 0, '--', 'Color', [0.6 0.6 0.6]);
ax2 = nexttile; hold(ax2, 'on');
title(ax2, 'ARC-reference corrected mismatch vs Phase I target');
xlabel(ax2, 'Frequency (GHz)'); ylabel(ax2, 'Corrected - target |S_{21}| (dB)'); grid(ax2, 'on'); yline(ax2, 0, '--', 'Color', [0.6 0.6 0.6]);

baseColor = [0.85 0.33 0.10];
expColor = [0 0.45 0.74];
keysStraight = {'P1P2_STRAIGHT', 'P3P4_STRAIGHT'};
keysArc = {'P1P3_ARC', 'P1P4_ARC', 'P2P3_ARC', 'P2P4_ARC'};

for idxBand = 1:numel(bandOrder)
    bandName = bandOrder{idxBand};
    baseBand = baseBands(bandName);
    expBand = expBands(bandName);
    for idxPair = 1:numel(baseBand.reference_pair_results)
        key = sprintf('%s_%s', upper(baseBand.reference_pair_results(idxPair).pair), upper(baseBand.reference_pair_results(idxPair).geometry));
        freqGHz = baseBand.reference_pair_results(idxPair).freq / 1e9;
        baseDelta = mag_db(squeeze(baseBand.reference_pair_results(idxPair).reference_corrected(2,1,:))) - mag_db(squeeze(baseBand.reference_pair_results(idxPair).reference_phase1_target.S(2,1,:)));
        expDelta = mag_db(squeeze(expBand.reference_pair_results(idxPair).reference_corrected(2,1,:))) - mag_db(squeeze(expBand.reference_pair_results(idxPair).reference_phase1_target.S(2,1,:)));
        if any(strcmpi(key, keysStraight))
            plot(ax1, freqGHz, baseDelta, '--', 'Color', baseColor, 'LineWidth', 1.0, 'HandleVisibility', conditional_vis(idxBand,1));
            plot(ax1, freqGHz, expDelta, '-', 'Color', expColor, 'LineWidth', 1.4, 'HandleVisibility', conditional_vis(idxBand,1));
        elseif any(strcmpi(key, keysArc))
            plot(ax2, freqGHz, baseDelta, '--', 'Color', baseColor, 'LineWidth', 1.0, 'HandleVisibility', conditional_vis(idxBand,1));
            plot(ax2, freqGHz, expDelta, '-', 'Color', expColor, 'LineWidth', 1.4, 'HandleVisibility', conditional_vis(idxBand,1));
        end
    end
end
h1 = plot(ax1, nan, nan, '--', 'Color', baseColor, 'LineWidth', 1.0, 'DisplayName', 'Baseline');
h2 = plot(ax1, nan, nan, '-', 'Color', expColor, 'LineWidth', 1.4, 'DisplayName', '67-115 normalized experiment');
legend(ax1, [h1 h2], 'Location', 'best');
h3 = plot(ax2, nan, nan, '--', 'Color', baseColor, 'LineWidth', 1.0, 'DisplayName', 'Baseline');
h4 = plot(ax2, nan, nan, '-', 'Color', expColor, 'LineWidth', 1.4, 'DisplayName', '67-115 normalized experiment');
legend(ax2, [h3 h4], 'Location', 'best');

save_debug_figure(fig, fullfile(expRoot, 'figures', 'PhaseII_MidbandNormalizationExperiment_Comparison'));
end

function mode = conditional_vis(idx, expected)
if idx == expected
    mode = 'on';
else
    mode = 'off';
end
end

function values = mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function save_debug_figure(fig, basePath)
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

function write_text_file(filePath, lines)
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
