function comparison = run_phase2_sol_reference_experiment()
%RUN_PHASE2_SOL_REFERENCE_EXPERIMENT Compare stitched vs bandwise Phase I SOL references.

baselineConfig = phase2_config();
baselineFile = fullfile(baselineConfig.mat_dir, 'PHASE2_RESULTS_ALL.mat');
if ~exist(baselineFile, 'file')
    error('Baseline Phase II results not found: %s', baselineFile);
end

loaded = load(baselineFile, 'results');
baselineResults = loaded.results;

expConfig = phase2_config();
expConfig.reference_standard_mode = 'bandwise_final';
expConfig.output_dir = fullfile(expConfig.interim_dir, 'solReferenceExperiment', 'outputs');
expConfig.touchstone_dir = fullfile(expConfig.output_dir, 'touchstone');
expConfig.mat_dir = fullfile(expConfig.output_dir, 'mat');
expConfig.figure_dir = fullfile(expConfig.output_dir, 'figures');
expConfig.note_dir = fullfile(expConfig.output_dir, 'notes');
expConfig.interim_dir = fullfile(expConfig.interim_dir, 'solReferenceExperiment', 'interim');
expConfig.interim_mat_dir = fullfile(expConfig.interim_dir, 'mat');
expConfig.interim_figure_dir = fullfile(expConfig.interim_dir, 'figures');
expConfig.interim_note_dir = fullfile(expConfig.interim_dir, 'notes');

experimentResults = run_phase2_solr(expConfig);
comparison = compare_reference_modes(baselineResults, experimentResults);

comparisonDir = fullfile(baselineConfig.interim_dir, 'solReferenceExperiment');
if ~exist(comparisonDir, 'dir')
    mkdir(comparisonDir);
end
save(fullfile(comparisonDir, 'PHASE2_SOL_REFERENCE_EXPERIMENT.mat'), 'baselineResults', 'experimentResults', 'comparison');
write_experiment_note(fullfile(comparisonDir, 'PhaseII_SOL_ReferenceExperiment_Summary.txt'), comparison);
plot_experiment_comparison(fullfile(comparisonDir, 'figures'), baselineResults, experimentResults, comparison);
end

function comparison = compare_reference_modes(baselineResults, experimentResults)
comparison = struct();
comparison.band_labels = cell(1, numel(baselineResults.bands));
comparison.pairs = struct();
for idxBand = 1:numel(baselineResults.bands)
    baseBand = baselineResults.bands{idxBand};
    expBand = experimentResults.bands{idxBand};
    comparison.band_labels{idxBand} = baseBand.band;

    comparison.dut(idxBand) = struct( ... %#ok<AGROW>
        'band', baseBand.band, ...
        'baseline_return_loss_db', baseBand.dut.summary.median_return_loss_db, ...
        'experiment_return_loss_db', expBand.dut.summary.median_return_loss_db, ...
        'baseline_recip', baseBand.dut.summary.median_reciprocity_mismatch, ...
        'experiment_recip', expBand.dut.summary.median_reciprocity_mismatch, ...
        'baseline_isolation_db', baseBand.dut.summary.isolation_db, ...
        'experiment_isolation_db', expBand.dut.summary.isolation_db);

    for idxPair = 1:numel(baseBand.reference_pair_results)
        basePair = baseBand.reference_pair_results(idxPair);
        expPair = expBand.reference_pair_results(idxPair);
        key = sprintf('%s_%s', basePair.pair, upper(basePair.geometry));
        if ~isfield(comparison.pairs, key)
            comparison.pairs.(key) = struct('bands', struct());
        end
        targetS21 = squeeze(basePair.reference_phase1_target.S(2, 1, :));
        targetS12 = squeeze(basePair.reference_phase1_target.S(1, 2, :));
        baseS21 = squeeze(basePair.reference_corrected(2, 1, :));
        expS21 = squeeze(expPair.reference_corrected(2, 1, :));
        baseS12 = squeeze(basePair.reference_corrected(1, 2, :));
        expS12 = squeeze(expPair.reference_corrected(1, 2, :));
        bandKey = matlab.lang.makeValidName(strrep(baseBand.band, '-', '_'));
        comparison.pairs.(key).bands.(bandKey) = struct( ...
            'band', baseBand.band, ...
            'freq', basePair.freq(:), ...
            'target_s21_db', mag_db(targetS21), ...
            'baseline_s21_db', mag_db(baseS21), ...
            'experiment_s21_db', mag_db(expS21), ...
            'target_s12_db', mag_db(targetS12), ...
            'baseline_s12_db', mag_db(baseS12), ...
            'experiment_s12_db', mag_db(expS12), ...
            'median_baseline_s21_target_delta_db', median(mag_db(baseS21) - mag_db(targetS21), 'omitnan'), ...
            'median_experiment_s21_target_delta_db', median(mag_db(expS21) - mag_db(targetS21), 'omitnan'), ...
            'median_baseline_s12_target_delta_db', median(mag_db(baseS12) - mag_db(targetS12), 'omitnan'), ...
            'median_experiment_s12_target_delta_db', median(mag_db(expS12) - mag_db(targetS12), 'omitnan'));
    end
end
end

function write_experiment_note(filePath, comparison)
lines = {
    'Phase II SOL reference experiment summary'
    '========================================'
    'Baseline: stitched Phase I SOL references'
    'Experiment: bandwise Phase I SOL winners'
    ''
    'DUT summaries by band'
    };
for idxBand = 1:numel(comparison.dut)
    item = comparison.dut(idxBand);
    lines{end + 1} = sprintf('%s', item.band); %#ok<AGROW>
    lines{end + 1} = sprintf('  Median DUT return loss baseline / experiment (dB): %.3f / %.3f', ...
        item.baseline_return_loss_db, item.experiment_return_loss_db);
    lines{end + 1} = sprintf('  Median DUT reciprocity baseline / experiment: %.5f / %.5f', ...
        item.baseline_recip, item.experiment_recip);
    lines{end + 1} = sprintf('  Isolation surrogate baseline / experiment (dB): %.3f / %.3f', ...
        item.baseline_isolation_db, item.experiment_isolation_db);
end
lines{end + 1} = '';
lines{end + 1} = 'Reference-pair median target mismatch (S21 / S12, dB)';

pairs = fieldnames(comparison.pairs);
for idxPair = 1:numel(pairs)
    pairData = comparison.pairs.(pairs{idxPair});
    lines{end + 1} = pairs{idxPair}; %#ok<AGROW>
    bandKeys = fieldnames(pairData.bands);
    for idxBand = 1:numel(bandKeys)
        bandData = pairData.bands.(bandKeys{idxBand});
        lines{end + 1} = sprintf('  %s baseline / experiment S21 delta: %.3f / %.3f', ...
            bandData.band, bandData.median_baseline_s21_target_delta_db, bandData.median_experiment_s21_target_delta_db);
        lines{end + 1} = sprintf('  %s baseline / experiment S12 delta: %.3f / %.3f', ...
            bandData.band, bandData.median_baseline_s12_target_delta_db, bandData.median_experiment_s12_target_delta_db);
    end
end

write_text_file(filePath, lines);
end

function plot_experiment_comparison(figDir, baselineResults, experimentResults, comparison)
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

fig = figure('Visible', 'on', 'Color', 'w', 'Position', [100 100 1400 900]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

freqBase = baselineResults.stitched.dut.freq / 1e9;
freqExp = experimentResults.stitched.dut.freq / 1e9;

ax1 = nexttile;
hold(ax1, 'on');
plot(ax1, freqBase, mag_db(squeeze(baselineResults.stitched.dut.S(2,1,:))), '-', 'LineWidth', 1.5, 'DisplayName', 'Baseline S21');
plot(ax1, freqExp, mag_db(squeeze(experimentResults.stitched.dut.S(2,1,:))), '--', 'LineWidth', 1.5, 'DisplayName', 'Experiment S21');
plot(ax1, freqBase, mag_db(squeeze(baselineResults.stitched.dut.S(3,1,:))), '-', 'LineWidth', 1.2, 'DisplayName', 'Baseline S31');
plot(ax1, freqExp, mag_db(squeeze(experimentResults.stitched.dut.S(3,1,:))), '--', 'LineWidth', 1.2, 'DisplayName', 'Experiment S31');
plot(ax1, freqBase, mag_db(squeeze(baselineResults.stitched.dut.S(4,1,:))), '-', 'LineWidth', 1.2, 'DisplayName', 'Baseline S41');
plot(ax1, freqExp, mag_db(squeeze(experimentResults.stitched.dut.S(4,1,:))), '--', 'LineWidth', 1.2, 'DisplayName', 'Experiment S41');
grid(ax1, 'on');
xlabel(ax1, 'Frequency (GHz)');
ylabel(ax1, '|S_{x1}| (dB)');
title(ax1, 'Stitched DUT comparison');
legend(ax1, 'Location', 'best');

pairKeys = {'P1P2_STRAIGHT', 'P3P4_STRAIGHT', 'P1P3_ARC'};
for idx = 1:min(3, numel(pairKeys))
    ax = nexttile;
    hold(ax, 'on');
    pairData = comparison.pairs.(pairKeys{idx});
    bandNames = fieldnames(pairData.bands);
    for idxBand = 1:numel(bandNames)
        bandData = pairData.bands.(bandNames{idxBand});
        freqGHz = bandData.freq / 1e9;
        plot(ax, freqGHz, bandData.target_s21_db, 'k:', 'LineWidth', 1.5, 'HandleVisibility', conditional_visibility(idxBand, 1));
        plot(ax, freqGHz, bandData.baseline_s21_db, '-', 'LineWidth', 1.4, 'HandleVisibility', conditional_visibility(idxBand, 1));
        plot(ax, freqGHz, bandData.experiment_s21_db, '--', 'LineWidth', 1.4, 'HandleVisibility', conditional_visibility(idxBand, 1));
    end
    grid(ax, 'on');
    xlabel(ax, 'Frequency (GHz)');
    ylabel(ax, '|S_{21}| (dB)');
    title(ax, strrep(pairKeys{idx}, '_', ' '));
    h1 = plot(ax, nan, nan, 'k:', 'LineWidth', 1.5, 'DisplayName', 'Target');
    h2 = plot(ax, nan, nan, '-', 'LineWidth', 1.4, 'DisplayName', 'Baseline');
    h3 = plot(ax, nan, nan, '--', 'LineWidth', 1.4, 'DisplayName', 'Bandwise SOL experiment');
    legend(ax, [h1 h2 h3], 'Location', 'best');
end

exportgraphics(fig, fullfile(figDir, 'PhaseII_SOL_ReferenceExperiment_Comparison.jpg'), 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, fullfile(figDir, 'PhaseII_SOL_ReferenceExperiment_Comparison.pdf'), 'BackgroundColor', 'white', 'ContentType', 'vector');
savefig(fig, fullfile(figDir, 'PhaseII_SOL_ReferenceExperiment_Comparison.fig'));
close(fig);
end

function textOut = conditional_visibility(idxBand, desiredBand)
if idxBand == desiredBand
    textOut = 'on';
else
    textOut = 'off';
end
end

function values = mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function write_text_file(filePath, lines)
fid = fopen(filePath, 'w');
cleanup = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end
