function generate_phase2_pair_waterfall_debug(config, bandResults)
%GENERATE_PHASE2_PAIR_WATERFALL_DEBUG Plot raw/switch/corrected thru stages.

pairs = config.reference_pairs;
debug = struct();
noteLines = {};
noteLines{end + 1} = 'Phase II Pairwise Thru Waterfall Summary';
noteLines{end + 1} = '======================================';
noteLines{end + 1} = '';

for idxPair = 1:numel(pairs)
    pairKey = upper(pairs(idxPair).pair);
    geomKey = upper(pairs(idxPair).geometry);
    pairLabel = sprintf('%s %s', pairKey, geomKey);

    fig = figure('Visible', config.figure_visible, 'Color', 'w');
    tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    stageColors = [0.45 0.45 0.45; 0.85 0.33 0.10; 0.00 0.45 0.74];
    stageStyles = {'--', '-.', '-'};
    stageWidths = [1.0, 1.1, 1.8];

    pairDebug = struct();
    pairDebug.pair = pairKey;
    pairDebug.geometry = geomKey;
    pairDebug.bands = cell(1, numel(bandResults));

    plot_panel(1, 2, 1, '|S_{21}| (dB)');
    plot_panel(2, 1, 2, '|S_{12}| (dB)');
    plot_panel(1, 1, 3, '|S_{11}| (dB)');
    plot_panel(2, 2, 4, '|S_{22}| (dB)');

    sgtitle(fig, sprintf('Phase II pairwise thru waterfall - %s', pairLabel));
    save_figure(fig, fullfile(config.interim_figure_dir, ...
        sprintf('PhaseII_Debug_PairWaterfall_%s_%s', pairKey, geomKey)));

    debug.(pairKey).(geomKey) = pairDebug; %#ok<STRNU>

    noteLines{end + 1} = pairLabel;
    noteLines{end + 1} = repmat('-', 1, numel(pairLabel));
    for idxBand = 1:numel(pairDebug.bands)
        bandSlice = pairDebug.bands{idxBand};
        if isempty(bandSlice)
            continue;
        end
        noteLines{end + 1} = sprintf('Band %s', bandSlice.band);
        noteLines{end + 1} = sprintf('  Median |S21| raw / switch / corrected (dB): %.3f / %.3f / %.3f', ...
            bandSlice.metrics.s21_db(1), bandSlice.metrics.s21_db(2), bandSlice.metrics.s21_db(3));
        noteLines{end + 1} = sprintf('  Median |S12| raw / switch / corrected (dB): %.3f / %.3f / %.3f', ...
            bandSlice.metrics.s12_db(1), bandSlice.metrics.s12_db(2), bandSlice.metrics.s12_db(3));
        noteLines{end + 1} = sprintf('  Median |S11| raw / switch / corrected (dB): %.3f / %.3f / %.3f', ...
            bandSlice.metrics.s11_db(1), bandSlice.metrics.s11_db(2), bandSlice.metrics.s11_db(3));
        noteLines{end + 1} = sprintf('  Median |S22| raw / switch / corrected (dB): %.3f / %.3f / %.3f', ...
            bandSlice.metrics.s22_db(1), bandSlice.metrics.s22_db(2), bandSlice.metrics.s22_db(3));
    end
    noteLines{end + 1} = '';
end

save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_PAIR_WATERFALL.mat'), 'debug');
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_PairWaterfall_Summary.txt'), noteLines);

    function plot_panel(row, col, tileIdx, yLabelText)
        nexttile(tileIdx);
        hold on;
        proxyHandles = gobjects(1, 3);
        for idxStage = 1:3
            proxyHandles(idxStage) = plot(nan, nan, stageStyles{idxStage}, ...
                'Color', stageColors(idxStage, :), 'LineWidth', stageWidths(idxStage));
        end

        for idxBand = 1:numel(bandResults)
            pairResult = find_pair_result(bandResults{idxBand}, pairKey, geomKey);
            if isempty(pairResult)
                continue;
            end
            freqGHz = pairResult.freq(:) / 1e9;
            rawTrace = squeeze(pairResult.reference_raw(row, col, :));
            switchTrace = squeeze(pairResult.reference_switch_corrected(row, col, :));
            corrTrace = squeeze(pairResult.reference_corrected(row, col, :));

            plot(freqGHz, mag_db(rawTrace), stageStyles{1}, 'Color', stageColors(1, :), 'LineWidth', stageWidths(1));
            plot(freqGHz, mag_db(switchTrace), stageStyles{2}, 'Color', stageColors(2, :), 'LineWidth', stageWidths(2));
            plot(freqGHz, mag_db(corrTrace), stageStyles{3}, 'Color', stageColors(3, :), 'LineWidth', stageWidths(3));

            pairDebug.bands{idxBand} = struct( ...
                'band', bandResults{idxBand}.band, ...
                'freq', pairResult.freq(:), ...
                'raw', pairResult.reference_raw, ...
                'switch_corrected', pairResult.reference_switch_corrected, ...
                'corrected', pairResult.reference_corrected, ...
                'metrics', summarize_stage_metrics(pairResult));
        end

        xline(67, ':', 'Color', [0.7 0.7 0.7]);
        xline(115, ':', 'Color', [0.7 0.7 0.7]);
        grid on;
        xlabel('Frequency (GHz)');
        ylabel(yLabelText);
        title(sprintf('%s %s', pairLabel, yLabelText));
        legend(proxyHandles, {'Raw', 'Switch-corrected', 'SOLR-corrected'}, 'Location', 'best');
    end
end

function pairResult = find_pair_result(bandResult, pairKey, geomKey)
pairResult = [];
for idx = 1:numel(bandResult.reference_pair_results)
    cand = bandResult.reference_pair_results(idx);
    if strcmpi(cand.pair, pairKey) && strcmpi(cand.geometry, geomKey)
        pairResult = cand;
        return;
    end
end
end

function metrics = summarize_stage_metrics(pairResult)
metrics = struct();
metrics.s21_db = local_stage_db(pairResult, 2, 1);
metrics.s12_db = local_stage_db(pairResult, 1, 2);
metrics.s11_db = local_stage_db(pairResult, 1, 1);
metrics.s22_db = local_stage_db(pairResult, 2, 2);
end

function values = local_stage_db(pairResult, row, col)
values = [ ...
    median(mag_db(squeeze(pairResult.reference_raw(row, col, :))), 'omitnan'), ...
    median(mag_db(squeeze(pairResult.reference_switch_corrected(row, col, :))), 'omitnan'), ...
    median(mag_db(squeeze(pairResult.reference_corrected(row, col, :))), 'omitnan')];
end

function values = mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function save_figure(fig, basePath)
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
fid = fopen(filePath, 'w');
cleanup = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end
