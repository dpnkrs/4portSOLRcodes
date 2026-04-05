function generate_phase2_geometry_consistency_debug(config, bandResults)
%GENERATE_PHASE2_GEOMETRY_CONSISTENCY_DEBUG Compare ARC reference vs DIAG holdout within each pair.

pairs = {'P1P3', 'P1P4', 'P2P3', 'P2P4'};
debug = struct();
noteLines = {};
noteLines{end + 1} = 'Phase II ARC vs DIAG Geometry Consistency Audit';
noteLines{end + 1} = '============================================';
noteLines{end + 1} = '';

for idxPair = 1:numel(pairs)
    pairKey = pairs{idxPair};
    fig = figure('Visible', config.figure_visible, 'Color', 'w');
    tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    pairDebug = struct('pair', pairKey, 'bands', {cell(1, numel(bandResults))});
    bandColors = lines(3);

    for idxBand = 1:numel(bandResults)
        bandResult = bandResults{idxBand};
        refResult = find_pair_result(bandResult.reference_pair_results, pairKey, 'ARC');
        holdResult = find_pair_result(bandResult.holdout_results, pairKey, 'DIAG');
        if isempty(refResult) || isempty(holdResult)
            continue;
        end

        freqGHz = refResult.freq(:) / 1e9;
        axMag = nexttile((idxBand - 1) * 2 + 1);
        hold(axMag, 'on');
        title(axMag, sprintf('%s | %s transmission', pairKey, bandResult.band));
        xlabel(axMag, 'Frequency (GHz)');
        ylabel(axMag, '|S_{21}| (dB)');
        grid(axMag, 'on');

        axDelta = nexttile((idxBand - 1) * 2 + 2);
        hold(axDelta, 'on');
        title(axDelta, sprintf('%s | %s DIAG-ARC delta', pairKey, bandResult.band));
        xlabel(axDelta, 'Frequency (GHz)');
        ylabel(axDelta, 'Delta (dB)');
        yline(axDelta, 0, '--', 'Color', [0.6 0.6 0.6]);
        grid(axDelta, 'on');

        refRaw = squeeze(refResult.reference_raw(2, 1, :));
        refSwitch = squeeze(refResult.reference_switch_corrected(2, 1, :));
        refCorr = squeeze(refResult.reference_corrected(2, 1, :));
        holdRaw = squeeze(holdResult.S_raw(2, 1, :));
        holdSwitch = squeeze(holdResult.S_switch(2, 1, :));
        holdCorr = squeeze(holdResult.S_corrected(2, 1, :));

        refRawDb = mag_db(refRaw);
        refSwitchDb = mag_db(refSwitch);
        refCorrDb = mag_db(refCorr);
        holdRawDb = mag_db(holdRaw);
        holdSwitchDb = mag_db(holdSwitch);
        holdCorrDb = mag_db(holdCorr);

        plot(axMag, freqGHz, refRawDb, '--', 'Color', bandColors(1, :), 'LineWidth', 1.0);
        plot(axMag, freqGHz, holdRawDb, ':', 'Color', bandColors(1, :), 'LineWidth', 1.2);
        plot(axMag, freqGHz, refSwitchDb, '--', 'Color', bandColors(2, :), 'LineWidth', 1.0);
        plot(axMag, freqGHz, holdSwitchDb, ':', 'Color', bandColors(2, :), 'LineWidth', 1.2);
        plot(axMag, freqGHz, refCorrDb, '-', 'Color', bandColors(3, :), 'LineWidth', 1.4);
        plot(axMag, freqGHz, holdCorrDb, '-', 'Color', 0.55 * bandColors(3, :) + 0.45, 'LineWidth', 1.4);

        plot(axDelta, freqGHz, holdRawDb - refRawDb, '--', 'Color', bandColors(1, :), 'LineWidth', 1.1);
        plot(axDelta, freqGHz, holdSwitchDb - refSwitchDb, '--', 'Color', bandColors(2, :), 'LineWidth', 1.1);
        plot(axDelta, freqGHz, holdCorrDb - refCorrDb, '-', 'Color', bandColors(3, :), 'LineWidth', 1.4);

        legend(axMag, {'ARC raw', 'DIAG raw', 'ARC switch', 'DIAG switch', 'ARC corrected', 'DIAG corrected'}, 'Location', 'best');
        legend(axDelta, {'Raw delta', 'Switch delta', 'Corrected delta'}, 'Location', 'best');

        bandDebug = struct();
        bandDebug.band = bandResult.band;
        bandDebug.freq = refResult.freq(:);
        bandDebug.metrics = struct( ...
            'median_arc_raw_s21_db', median(refRawDb, 'omitnan'), ...
            'median_diag_raw_s21_db', median(holdRawDb, 'omitnan'), ...
            'median_arc_switch_s21_db', median(refSwitchDb, 'omitnan'), ...
            'median_diag_switch_s21_db', median(holdSwitchDb, 'omitnan'), ...
            'median_arc_corrected_s21_db', median(refCorrDb, 'omitnan'), ...
            'median_diag_corrected_s21_db', median(holdCorrDb, 'omitnan'), ...
            'median_diag_minus_arc_raw_db', median(holdRawDb - refRawDb, 'omitnan'), ...
            'median_diag_minus_arc_switch_db', median(holdSwitchDb - refSwitchDb, 'omitnan'), ...
            'median_diag_minus_arc_corrected_db', median(holdCorrDb - refCorrDb, 'omitnan'));
        pairDebug.bands{idxBand} = bandDebug;

        noteLines{end + 1} = sprintf('%s | %s', pairKey, bandResult.band);
        noteLines{end + 1} = sprintf('  ARC raw/switch/corrected |S21| (dB): %.3f / %.3f / %.3f', ...
            bandDebug.metrics.median_arc_raw_s21_db, ...
            bandDebug.metrics.median_arc_switch_s21_db, ...
            bandDebug.metrics.median_arc_corrected_s21_db);
        noteLines{end + 1} = sprintf('  DIAG raw/switch/corrected |S21| (dB): %.3f / %.3f / %.3f', ...
            bandDebug.metrics.median_diag_raw_s21_db, ...
            bandDebug.metrics.median_diag_switch_s21_db, ...
            bandDebug.metrics.median_diag_corrected_s21_db);
        noteLines{end + 1} = sprintf('  DIAG-ARC raw/switch/corrected delta (dB): %.3f / %.3f / %.3f', ...
            bandDebug.metrics.median_diag_minus_arc_raw_db, ...
            bandDebug.metrics.median_diag_minus_arc_switch_db, ...
            bandDebug.metrics.median_diag_minus_arc_corrected_db);
        noteLines{end + 1} = '';
    end

    sgtitle(fig, sprintf('Phase II ARC vs DIAG geometry audit - %s', pairKey));
    save_debug_figure(fig, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_GeometryAudit_%s', pairKey)));
    debug.(pairKey) = pairDebug; %#ok<STRNU>
end

save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_GEOMETRY_CONSISTENCY.mat'), 'debug');
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_GeometryConsistency_Summary.txt'), noteLines);
end

function result = find_pair_result(results, pairKey, geometry)
result = [];
for idx = 1:numel(results)
    cand = results(idx);
    if strcmpi(cand.pair, pairKey) && strcmpi(cand.geometry, geometry)
        result = cand;
        return;
    end
end
end

function values = mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function save_debug_figure(fig, basePath)
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
