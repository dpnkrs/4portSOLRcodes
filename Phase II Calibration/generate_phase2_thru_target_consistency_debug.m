function generate_phase2_thru_target_consistency_debug(config, bandResults)
%GENERATE_PHASE2_THRU_TARGET_CONSISTENCY_DEBUG Compare raw/switched references to Phase I thru targets.

pairs = config.reference_pairs;
debug = struct();
noteLines = {};
noteLines{end + 1} = 'Phase II Reference-Thru Target Consistency Summary';
noteLines{end + 1} = '===============================================';
noteLines{end + 1} = '';

for idxPair = 1:numel(pairs)
    pairKey = upper(pairs(idxPair).pair);
    geomKey = upper(pairs(idxPair).geometry);
    pairLabel = sprintf('%s %s', pairKey, geomKey);

    fig = figure('Visible', config.figure_visible, 'Color', 'w');
    tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold(ax1, 'on');
    title(ax1, sprintf('%s S_{21} target consistency', pairLabel));
    xlabel(ax1, 'Frequency (GHz)');
    ylabel(ax1, '|S_{21}| (dB)');
    grid(ax1, 'on');
    xline(ax1, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax1, 115, ':', 'Color', [0.7 0.7 0.7]);

    ax2 = nexttile;
    hold(ax2, 'on');
    title(ax2, sprintf('%s S_{12} target consistency', pairLabel));
    xlabel(ax2, 'Frequency (GHz)');
    ylabel(ax2, '|S_{12}| (dB)');
    grid(ax2, 'on');
    xline(ax2, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax2, 115, ':', 'Color', [0.7 0.7 0.7]);

    ax3 = nexttile;
    hold(ax3, 'on');
    title(ax3, sprintf('%s S_{21} phase', pairLabel));
    xlabel(ax3, 'Frequency (GHz)');
    ylabel(ax3, 'Phase (deg)');
    grid(ax3, 'on');
    xline(ax3, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax3, 115, ':', 'Color', [0.7 0.7 0.7]);

    ax4 = nexttile;
    hold(ax4, 'on');
    title(ax4, sprintf('%s target mismatch', pairLabel));
    xlabel(ax4, 'Frequency (GHz)');
    ylabel(ax4, 'Delta (dB)');
    grid(ax4, 'on');
    yline(ax4, 0, '--', 'Color', [0.6 0.6 0.6]);
    xline(ax4, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax4, 115, ':', 'Color', [0.7 0.7 0.7]);

    proxy = gobjects(1, 4);
    proxy(1) = plot(ax1, nan, nan, ':', 'Color', [0 0 0], 'LineWidth', 1.8);
    proxy(2) = plot(ax1, nan, nan, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
    proxy(3) = plot(ax1, nan, nan, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3);
    proxy(4) = plot(ax1, nan, nan, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
    legend(ax1, proxy, {'Phase I target', 'Raw downstream', 'Switch-corrected', 'Raw-Switch bridge'}, 'Location', 'best');
    legend(ax2, proxy, {'Phase I target', 'Raw downstream', 'Switch-corrected', 'Raw-Switch bridge'}, 'Location', 'best');

    proxyPhase = gobjects(1, 3);
    proxyPhase(1) = plot(ax3, nan, nan, ':', 'Color', [0 0 0], 'LineWidth', 1.8);
    proxyPhase(2) = plot(ax3, nan, nan, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
    proxyPhase(3) = plot(ax3, nan, nan, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3);
    legend(ax3, proxyPhase, {'Phase I target', 'Raw downstream', 'Switch-corrected'}, 'Location', 'best');

    proxyDelta = gobjects(1, 4);
    proxyDelta(1) = plot(ax4, nan, nan, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
    proxyDelta(2) = plot(ax4, nan, nan, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3);
    proxyDelta(3) = plot(ax4, nan, nan, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
    proxyDelta(4) = plot(ax4, nan, nan, '-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);
    legend(ax4, proxyDelta, {'Raw-target S_{21}', 'Switch-target S_{21}', 'Raw-target S_{12}', 'Switch-target S_{12}'}, 'Location', 'best');

    pairDebug = struct('pair', pairKey, 'geometry', geomKey, 'bands', {cell(1, numel(bandResults))});

    for idxBand = 1:numel(bandResults)
        pairResult = find_pair_result(bandResults{idxBand}, pairKey, geomKey);
        if isempty(pairResult)
            continue;
        end

        freq = pairResult.freq(:);
        freqGHz = freq / 1e9;
        target = pairResult.reference_phase1_target;
        rawS21 = squeeze(pairResult.reference_raw(2, 1, :));
        rawS12 = squeeze(pairResult.reference_raw(1, 2, :));
        swS21 = squeeze(pairResult.reference_switch_corrected(2, 1, :));
        swS12 = squeeze(pairResult.reference_switch_corrected(1, 2, :));
        targetS21 = squeeze(target.S(2, 1, :));
        targetS12 = squeeze(target.S(1, 2, :));

        rawS21db = local_mag_db(rawS21);
        rawS12db = local_mag_db(rawS12);
        swS21db = local_mag_db(swS21);
        swS12db = local_mag_db(swS12);
        targetS21db = local_mag_db(targetS21);
        targetS12db = local_mag_db(targetS12);

        plot(ax1, freqGHz, targetS21db, ':', 'Color', [0 0 0], 'LineWidth', 1.8);
        plot(ax1, freqGHz, rawS21db, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
        plot(ax1, freqGHz, swS21db, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3);
        plot(ax2, freqGHz, targetS12db, ':', 'Color', [0 0 0], 'LineWidth', 1.8);
        plot(ax2, freqGHz, rawS12db, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
        plot(ax2, freqGHz, swS12db, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3);

        targetPhase21 = unwrap(angle(targetS21)) * 180 / pi;
        rawPhase21 = align_phase_branch(unwrap(angle(rawS21)) * 180 / pi, targetPhase21);
        swPhase21 = align_phase_branch(unwrap(angle(swS21)) * 180 / pi, targetPhase21);
        plot(ax3, freqGHz, targetPhase21, ':', 'Color', [0 0 0], 'LineWidth', 1.8);
        plot(ax3, freqGHz, rawPhase21, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
        plot(ax3, freqGHz, swPhase21, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3);

        rawDelta21 = rawS21db - targetS21db;
        swDelta21 = swS21db - targetS21db;
        rawDelta12 = rawS12db - targetS12db;
        swDelta12 = swS12db - targetS12db;
        plot(ax4, freqGHz, rawDelta21, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.2);
        plot(ax4, freqGHz, swDelta21, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.3);
        plot(ax4, freqGHz, rawDelta12, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
        plot(ax4, freqGHz, swDelta12, '-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);

        bandDebug = struct();
        bandDebug.band = bandResults{idxBand}.band;
        bandDebug.freq = freq;
        bandDebug.target = target;
        bandDebug.reference_raw = pairResult.reference_raw;
        bandDebug.reference_switch_corrected = pairResult.reference_switch_corrected;
        bandDebug.metrics = struct( ...
            'median_target_s21_db', median(targetS21db, 'omitnan'), ...
            'median_raw_s21_db', median(rawS21db, 'omitnan'), ...
            'median_switch_s21_db', median(swS21db, 'omitnan'), ...
            'median_raw_target_s21_delta_db', median(rawDelta21, 'omitnan'), ...
            'median_switch_target_s21_delta_db', median(swDelta21, 'omitnan'), ...
            'median_target_s12_db', median(targetS12db, 'omitnan'), ...
            'median_raw_s12_db', median(rawS12db, 'omitnan'), ...
            'median_switch_s12_db', median(swS12db, 'omitnan'), ...
            'median_raw_target_s12_delta_db', median(rawDelta12, 'omitnan'), ...
            'median_switch_target_s12_delta_db', median(swDelta12, 'omitnan'));
        pairDebug.bands{idxBand} = bandDebug;

        noteLines{end + 1} = sprintf('%s | %s', pairLabel, bandResults{idxBand}.band);
        noteLines{end + 1} = sprintf('  Median target/raw/switch |S21| (dB): %.3f / %.3f / %.3f', ...
            bandDebug.metrics.median_target_s21_db, bandDebug.metrics.median_raw_s21_db, bandDebug.metrics.median_switch_s21_db);
        noteLines{end + 1} = sprintf('  Median raw-target / switch-target |S21| delta (dB): %.3f / %.3f', ...
            bandDebug.metrics.median_raw_target_s21_delta_db, bandDebug.metrics.median_switch_target_s21_delta_db);
        noteLines{end + 1} = sprintf('  Median target/raw/switch |S12| (dB): %.3f / %.3f / %.3f', ...
            bandDebug.metrics.median_target_s12_db, bandDebug.metrics.median_raw_s12_db, bandDebug.metrics.median_switch_s12_db);
        noteLines{end + 1} = sprintf('  Median raw-target / switch-target |S12| delta (dB): %.3f / %.3f', ...
            bandDebug.metrics.median_raw_target_s12_delta_db, bandDebug.metrics.median_switch_target_s12_delta_db);
        noteLines{end + 1} = '';
    end

    sgtitle(fig, sprintf('Phase II thru target consistency - %s', pairLabel));
    save_debug_figure(fig, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_TargetConsistency_%s_%s', pairKey, geomKey)));
    debug.(pairKey).(geomKey) = pairDebug; %#ok<STRNU>
end

save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_TARGET_CONSISTENCY.mat'), 'debug');
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_TargetConsistency_Summary.txt'), noteLines);
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

function aligned = align_phase_branch(valuesDeg, targetDeg)
if isempty(valuesDeg) || isempty(targetDeg)
    aligned = valuesDeg;
    return;
end
delta = median(targetDeg - valuesDeg, 'omitnan');
if ~isfinite(delta)
    delta = 0;
end
aligned = valuesDeg + 360 * round(delta / 360);
end

function values = local_mag_db(trace)
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
