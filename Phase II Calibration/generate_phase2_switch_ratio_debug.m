function generate_phase2_switch_ratio_debug(config, bandResults)
%GENERATE_PHASE2_SWITCH_RATIO_DEBUG Inspect raw switch-wave magnitudes and ratios.

debug = struct();
noteLines = {};
noteLines{end + 1} = 'Phase II Switch-Ratio Debug Summary';
noteLines{end + 1} = '==================================';
noteLines{end + 1} = sprintf('Active ratio mode: %s', config.switch_ratio_mode);
noteLines{end + 1} = '';

for idxPair = 1:numel(config.switch_debug_pairs)
    pairKey = upper(config.switch_debug_pairs(idxPair).pair);
    geomKey = upper(config.switch_debug_pairs(idxPair).geometry);
    pairLabel = sprintf('%s %s', pairKey, geomKey);

    pairDebug = collect_pair_debug(bandResults, pairKey, geomKey);

    fig = figure('Visible', config.figure_visible, 'Color', 'w');
    tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax = nexttile;
    plot_wave_panel(ax, pairDebug.bands, 'a_forward_raw', 'b_forward_raw', 'a_F', 'b_F', ...
        sprintf('%s forward receiver waves', pairLabel), 'Magnitude (dB)');

    ax = nexttile;
    plot_wave_panel(ax, pairDebug.bands, 'a_reverse_raw', 'b_reverse_raw', 'a_R', 'b_R', ...
        sprintf('%s reverse receiver waves', pairLabel), 'Magnitude (dB)');

    ax = nexttile;
    plot_ratio_panel(ax, pairDebug.bands, 'gamma_forward_raw', '\Gamma_F', ...
        sprintf('%s |\\Gamma_F|', pairLabel), 'Magnitude (dB)', false);

    ax = nexttile;
    plot_ratio_panel(ax, pairDebug.bands, 'gamma_reverse_raw', '\Gamma_R', ...
        sprintf('%s |\\Gamma_R|', pairLabel), 'Magnitude (dB)', false);

    ax = nexttile;
    plot_ratio_panel(ax, pairDebug.bands, 'gamma_forward_raw', '\Gamma_F', ...
        sprintf('%s \\angle\\Gamma_F', pairLabel), 'Phase (deg)', true);
    xlabel(ax, 'Frequency (GHz)');

    ax = nexttile;
    plot_ratio_panel(ax, pairDebug.bands, 'gamma_reverse_raw', '\Gamma_R', ...
        sprintf('%s \\angle\\Gamma_R', pairLabel), 'Phase (deg)', true);
    xlabel(ax, 'Frequency (GHz)');

    style_axes(findall(fig, 'type', 'axes'));
    sgtitle(fig, sprintf('Phase II switch-ratio inspection - %s (%s)', pairLabel, config.switch_ratio_mode));
    save_figure(fig, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_SwitchRatio_%s_%s', pairKey, geomKey)));

    debug.(pairKey).(geomKey) = pairDebug; %#ok<STRNU>

    noteLines{end + 1} = pairLabel;
    noteLines{end + 1} = repmat('-', 1, numel(pairLabel));
    for idxBand = 1:numel(pairDebug.bands)
        bandSlice = pairDebug.bands{idxBand};
        if isempty(bandSlice)
            continue;
        end
        noteLines{end + 1} = sprintf('Band %s', bandSlice.band);
        noteLines{end + 1} = sprintf('  Median |Gamma_F| (dB): %.3f', bandSlice.metrics.gammaF_db);
        noteLines{end + 1} = sprintf('  Median |Gamma_R| (dB): %.3f', bandSlice.metrics.gammaR_db);
        noteLines{end + 1} = sprintf('  Median |a_F| / |b_F| (dB): %.3f / %.3f', bandSlice.metrics.aF_db, bandSlice.metrics.bF_db);
        noteLines{end + 1} = sprintf('  Median |a_R| / |b_R| (dB): %.3f / %.3f', bandSlice.metrics.aR_db, bandSlice.metrics.bR_db);
    end
    noteLines{end + 1} = '';
end

save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_SWITCH_RATIOS.mat'), 'debug');
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_SwitchRatio_Summary.txt'), noteLines);
end

function pairDebug = collect_pair_debug(bandResults, pairKey, geomKey)
pairDebug = struct();
pairDebug.pair = pairKey;
pairDebug.geometry = geomKey;
pairDebug.bands = cell(1, numel(bandResults));

for idxBand = 1:numel(bandResults)
    pairResult = find_pair_result(bandResults{idxBand}, pairKey, geomKey);
    if isempty(pairResult)
        continue;
    end
    pairDebug.bands{idxBand} = build_band_debug(bandResults{idxBand}.band, pairResult.switch_terms);
end
end

function plot_wave_panel(ax, bands, fieldA, fieldB, labelA, labelB, titleText, yLabelText)
hold(ax, 'on');
for idxBand = 1:numel(bands)
    bandSlice = bands{idxBand};
    if isempty(bandSlice)
        continue;
    end
    plot(ax, bandSlice.freqA / 1e9, mag_db(bandSlice.(fieldA)), 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%s %s', labelA, bandSlice.band));
    plot(ax, bandSlice.freqB / 1e9, mag_db(bandSlice.(fieldB)), '--', 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%s %s', labelB, bandSlice.band));
end
title(ax, titleText);
ylabel(ax, yLabelText);
legend(ax, 'Location', 'best');
end

function plot_ratio_panel(ax, bands, fieldName, labelPrefix, titleText, yLabelText, isPhase)
hold(ax, 'on');
for idxBand = 1:numel(bands)
    bandSlice = bands{idxBand};
    if isempty(bandSlice)
        continue;
    end
    values = bandSlice.(fieldName);
    if isPhase
        y = unwrap(angle(values)) * 180 / pi;
    else
        y = mag_db(values);
    end
    plot(ax, bandSlice.freqGamma / 1e9, y, 'LineWidth', 1.3, ...
        'DisplayName', sprintf('%s %s', labelPrefix, bandSlice.band));
end
title(ax, titleText);
ylabel(ax, yLabelText);
legend(ax, 'Location', 'best');
end

function style_axes(ax)
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
for idxAx = 1:numel(ax)
    xline(ax(idxAx), 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax(idxAx), 115, ':', 'Color', [0.7 0.7 0.7]);
    grid(ax(idxAx), 'on');
end
end

function bandDebug = build_band_debug(bandName, switchTerms)
bandDebug = struct();
bandDebug.band = bandName;
bandDebug.freqA = switchTerms.freq_forward_raw;
bandDebug.freqB = switchTerms.freq_reverse_raw;
bandDebug.freqGamma = switchTerms.freq_forward_raw;
bandDebug.a_forward_raw = switchTerms.a_forward_raw;
bandDebug.b_forward_raw = switchTerms.b_forward_raw;
bandDebug.a_reverse_raw = switchTerms.a_reverse_raw;
bandDebug.b_reverse_raw = switchTerms.b_reverse_raw;
bandDebug.gamma_forward_raw = switchTerms.gamma_forward_raw;
bandDebug.gamma_reverse_raw = switchTerms.gamma_reverse_raw;
bandDebug.metrics = struct( ...
    'gammaF_db', median(mag_db(switchTerms.gamma_forward_raw), 'omitnan'), ...
    'gammaR_db', median(mag_db(switchTerms.gamma_reverse_raw), 'omitnan'), ...
    'aF_db', median(mag_db(switchTerms.a_forward_raw), 'omitnan'), ...
    'bF_db', median(mag_db(switchTerms.b_forward_raw), 'omitnan'), ...
    'aR_db', median(mag_db(switchTerms.a_reverse_raw), 'omitnan'), ...
    'bR_db', median(mag_db(switchTerms.b_reverse_raw), 'omitnan'));
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
