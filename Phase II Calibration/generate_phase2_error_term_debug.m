function generate_phase2_error_term_debug(config, bandResults)
%GENERATE_PHASE2_ERROR_TERM_DEBUG Inspect pairwise SOLR error terms and thru targets.

pairs = config.reference_pairs;
debug = struct();
noteLines = {};
noteLines{end + 1} = 'Phase II Error-Term and Thru-Target Debug Summary';
noteLines{end + 1} = '===============================================';
noteLines{end + 1} = '';

for idxPair = 1:numel(pairs)
    pairKey = upper(pairs(idxPair).pair);
    geomKey = upper(pairs(idxPair).geometry);
    pairLabel = sprintf('%s %s', pairKey, geomKey);

    fig = figure('Visible', config.figure_visible, 'Color', 'w');
    tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    hold on;
    title(sprintf('%s |t| terms', pairLabel));
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    grid on;
    xline(67, ':', 'Color', [0.7 0.7 0.7]);
    xline(115, ':', 'Color', [0.7 0.7 0.7]);

    nexttile;
    hold on;
    title(sprintf('%s \x2220 t_{21}, t_{12}', pairLabel));
    xlabel('Frequency (GHz)');
    ylabel('Phase (deg)');
    grid on;
    xline(67, ':', 'Color', [0.7 0.7 0.7]);
    xline(115, ':', 'Color', [0.7 0.7 0.7]);

    nexttile;
    hold on;
    title(sprintf('%s correction denominator', pairLabel));
    xlabel('Frequency (GHz)');
    ylabel('|den|');
    grid on;
    xline(67, ':', 'Color', [0.7 0.7 0.7]);
    xline(115, ':', 'Color', [0.7 0.7 0.7]);

    nexttile;
    hold on;
    title(sprintf('%s S_{21} vs Phase I target', pairLabel));
    xlabel('Frequency (GHz)');
    ylabel('|S_{21}| (dB)');
    grid on;
    xline(67, ':', 'Color', [0.7 0.7 0.7]);
    xline(115, ':', 'Color', [0.7 0.7 0.7]);

    nexttile;
    hold on;
    title(sprintf('%s S_{12} vs Phase I target', pairLabel));
    xlabel('Frequency (GHz)');
    ylabel('|S_{12}| (dB)');
    grid on;
    xline(67, ':', 'Color', [0.7 0.7 0.7]);
    xline(115, ':', 'Color', [0.7 0.7 0.7]);

    nexttile;
    hold on;
    title(sprintf('%s return loss vs Phase I target', pairLabel));
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    grid on;
    xline(67, ':', 'Color', [0.7 0.7 0.7]);
    xline(115, ':', 'Color', [0.7 0.7 0.7]);

    proxyT = gobjects(1, 4);
    proxyT(1) = plot(nan, nan, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
    proxyT(2) = plot(nan, nan, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
    proxyT(3) = plot(nan, nan, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.5);
    proxyT(4) = plot(nan, nan, '-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.5);
    nexttile(1); legend(proxyT, {'|t_{11}|', '|t_{22}|', '|t_{21}|', '|t_{12}|'}, 'Location', 'best');

    proxyComp = gobjects(1, 4);
    proxyComp(1) = plot(nan, nan, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
    proxyComp(2) = plot(nan, nan, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
    proxyComp(3) = plot(nan, nan, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8);
    proxyComp(4) = plot(nan, nan, ':', 'Color', [0 0 0], 'LineWidth', 1.8);
    nexttile(4); legend(proxyComp, {'Raw', 'Switch-corrected', 'SOLR-corrected', 'Phase I target'}, 'Location', 'best');
    nexttile(5); legend(proxyComp, {'Raw', 'Switch-corrected', 'SOLR-corrected', 'Phase I target'}, 'Location', 'best');
    nexttile(6); legend(proxyComp, {'S_{11} corrected', 'S_{11} target', 'S_{22} corrected', 'S_{22} target'}, 'Location', 'best');

    pairDebug = struct('pair', pairKey, 'geometry', geomKey, 'bands', {cell(1, numel(bandResults))});

    for idxBand = 1:numel(bandResults)
        pairResult = find_pair_result(bandResults{idxBand}, pairKey, geomKey);
        if isempty(pairResult)
            continue;
        end

        phase1Target = load_phase1_thru_target(config, geomKey, pairKey, bandResults{idxBand}.band);
        freq = pairResult.freq(:);
        freqGHz = freq / 1e9;
        e = pairResult.error_terms;

        nexttile(1);
        plot(freqGHz, local_mag_db(e.t11), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
        plot(freqGHz, local_mag_db(e.t22), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
        plot(freqGHz, local_mag_db(e.t21), '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.4);
        plot(freqGHz, local_mag_db(e.t12), '-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.4);

        nexttile(2);
        plot(freqGHz, unwrap(angle(e.t21)) * 180 / pi, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.4);
        plot(freqGHz, unwrap(angle(e.t12)) * 180 / pi, '-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.4);

        denVals = compute_correction_denominator(pairResult.reference_switch_corrected, e, config.cal_den_floor);
        nexttile(3);
        plot(freqGHz, denVals, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.4);

        rawS21 = squeeze(pairResult.reference_raw(2, 1, :));
        swS21 = squeeze(pairResult.reference_switch_corrected(2, 1, :));
        corrS21 = squeeze(pairResult.reference_corrected(2, 1, :));
        rawS12 = squeeze(pairResult.reference_raw(1, 2, :));
        swS12 = squeeze(pairResult.reference_switch_corrected(1, 2, :));
        corrS12 = squeeze(pairResult.reference_corrected(1, 2, :));

        nexttile(4);
        plot(freqGHz, local_mag_db(rawS21), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
        plot(freqGHz, local_mag_db(swS21), '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
        plot(freqGHz, local_mag_db(corrS21), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8);
        plot(freqGHz, local_mag_db(squeeze(phase1Target.S(2, 1, :))), ':', 'Color', [0 0 0], 'LineWidth', 1.8);

        nexttile(5);
        plot(freqGHz, local_mag_db(rawS12), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
        plot(freqGHz, local_mag_db(swS12), '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
        plot(freqGHz, local_mag_db(corrS12), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8);
        plot(freqGHz, local_mag_db(squeeze(phase1Target.S(1, 2, :))), ':', 'Color', [0 0 0], 'LineWidth', 1.8);

        nexttile(6);
        plot(freqGHz, local_mag_db(squeeze(pairResult.reference_corrected(1, 1, :))), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
        plot(freqGHz, local_mag_db(squeeze(phase1Target.S(1, 1, :))), '--', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
        plot(freqGHz, local_mag_db(squeeze(pairResult.reference_corrected(2, 2, :))), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
        plot(freqGHz, local_mag_db(squeeze(phase1Target.S(2, 2, :))), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);

        bandDebug = struct();
        bandDebug.band = bandResults{idxBand}.band;
        bandDebug.freq = freq;
        bandDebug.error_terms = e;
        bandDebug.denominator = denVals;
        bandDebug.phase1_target = phase1Target;
        bandDebug.reference_raw = pairResult.reference_raw;
        bandDebug.reference_switch_corrected = pairResult.reference_switch_corrected;
        bandDebug.reference_corrected = pairResult.reference_corrected;
        bandDebug.metrics = struct( ...
            'median_target_s21_db', median(local_mag_db(squeeze(phase1Target.S(2, 1, :))), 'omitnan'), ...
            'median_corr_s21_db', median(local_mag_db(corrS21), 'omitnan'), ...
            'median_target_s12_db', median(local_mag_db(squeeze(phase1Target.S(1, 2, :))), 'omitnan'), ...
            'median_corr_s12_db', median(local_mag_db(corrS12), 'omitnan'), ...
            'median_t21_db', median(local_mag_db(e.t21), 'omitnan'), ...
            'median_t12_db', median(local_mag_db(e.t12), 'omitnan'), ...
            'median_t11_db', median(local_mag_db(e.t11), 'omitnan'), ...
            'median_t22_db', median(local_mag_db(e.t22), 'omitnan'), ...
            'min_den', min(denVals, [], 'omitnan'), ...
            'median_den', median(denVals, 'omitnan'), ...
            'anchor_ratio_db', pairResult.transmission_anchor.anchor_ratio_db, ...
            'anchor_ratio_deg', pairResult.transmission_anchor.anchor_ratio_deg);
        pairDebug.bands{idxBand} = bandDebug;

        noteLines{end + 1} = sprintf('%s | %s', pairLabel, bandResults{idxBand}.band);
        noteLines{end + 1} = sprintf('  Median target/corrected |S21| (dB): %.3f / %.3f', ...
            bandDebug.metrics.median_target_s21_db, bandDebug.metrics.median_corr_s21_db);
        noteLines{end + 1} = sprintf('  Median target/corrected |S12| (dB): %.3f / %.3f', ...
            bandDebug.metrics.median_target_s12_db, bandDebug.metrics.median_corr_s12_db);
        noteLines{end + 1} = sprintf('  Median |t11|, |t22|, |t21|, |t12| (dB): %.3f, %.3f, %.3f, %.3f', ...
            bandDebug.metrics.median_t11_db, bandDebug.metrics.median_t22_db, ...
            bandDebug.metrics.median_t21_db, bandDebug.metrics.median_t12_db);
        noteLines{end + 1} = sprintf('  Transmission-anchor ratio (dB / deg): %.3f / %.3f', ...
            bandDebug.metrics.anchor_ratio_db, bandDebug.metrics.anchor_ratio_deg);
        noteLines{end + 1} = sprintf('  Min / median correction denominator magnitude: %.4e / %.4e', ...
            bandDebug.metrics.min_den, bandDebug.metrics.median_den);
    end

    sgtitle(fig, sprintf('Phase II error-term debug - %s', pairLabel));
    save_debug_figure(fig, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_ErrorTerms_%s_%s', pairKey, geomKey)));

    debug.(pairKey).(geomKey) = pairDebug; %#ok<STRNU>
    noteLines{end + 1} = '';
end

save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_ERROR_TERMS.mat'), 'debug');
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_ErrorTerms_Summary.txt'), noteLines);
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

function target = load_phase1_thru_target(config, geomKey, pairKey, bandName)
switch upper(geomKey)
    case 'STRAIGHT'
        if strcmpi(pairKey, 'P3P4')
            baseName = sprintf('EXTRACTED_THRU_P3P4_STRAIGHT_%s.mat', bandName);
        else
            baseName = sprintf('EXTRACTED_THRU_P1P2_STRAIGHT_%s.mat', bandName);
        end
    case 'ARC'
        baseName = sprintf('EXTRACTED_THRU_P1P4_ARC_%s.mat', bandName);
    case {'DIAG', 'DIAGONAL'}
        baseName = sprintf('EXTRACTED_THRU_P1P4_DIAGONAL_%s.mat', bandName);
    otherwise
        error('Unsupported geometry %s for Phase I target lookup.', geomKey);
end

loaded = load(fullfile(config.phase1_output_mat_dir, baseName));
target = loaded.reciprocalResult;
end

function denVals = compute_correction_denominator(Sraw, errorTerms, denominatorFloor)
nFreq = size(Sraw, 3);
denVals = nan(nFreq, 1);
for idx = 1:nFreq
    e1_00 = errorTerms.e1_00(idx);
    e1_11 = errorTerms.e1_11(idx);
    t11 = errorTerms.t11(idx);
    e2_00 = errorTerms.e2_00(idx);
    e2_11 = errorTerms.e2_11(idx);
    t22 = errorTerms.t22(idx);
    t21 = errorTerms.t21(idx);
    t12 = errorTerms.t12(idx);
    if any(~isfinite([e1_00, e1_11, t11, e2_00, e2_11, t22, t21, t12])) || ...
            any(abs([t11, t22, t21, t12]) < denominatorFloor)
        continue;
    end
    S11m = Sraw(1, 1, idx);
    S12m = Sraw(1, 2, idx);
    S21m = Sraw(2, 1, idx);
    S22m = Sraw(2, 2, idx);
    term1 = (S11m - e1_00) / t11;
    term2 = (S22m - e2_00) / t22;
    den = (1 + term1 * e1_11) * (1 + term2 * e2_11) - (S21m / t21) * (S12m / t12) * e1_11 * e2_11;
    denVals(idx) = abs(den);
end
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
