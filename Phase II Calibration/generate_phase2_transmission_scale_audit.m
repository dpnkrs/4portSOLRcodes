function generate_phase2_transmission_scale_audit(config, bandResults, referenceStandards)
%GENERATE_PHASE2_TRANSMISSION_SCALE_AUDIT Compare local and legacy-form pair solves.

pairs = config.reference_pairs;
debug = struct();
noteLines = {};
noteLines{end + 1} = 'Phase II Transmission-Scale Audit Summary';
noteLines{end + 1} = '======================================';
noteLines{end + 1} = '';

for idxPair = 1:numel(pairs)
    pairKey = upper(pairs(idxPair).pair);
    geomKey = upper(pairs(idxPair).geometry);
    pairLabel = sprintf('%s %s', pairKey, geomKey);

    fig = figure('Visible', config.figure_visible, 'Color', 'w');
    tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile;
    hold(ax1, 'on');
    title(ax1, sprintf('%s S_{21} transmission audit', pairLabel));
    xlabel(ax1, 'Frequency (GHz)');
    ylabel(ax1, '|S_{21}| (dB)');
    grid(ax1, 'on');
    xline(ax1, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax1, 115, ':', 'Color', [0.7 0.7 0.7]);

    ax2 = nexttile;
    hold(ax2, 'on');
    title(ax2, sprintf('%s S_{12} transmission audit', pairLabel));
    xlabel(ax2, 'Frequency (GHz)');
    ylabel(ax2, '|S_{12}| (dB)');
    grid(ax2, 'on');
    xline(ax2, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax2, 115, ':', 'Color', [0.7 0.7 0.7]);

    ax3 = nexttile;
    hold(ax3, 'on');
    title(ax3, sprintf('%s local vs legacy |t_{21}|, |t_{12}|', pairLabel));
    xlabel(ax3, 'Frequency (GHz)');
    ylabel(ax3, 'Magnitude (dB)');
    grid(ax3, 'on');
    xline(ax3, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax3, 115, ':', 'Color', [0.7 0.7 0.7]);

    ax4 = nexttile;
    hold(ax4, 'on');
    title(ax4, sprintf('%s local vs legacy |t_{11}|, |t_{22}|', pairLabel));
    xlabel(ax4, 'Frequency (GHz)');
    ylabel(ax4, 'Magnitude (dB)');
    grid(ax4, 'on');
    xline(ax4, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax4, 115, ':', 'Color', [0.7 0.7 0.7]);

    ax5 = nexttile;
    hold(ax5, 'on');
    title(ax5, sprintf('%s corrected-target mismatch', pairLabel));
    xlabel(ax5, 'Frequency (GHz)');
    ylabel(ax5, 'Delta (dB)');
    grid(ax5, 'on');
    xline(ax5, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax5, 115, ':', 'Color', [0.7 0.7 0.7]);

    ax6 = nexttile;
    hold(ax6, 'on');
    title(ax6, sprintf('%s local vs legacy corrected delta', pairLabel));
    xlabel(ax6, 'Frequency (GHz)');
    ylabel(ax6, 'Delta (dB)');
    grid(ax6, 'on');
    xline(ax6, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax6, 115, ':', 'Color', [0.7 0.7 0.7]);

    proxyTx = gobjects(1, 4);
    proxyTx(1) = plot(ax1, nan, nan, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
    proxyTx(2) = plot(ax1, nan, nan, '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
    proxyTx(3) = plot(ax1, nan, nan, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8);
    proxyTx(4) = plot(ax1, nan, nan, ':', 'Color', [0 0 0], 'LineWidth', 1.8);
    legend(ax1, proxyTx, {'Switch-corrected', 'Legacy-form corrected', 'Local corrected', 'Phase I target'}, 'Location', 'best');
    legend(ax2, proxyTx, {'Switch-corrected', 'Legacy-form corrected', 'Local corrected', 'Phase I target'}, 'Location', 'best');

    proxyTerms = gobjects(1, 4);
    proxyTerms(1) = plot(ax3, nan, nan, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
    proxyTerms(2) = plot(ax3, nan, nan, '--', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
    proxyTerms(3) = plot(ax3, nan, nan, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
    proxyTerms(4) = plot(ax3, nan, nan, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
    legend(ax3, proxyTerms, {'Local |t_{21}|', 'Legacy |t_{21}|', 'Local |t_{12}|', 'Legacy |t_{12}|'}, 'Location', 'best');

    proxyRefl = gobjects(1, 4);
    proxyRefl(1) = plot(ax4, nan, nan, '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.5);
    proxyRefl(2) = plot(ax4, nan, nan, '--', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.2);
    proxyRefl(3) = plot(ax4, nan, nan, '-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.5);
    proxyRefl(4) = plot(ax4, nan, nan, '--', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);
    legend(ax4, proxyRefl, {'Local |t_{11}|', 'Legacy |t_{11}|', 'Local |t_{22}|', 'Legacy |t_{22}|'}, 'Location', 'best');

    proxyMismatch = gobjects(1, 4);
    proxyMismatch(1) = plot(ax5, nan, nan, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
    proxyMismatch(2) = plot(ax5, nan, nan, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
    proxyMismatch(3) = plot(ax5, nan, nan, '-.', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.5);
    proxyMismatch(4) = plot(ax5, nan, nan, ':', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);
    legend(ax5, proxyMismatch, {'Local S_{21}-target', 'Legacy S_{21}-target', 'Local S_{12}-target', 'Legacy S_{12}-target'}, 'Location', 'best');
    legend(ax6, proxyMismatch, {'S_{21} local-legacy', 'S_{12} local-legacy', '|t_{21}| local-legacy', '|t_{11}| local-legacy'}, 'Location', 'best');

    pairDebug = struct('pair', pairKey, 'geometry', geomKey, 'bands', {cell(1, numel(bandResults))});

    for idxBand = 1:numel(bandResults)
        pairResult = find_pair_result(bandResults{idxBand}, pairKey, geomKey);
        if isempty(pairResult)
            continue;
        end

        freq = pairResult.freq(:);
        freqGHz = freq / 1e9;
        gammaShort = interpolate_reference_gamma(referenceStandards.Short, freq);
        gammaOpen = interpolate_reference_gamma(referenceStandards.Open, freq);
        gammaLoad = interpolate_reference_gamma(referenceStandards.Load, freq);
        legacyTerms = calculate_error_terms_solr_legacy_form_from_gamma(freq, ...
            pairResult.standards_measured_port1.Short(:), pairResult.standards_measured_port1.Open(:), pairResult.standards_measured_port1.Load(:), ...
            pairResult.standards_measured_port2.Short(:), pairResult.standards_measured_port2.Open(:), pairResult.standards_measured_port2.Load(:), ...
            pairResult.reference_switch_corrected, gammaShort, gammaOpen, gammaLoad, 0);
        legacyCorrected = correct_network_from_error_terms(pairResult.reference_switch_corrected, legacyTerms, config.cal_den_floor);

        targetS21 = squeeze(pairResult.reference_phase1_target.S(2, 1, :));
        targetS12 = squeeze(pairResult.reference_phase1_target.S(1, 2, :));
        switchS21 = squeeze(pairResult.reference_switch_corrected(2, 1, :));
        switchS12 = squeeze(pairResult.reference_switch_corrected(1, 2, :));
        localS21 = squeeze(pairResult.reference_corrected(2, 1, :));
        localS12 = squeeze(pairResult.reference_corrected(1, 2, :));
        legacyS21 = squeeze(legacyCorrected(2, 1, :));
        legacyS12 = squeeze(legacyCorrected(1, 2, :));

        plot(ax1, freqGHz, local_mag_db(switchS21), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
        plot(ax1, freqGHz, local_mag_db(legacyS21), '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
        plot(ax1, freqGHz, local_mag_db(localS21), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8);
        plot(ax1, freqGHz, local_mag_db(targetS21), ':', 'Color', [0 0 0], 'LineWidth', 1.8);

        plot(ax2, freqGHz, local_mag_db(switchS12), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
        plot(ax2, freqGHz, local_mag_db(legacyS12), '-.', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
        plot(ax2, freqGHz, local_mag_db(localS12), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.8);
        plot(ax2, freqGHz, local_mag_db(targetS12), ':', 'Color', [0 0 0], 'LineWidth', 1.8);

        plot(ax3, freqGHz, local_mag_db(pairResult.error_terms.t21), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
        plot(ax3, freqGHz, local_mag_db(legacyTerms.t21), '--', 'Color', [0 0.45 0.74], 'LineWidth', 1.2);
        plot(ax3, freqGHz, local_mag_db(pairResult.error_terms.t12), '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
        plot(ax3, freqGHz, local_mag_db(legacyTerms.t12), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);

        plot(ax4, freqGHz, local_mag_db(pairResult.error_terms.t11), '-', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.5);
        plot(ax4, freqGHz, local_mag_db(legacyTerms.t11), '--', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.2);
        plot(ax4, freqGHz, local_mag_db(pairResult.error_terms.t22), '-', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.5);
        plot(ax4, freqGHz, local_mag_db(legacyTerms.t22), '--', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);

        localMismatch21 = local_mag_db(localS21) - local_mag_db(targetS21);
        legacyMismatch21 = local_mag_db(legacyS21) - local_mag_db(targetS21);
        localMismatch12 = local_mag_db(localS12) - local_mag_db(targetS12);
        legacyMismatch12 = local_mag_db(legacyS12) - local_mag_db(targetS12);
        plot(ax5, freqGHz, localMismatch21, '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
        plot(ax5, freqGHz, legacyMismatch21, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
        plot(ax5, freqGHz, localMismatch12, '-.', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.5);
        plot(ax5, freqGHz, legacyMismatch12, ':', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);

        plot(ax6, freqGHz, local_mag_db(localS21) - local_mag_db(legacyS21), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
        plot(ax6, freqGHz, local_mag_db(localS12) - local_mag_db(legacyS12), '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.2);
        plot(ax6, freqGHz, local_mag_db(pairResult.error_terms.t21) - local_mag_db(legacyTerms.t21), '-.', 'Color', [0.47 0.67 0.19], 'LineWidth', 1.5);
        plot(ax6, freqGHz, local_mag_db(pairResult.error_terms.t11) - local_mag_db(legacyTerms.t11), ':', 'Color', [0.49 0.18 0.56], 'LineWidth', 1.2);

        bandDebug = struct();
        bandDebug.band = bandResults{idxBand}.band;
        bandDebug.freq = freq;
        bandDebug.local_terms = pairResult.error_terms;
        bandDebug.legacy_terms = legacyTerms;
        bandDebug.local_corrected = pairResult.reference_corrected;
        bandDebug.legacy_corrected = legacyCorrected;
        bandDebug.target = pairResult.reference_phase1_target;
        bandDebug.metrics = struct( ...
            'median_switch_s21_db', median(local_mag_db(switchS21), 'omitnan'), ...
            'median_target_s21_db', median(local_mag_db(targetS21), 'omitnan'), ...
            'median_local_s21_db', median(local_mag_db(localS21), 'omitnan'), ...
            'median_legacy_s21_db', median(local_mag_db(legacyS21), 'omitnan'), ...
            'median_switch_s12_db', median(local_mag_db(switchS12), 'omitnan'), ...
            'median_target_s12_db', median(local_mag_db(targetS12), 'omitnan'), ...
            'median_local_s12_db', median(local_mag_db(localS12), 'omitnan'), ...
            'median_legacy_s12_db', median(local_mag_db(legacyS12), 'omitnan'), ...
            'median_local_legacy_t21_db_delta', median(local_mag_db(pairResult.error_terms.t21) - local_mag_db(legacyTerms.t21), 'omitnan'), ...
            'median_local_legacy_t11_db_delta', median(local_mag_db(pairResult.error_terms.t11) - local_mag_db(legacyTerms.t11), 'omitnan'), ...
            'median_local_legacy_s21_db_delta', median(local_mag_db(localS21) - local_mag_db(legacyS21), 'omitnan'), ...
            'median_local_legacy_s12_db_delta', median(local_mag_db(localS12) - local_mag_db(legacyS12), 'omitnan'));
        pairDebug.bands{idxBand} = bandDebug;

        noteLines{end + 1} = sprintf('%s | %s', pairLabel, bandResults{idxBand}.band);
        noteLines{end + 1} = sprintf('  Median target / switch / local / legacy |S21| (dB): %.3f / %.3f / %.3f / %.3f', ...
            bandDebug.metrics.median_target_s21_db, bandDebug.metrics.median_switch_s21_db, ...
            bandDebug.metrics.median_local_s21_db, bandDebug.metrics.median_legacy_s21_db);
        noteLines{end + 1} = sprintf('  Median target / switch / local / legacy |S12| (dB): %.3f / %.3f / %.3f / %.3f', ...
            bandDebug.metrics.median_target_s12_db, bandDebug.metrics.median_switch_s12_db, ...
            bandDebug.metrics.median_local_s12_db, bandDebug.metrics.median_legacy_s12_db);
        noteLines{end + 1} = sprintf('  Median local-legacy delta |t21| / |t11| (dB): %.4f / %.4f', ...
            bandDebug.metrics.median_local_legacy_t21_db_delta, bandDebug.metrics.median_local_legacy_t11_db_delta);
        noteLines{end + 1} = sprintf('  Median local-legacy delta corrected |S21| / |S12| (dB): %.4f / %.4f', ...
            bandDebug.metrics.median_local_legacy_s21_db_delta, bandDebug.metrics.median_local_legacy_s12_db_delta);
        noteLines{end + 1} = '';
    end

    sgtitle(fig, sprintf('Phase II transmission-scale audit - %s', pairLabel));
    save_debug_figure(fig, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_TransmissionScaleAudit_%s_%s', pairKey, geomKey)));
    debug.(pairKey).(geomKey) = pairDebug; %#ok<STRNU>
end

save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_TRANSMISSION_SCALE_AUDIT.mat'), 'debug');
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_TransmissionScaleAudit_Summary.txt'), noteLines);
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

function gamma = interpolate_reference_gamma(reference, freq)
gamma = interp1(reference.freq, reference.gamma, freq, 'pchip', 'extrap');
end

function corrected = correct_network_from_error_terms(Sraw, errorTerms, denominatorFloor)
corrected = nan(size(Sraw));
for idxFreq = 1:size(Sraw, 3)
    corrected(:, :, idxFreq) = apply_error_correction_8term_local(Sraw(:, :, idxFreq), errorTerms, idxFreq, denominatorFloor);
end
end

function errorTerms = calculate_error_terms_solr_legacy_form_from_gamma(freq, measShortP1, measOpenP1, measLoadP1, measShortP2, measOpenP2, measLoadP2, reciprocalMeas, gammaShort, gammaOpen, gammaLoad, thruDelay)
% Legacy-form reproduction of SOLR Cal\calculate_error_terms_solr_NonIdeal_PhAmbg

if nargin < 12 || isempty(thruDelay)
    thruDelay = 0;
end

nFreq = numel(freq);
fields = {'e1_00', 'e1_11', 't11', 'e2_00', 'e2_11', 't22', 't21', 't12'};
for idx = 1:numel(fields)
    errorTerms.(fields{idx}) = nan(nFreq, 1);
end

for idx = 1:nFreq
    A1 = [ ...
        1, measShortP1(idx) * gammaShort(idx), gammaShort(idx); ...
        1, measOpenP1(idx) * gammaOpen(idx), gammaOpen(idx); ...
        1, measLoadP1(idx) * gammaLoad(idx), gammaLoad(idx)];
    b1 = [measShortP1(idx); measOpenP1(idx); measLoadP1(idx)];

    A2 = [ ...
        1, measShortP2(idx) * gammaShort(idx), gammaShort(idx); ...
        1, measOpenP2(idx) * gammaOpen(idx), gammaOpen(idx); ...
        1, measLoadP2(idx) * gammaLoad(idx), gammaLoad(idx)];
    b2 = [measShortP2(idx); measOpenP2(idx); measLoadP2(idx)];

    if rcond(A1) < 1e-12 || rcond(A2) < 1e-12
        continue;
    end

    sol1 = A1 \ b1;
    sol2 = A2 \ b2;
    errorTerms.e1_00(idx) = sol1(1);
    errorTerms.e1_11(idx) = sol1(2);
    errorTerms.t11(idx) = sol1(3) + sol1(1) * sol1(2);
    errorTerms.e2_00(idx) = sol2(1);
    errorTerms.e2_11(idx) = sol2(2);
    errorTerms.t22(idx) = sol2(3) + sol2(1) * sol2(2);

    s21m = reciprocalMeas(2, 1, idx);
    s12m = reciprocalMeas(1, 2, idx);
    if abs(s12m) < 1e-15
        continue;
    end

    t21sq = s21m * errorTerms.t11(idx) * errorTerms.t22(idx) / s12m;
    cand1 = sqrt(t21sq);
    cand2 = -cand1;

    expectedThruPhase = -2 * pi * freq(idx) * thruDelay;
    measuredThruPhase = angle(s21m);
    targetPhase = measuredThruPhase - expectedThruPhase;

    diff1 = abs(angle(exp(1i * (angle(cand1) - targetPhase))));
    diff2 = abs(angle(exp(1i * (angle(cand2) - targetPhase))));
    if diff1 <= diff2
        errorTerms.t21(idx) = cand1;
    else
        errorTerms.t21(idx) = cand2;
    end

    if abs(errorTerms.t21(idx)) >= 1e-15
        errorTerms.t12(idx) = (errorTerms.t11(idx) * errorTerms.t22(idx)) / errorTerms.t21(idx);
    end
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
