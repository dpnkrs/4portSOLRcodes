function generate_phase2_midband_gain_decomposition_debug(config, bandResults)
%GENERATE_PHASE2_MIDBAND_GAIN_DECOMPOSITION_DEBUG Decompose 67-115 corrected gain into switch, t, and denominator terms.

bandIdx = find(strcmpi(config.bands, '67-115'), 1, 'first');
if isempty(bandIdx)
    error('67-115 band not found in config.bands');
end
bandResult = bandResults{bandIdx};

pairs = bandResult.reference_pair_results;
debug = struct();
noteLines = {};
noteLines{end + 1} = 'Phase II 67-115 GHz gain decomposition audit';
noteLines{end + 1} = '==========================================';
noteLines{end + 1} = '';
noteLines{end + 1} = 'Convention:';
noteLines{end + 1} = '  |S21_corr|_dB = |S21_switch|_dB - |t21|_dB - |den|_dB';
noteLines{end + 1} = '  |S12_corr|_dB = |S12_switch|_dB - |t12|_dB - |den|_dB';
noteLines{end + 1} = '';

fig = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile; hold(ax1, 'on'); grid(ax1, 'on');
title(ax1, '67-115 GHz S_{21}: switch, -|t_{21}|, -|den|, corrected');
xlabel(ax1, 'Frequency (GHz)'); ylabel(ax1, 'Magnitude (dB)');

ax2 = nexttile; hold(ax2, 'on'); grid(ax2, 'on');
title(ax2, '67-115 GHz S_{12}: switch, -|t_{12}|, -|den|, corrected');
xlabel(ax2, 'Frequency (GHz)'); ylabel(ax2, 'Magnitude (dB)');

ax3 = nexttile; hold(ax3, 'on'); grid(ax3, 'on');
title(ax3, '67-115 GHz corrected minus switch contribution');
xlabel(ax3, 'Frequency (GHz)'); ylabel(ax3, 'Delta (dB)');
yline(ax3, 0, '--', 'Color', [0.6 0.6 0.6]);

ax4 = nexttile; hold(ax4, 'on'); grid(ax4, 'on');
title(ax4, '67-115 GHz denominator contribution');
xlabel(ax4, 'Frequency (GHz)'); ylabel(ax4, '-|den| (dB)');
yline(ax4, 0, '--', 'Color', [0.6 0.6 0.6]);

ax5 = nexttile; hold(ax5, 'on'); grid(ax5, 'on');
title(ax5, '67-115 GHz port participation: t_{11}, t_{22}');
xlabel(ax5, 'Frequency (GHz)'); ylabel(ax5, 'Magnitude (dB)');

ax6 = nexttile; hold(ax6, 'on'); grid(ax6, 'on');
title(ax6, '67-115 GHz predicted vs actual corrected |S21|');
xlabel(ax6, 'Frequency (GHz)'); ylabel(ax6, 'Magnitude (dB)');

colors = lines(numel(pairs));

for idxPair = 1:numel(pairs)
    pr = pairs(idxPair);
    freqGHz = pr.freq(:) / 1e9;
    s21sw = squeeze(pr.reference_switch_corrected(2, 1, :));
    s12sw = squeeze(pr.reference_switch_corrected(1, 2, :));
    s21corr = squeeze(pr.reference_corrected(2, 1, :));
    s12corr = squeeze(pr.reference_corrected(1, 2, :));
    t21 = pr.error_terms.t21(:);
    t12 = pr.error_terms.t12(:);
    t11 = pr.error_terms.t11(:);
    t22 = pr.error_terms.t22(:);
    den = compute_correction_denominator(pr.reference_switch_corrected, pr.error_terms, config.cal_den_floor);

    s21swDb = mag_db(s21sw);
    s12swDb = mag_db(s12sw);
    s21corrDb = mag_db(s21corr);
    s12corrDb = mag_db(s12corr);
    t21Db = mag_db(t21);
    t12Db = mag_db(t12);
    t11Db = mag_db(t11);
    t22Db = mag_db(t22);
    denDb = mag_db(den);
    predicted21Db = s21swDb - t21Db - denDb;
    predicted12Db = s12swDb - t12Db - denDb;

    c = colors(idxPair, :);
    label = sprintf('%s %s', pr.pair, pr.geometry);

    plot(ax1, freqGHz, s21swDb, '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax1, freqGHz, -t21Db, ':', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax1, freqGHz, -denDb, '-.', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax1, freqGHz, s21corrDb, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', label);

    plot(ax2, freqGHz, s12swDb, '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax2, freqGHz, -t12Db, ':', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax2, freqGHz, -denDb, '-.', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax2, freqGHz, s12corrDb, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', label);

    plot(ax3, freqGHz, s21corrDb - s21swDb, '-', 'Color', c, 'LineWidth', 1.4, 'DisplayName', label);
    plot(ax4, freqGHz, -denDb, '-', 'Color', c, 'LineWidth', 1.4, 'DisplayName', label);
    plot(ax5, freqGHz, t11Db, '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax5, freqGHz, t22Db, '-', 'Color', c, 'LineWidth', 1.4, 'DisplayName', label);
    plot(ax6, freqGHz, predicted21Db, '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax6, freqGHz, s21corrDb, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', label);

    bandDebug = struct();
    bandDebug.pair = pr.pair;
    bandDebug.geometry = pr.geometry;
    bandDebug.freq = pr.freq(:);
    bandDebug.metrics = struct( ...
        'median_switch_s21_db', median(s21swDb, 'omitnan'), ...
        'median_corrected_s21_db', median(s21corrDb, 'omitnan'), ...
        'median_switch_s12_db', median(s12swDb, 'omitnan'), ...
        'median_corrected_s12_db', median(s12corrDb, 'omitnan'), ...
        'median_minus_t21_db', median(-t21Db, 'omitnan'), ...
        'median_minus_t12_db', median(-t12Db, 'omitnan'), ...
        'median_minus_den_db', median(-denDb, 'omitnan'), ...
        'median_gain_change_s21_db', median(s21corrDb - s21swDb, 'omitnan'), ...
        'median_gain_change_s12_db', median(s12corrDb - s12swDb, 'omitnan'), ...
        'median_t11_db', median(t11Db, 'omitnan'), ...
        'median_t22_db', median(t22Db, 'omitnan'));
    debug.(pr.pair).(pr.geometry) = bandDebug; %#ok<STRNU>

    noteLines{end + 1} = label;
    noteLines{end + 1} = sprintf('  Median switch/corrected |S21| (dB): %.3f / %.3f', ...
        bandDebug.metrics.median_switch_s21_db, bandDebug.metrics.median_corrected_s21_db);
    noteLines{end + 1} = sprintf('  Median switch/corrected |S12| (dB): %.3f / %.3f', ...
        bandDebug.metrics.median_switch_s12_db, bandDebug.metrics.median_corrected_s12_db);
    noteLines{end + 1} = sprintf('  Median correction gain change S21 / S12 (dB): %.3f / %.3f', ...
        bandDebug.metrics.median_gain_change_s21_db, bandDebug.metrics.median_gain_change_s12_db);
    noteLines{end + 1} = sprintf('  Median contributions -|t21| / -|t12| / -|den| (dB): %.3f / %.3f / %.3f', ...
        bandDebug.metrics.median_minus_t21_db, bandDebug.metrics.median_minus_t12_db, bandDebug.metrics.median_minus_den_db);
    noteLines{end + 1} = sprintf('  Median t11 / t22 (dB): %.3f / %.3f', ...
        bandDebug.metrics.median_t11_db, bandDebug.metrics.median_t22_db);
    noteLines{end + 1} = '';
end

legend(ax1, 'Location', 'bestoutside');
legend(ax2, 'Location', 'bestoutside');
legend(ax3, 'Location', 'bestoutside');
legend(ax4, 'Location', 'bestoutside');
legend(ax5, 'Location', 'bestoutside');
legend(ax6, 'Location', 'bestoutside');

sgtitle(fig, 'Phase II 67-115 GHz transmission-gain decomposition');
save_debug_figure(fig, fullfile(config.interim_figure_dir, 'PhaseII_Debug_67115_GainDecomposition'));
save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_67115_GAIN_DECOMPOSITION.mat'), 'debug');
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_67115_GainDecomposition_Summary.txt'), noteLines);
end

function denVals = compute_correction_denominator(Sraw, errorTerms, denominatorFloor)
denVals = nan(size(Sraw, 3), 1);
for idx = 1:size(Sraw, 3)
    S11m = Sraw(1, 1, idx);
    S12m = Sraw(1, 2, idx);
    S21m = Sraw(2, 1, idx);
    S22m = Sraw(2, 2, idx);
    term1 = (S11m - errorTerms.e1_00(idx)) / errorTerms.t11(idx);
    term2 = (S22m - errorTerms.e2_00(idx)) / errorTerms.t22(idx);
    den = (1 + term1 * errorTerms.e1_11(idx)) * (1 + term2 * errorTerms.e2_11(idx)) ...
        - (S21m / errorTerms.t21(idx)) * (S12m / errorTerms.t12(idx)) * errorTerms.e1_11(idx) * errorTerms.e2_11(idx);
    if abs(den) < denominatorFloor
        den = denominatorFloor * exp(1i * angle(den + (den == 0)));
    end
    denVals(idx) = abs(den);
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
