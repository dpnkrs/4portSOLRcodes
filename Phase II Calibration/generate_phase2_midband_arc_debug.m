function generate_phase2_midband_arc_debug(config, bandResults)
%GENERATE_PHASE2_MIDBAND_ARC_DEBUG Focused 67-115 GHz ARC-reference consistency debug.

bandIdx = find(strcmpi(config.bands, '67-115'), 1, 'first');
if isempty(bandIdx)
    error('67-115 band not found in config.bands');
end
bandResult = bandResults{bandIdx};
if isempty(bandResult)
    error('67-115 band result is empty.');
end

pairs = bandResult.reference_pair_results;
freq = bandResult.freq(:);
freqGHz = freq / 1e9;

straightKeys = {'P1P2_STRAIGHT', 'P3P4_STRAIGHT'};
arcKeys = {'P1P3_ARC', 'P1P4_ARC', 'P2P3_ARC', 'P2P4_ARC'};

pairMap = struct();
for idx = 1:numel(pairs)
    key = sprintf('%s_%s', upper(pairs(idx).pair), upper(pairs(idx).geometry));
    pairMap.(key) = pairs(idx);
end

fig1 = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig1, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile; hold(ax1, 'on');
title(ax1, '67-115 GHz switch-corrected reciprocal references vs Phase I targets: |S_{21}|');
xlabel(ax1, 'Frequency (GHz)'); ylabel(ax1, '|S_{21}| (dB)'); grid(ax1, 'on');
ax2 = nexttile; hold(ax2, 'on');
title(ax2, '67-115 GHz switch-corrected reciprocal references vs Phase I targets: |S_{12}|');
xlabel(ax2, 'Frequency (GHz)'); ylabel(ax2, '|S_{12}| (dB)'); grid(ax2, 'on');

colors = lines(numel(straightKeys) + numel(arcKeys));
allKeys = [straightKeys, arcKeys];
for idx = 1:numel(allKeys)
    key = allKeys{idx};
    pr = pairMap.(key);
    sw21 = squeeze(pr.reference_switch_corrected(2,1,:));
    sw12 = squeeze(pr.reference_switch_corrected(1,2,:));
    tgt = load_phase1_target_for_debug(config, upper(pr.geometry), upper(pr.pair), '67-115');
    tgt21 = squeeze(tgt.S(2,1,:));
    tgt12 = squeeze(tgt.S(1,2,:));
    c = colors(idx,:);
    plot(ax1, freqGHz, mag_db(sw21), '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', [strrep(key,'_',' ') ' switch']);
    plot(ax1, freqGHz, mag_db(tgt21), '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax2, freqGHz, mag_db(sw12), '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', [strrep(key,'_',' ') ' switch']);
    plot(ax2, freqGHz, mag_db(tgt12), '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
end
legend(ax1, 'Location', 'bestoutside');
legend(ax2, 'Location', 'bestoutside');

fig2 = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig2, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax3 = nexttile; hold(ax3, 'on');
title(ax3, '67-115 GHz target mismatch after switch correction: |S_{21}|');
xlabel(ax3, 'Frequency (GHz)'); ylabel(ax3, 'Switch - target (dB)'); grid(ax3, 'on'); yline(ax3, 0, '--', 'Color', [0.6 0.6 0.6]);
ax4 = nexttile; hold(ax4, 'on');
title(ax4, '67-115 GHz target mismatch after switch correction: |S_{12}|');
xlabel(ax4, 'Frequency (GHz)'); ylabel(ax4, 'Switch - target (dB)'); grid(ax4, 'on'); yline(ax4, 0, '--', 'Color', [0.6 0.6 0.6]);

for idx = 1:numel(allKeys)
    key = allKeys{idx};
    pr = pairMap.(key);
    sw21 = squeeze(pr.reference_switch_corrected(2,1,:));
    sw12 = squeeze(pr.reference_switch_corrected(1,2,:));
    tgt = load_phase1_target_for_debug(config, upper(pr.geometry), upper(pr.pair), '67-115');
    tgt21 = squeeze(tgt.S(2,1,:));
    tgt12 = squeeze(tgt.S(1,2,:));
    c = colors(idx,:);
    plot(ax3, freqGHz, mag_db(sw21) - mag_db(tgt21), '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', strrep(key,'_',' '));
    plot(ax4, freqGHz, mag_db(sw12) - mag_db(tgt12), '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', strrep(key,'_',' '));
end
legend(ax3, 'Location', 'bestoutside');
legend(ax4, 'Location', 'bestoutside');

fig3 = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig3, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax5 = nexttile; hold(ax5, 'on');
title(ax5, '67-115 GHz ARC excess loss relative to straight-reference average: |S_{21}|');
xlabel(ax5, 'Frequency (GHz)'); ylabel(ax5, 'ARC - straight average (dB)'); grid(ax5, 'on');
yline(ax5, 0, '--', 'Color', [0.6 0.6 0.6]);
ax6 = nexttile; hold(ax6, 'on');
title(ax6, '67-115 GHz ARC excess loss relative to straight-reference average: |S_{12}|');
xlabel(ax6, 'Frequency (GHz)'); ylabel(ax6, 'ARC - straight average (dB)'); grid(ax6, 'on');
yline(ax6, 0, '--', 'Color', [0.6 0.6 0.6]);

straightSw21 = []; straightSw12 = []; straightTgt21 = []; straightTgt12 = [];
for idx = 1:numel(straightKeys)
    pr = pairMap.(straightKeys{idx});
    tgt = load_phase1_target_for_debug(config, upper(pr.geometry), upper(pr.pair), '67-115');
    straightSw21(:,idx) = mag_db(squeeze(pr.reference_switch_corrected(2,1,:))); %#ok<AGROW>
    straightSw12(:,idx) = mag_db(squeeze(pr.reference_switch_corrected(1,2,:))); %#ok<AGROW>
    straightTgt21(:,idx) = mag_db(squeeze(tgt.S(2,1,:))); %#ok<AGROW>
    straightTgt12(:,idx) = mag_db(squeeze(tgt.S(1,2,:))); %#ok<AGROW>
end
straightSw21Avg = mean(straightSw21, 2, 'omitnan');
straightSw12Avg = mean(straightSw12, 2, 'omitnan');
straightTgt21Avg = mean(straightTgt21, 2, 'omitnan');
straightTgt12Avg = mean(straightTgt12, 2, 'omitnan');

arcSummary = struct();
for idx = 1:numel(arcKeys)
    key = arcKeys{idx};
    pr = pairMap.(key);
    c = colors(numel(straightKeys) + idx,:);
    tgt = load_phase1_target_for_debug(config, upper(pr.geometry), upper(pr.pair), '67-115');
    arcSw21 = mag_db(squeeze(pr.reference_switch_corrected(2,1,:)));
    arcSw12 = mag_db(squeeze(pr.reference_switch_corrected(1,2,:)));
    arcTgt21 = mag_db(squeeze(tgt.S(2,1,:)));
    arcTgt12 = mag_db(squeeze(tgt.S(1,2,:)));

    deltaSw21 = arcSw21 - straightSw21Avg;
    deltaSw12 = arcSw12 - straightSw12Avg;
    deltaTgt21 = arcTgt21 - straightTgt21Avg;
    deltaTgt12 = arcTgt12 - straightTgt12Avg;

    plot(ax5, freqGHz, deltaSw21, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', [strrep(key,'_',' ') ' switch']);
    plot(ax5, freqGHz, deltaTgt21, '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
    plot(ax6, freqGHz, deltaSw12, '-', 'Color', c, 'LineWidth', 1.5, 'DisplayName', [strrep(key,'_',' ') ' switch']);
    plot(ax6, freqGHz, deltaTgt12, '--', 'Color', c, 'LineWidth', 1.0, 'HandleVisibility', 'off');

    arcSummary.(key) = struct( ...
        'median_switch_target_s21_delta_db', median(arcSw21 - arcTgt21, 'omitnan'), ...
        'median_switch_target_s12_delta_db', median(arcSw12 - arcTgt12, 'omitnan'), ...
        'median_switch_minus_straightavg_s21_db', median(deltaSw21, 'omitnan'), ...
        'median_target_minus_straightavg_s21_db', median(deltaTgt21, 'omitnan'), ...
        'median_switch_minus_straightavg_s12_db', median(deltaSw12, 'omitnan'), ...
        'median_target_minus_straightavg_s12_db', median(deltaTgt12, 'omitnan'));
end
legend(ax5, 'Location', 'bestoutside');
legend(ax6, 'Location', 'bestoutside');

save_debug_figure(fig1, fullfile(config.interim_figure_dir, 'PhaseII_Debug_67115_RefComparison'));
save_debug_figure(fig2, fullfile(config.interim_figure_dir, 'PhaseII_Debug_67115_TargetMismatch'));
save_debug_figure(fig3, fullfile(config.interim_figure_dir, 'PhaseII_Debug_67115_ARCvsStraight'));

debug = struct();
debug.band = '67-115';
debug.freq = freq;
debug.arcSummary = arcSummary;
debug.straightAverage = struct('switchS21dB', straightSw21Avg, 'switchS12dB', straightSw12Avg, 'targetS21dB', straightTgt21Avg, 'targetS12dB', straightTgt12Avg);
save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_67115_ARC_REFERENCE_AUDIT.mat'), 'debug');

noteLines = {};
noteLines{end+1} = 'Phase II 67-115 GHz ARC-reference audit';
noteLines{end+1} = '====================================';
noteLines{end+1} = '';
for idx = 1:numel(arcKeys)
    key = arcKeys{idx};
    m = arcSummary.(key);
    noteLines{end+1} = strrep(key, '_', ' ');
    noteLines{end+1} = sprintf('  Median switch-target delta |S21| / |S12| (dB): %.3f / %.3f', m.median_switch_target_s21_delta_db, m.median_switch_target_s12_delta_db);
    noteLines{end+1} = sprintf('  Median switch minus straight-average |S21| / |S12| (dB): %.3f / %.3f', m.median_switch_minus_straightavg_s21_db, m.median_switch_minus_straightavg_s12_db);
    noteLines{end+1} = sprintf('  Median target minus straight-average |S21| / |S12| (dB): %.3f / %.3f', m.median_target_minus_straightavg_s21_db, m.median_target_minus_straightavg_s12_db);
    noteLines{end+1} = '';
end
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_67115_ARC_ReferenceAudit_Summary.txt'), noteLines);
end

function values = mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function target = load_phase1_target_for_debug(config, geomKey, pairKey, bandName)
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
