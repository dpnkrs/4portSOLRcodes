function audit = generate_phase2_reciprocal_scale_audit(config, bandResults)
%GENERATE_PHASE2_RECIPROCAL_SCALE_AUDIT Audit downstream-only reciprocal transmission scaling.

if nargin < 1 || isempty(config)
    config = phase2_config();
end
if nargin < 2 || isempty(bandResults)
    loaded = load(fullfile(config.mat_dir, 'PHASE2_RESULTS_ALL.mat'), 'results');
    bandResults = loaded.results.bands;
end

noteLines = {
    'Phase II Reciprocal Transmission-Scale Audit'
    '==========================================='
    ''
    'This audit is downstream-only: it does not use any Phase I thru target overlay.'
    'It decomposes the switch-to-corrected transmission lift into the solved tracking terms'
    'and groups the results by pair, geometry, and receive-port participation.'
    ''
    };

audit = struct();
allRows = struct([]);
rowCount = 0;

for idxBand = 1:numel(bandResults)
    bandResult = bandResults{idxBand};
    pairRows = struct([]);

    for idxPair = 1:numel(bandResult.reference_pair_results)
        pr = bandResult.reference_pair_results(idxPair);
        row = summarize_pair_scale_row(pr, config, bandResult.band);
        rowCount = rowCount + 1;
        if isempty(allRows)
            allRows = row;
        else
            allRows(rowCount) = row; %#ok<AGROW>
        end
        if isempty(pairRows)
            pairRows = row;
        else
            pairRows(end + 1) = row; %#ok<AGROW>
        end
    end

    audit.(matlab.lang.makeValidName(strrep(bandResult.band, '-', '_'))) = struct( ...
        'band', bandResult.band, ...
        'rows', pairRows, ...
        'receive_port_summary_s21', summarize_by_receive_port(pairRows, 's21'), ...
        'receive_port_summary_s12', summarize_by_receive_port(pairRows, 's12'), ...
        'geometry_summary', summarize_by_geometry(pairRows));

    noteLines = [noteLines; build_band_note_block(audit.(matlab.lang.makeValidName(strrep(bandResult.band, '-', '_')))); {''}]; %#ok<AGROW>
end

audit.all_rows = allRows;
save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_RECIPROCAL_SCALE_AUDIT.mat'), 'audit');

fig = create_reciprocal_scale_figure(config, allRows);
save_debug_figure(fig, fullfile(config.interim_figure_dir, 'PhaseII_Debug_ReciprocalScaleAudit'));

write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_ReciprocalScaleAudit_Summary.txt'), noteLines);
end

function row = summarize_pair_scale_row(pr, config, bandName)
freq = pr.freq(:);
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
minusT21Db = -mag_db(t21);
minusT12Db = -mag_db(t12);
minusDenDb = -mag_db(den);

row = struct();
row.band = bandName;
row.pair = pr.pair;
row.geometry = pr.geometry;
row.port1 = pr.ports(1);
row.port2 = pr.ports(2);
row.freq = freq;
row.median_switch_s21_db = median(s21swDb, 'omitnan');
row.median_switch_s12_db = median(s12swDb, 'omitnan');
row.median_corrected_s21_db = median(s21corrDb, 'omitnan');
row.median_corrected_s12_db = median(s12corrDb, 'omitnan');
row.median_gain_s21_db = median(s21corrDb - s21swDb, 'omitnan');
row.median_gain_s12_db = median(s12corrDb - s12swDb, 'omitnan');
row.median_minus_t21_db = median(minusT21Db, 'omitnan');
row.median_minus_t12_db = median(minusT12Db, 'omitnan');
row.median_minus_den_db = median(minusDenDb, 'omitnan');
row.median_t11_db = median(mag_db(t11), 'omitnan');
row.median_t22_db = median(mag_db(t22), 'omitnan');
row.median_return_radius = median(max(abs([squeeze(pr.reference_corrected(1, 1, :)), squeeze(pr.reference_corrected(2, 2, :))]), [], 2), 'omitnan');
row.above_zero_s21 = nnz(s21corrDb > 0);
row.above_zero_s12 = nnz(s12corrDb > 0);
end

function summary = summarize_by_receive_port(rows, mode)
summary = struct();
if strcmpi(mode, 's21')
    portField = 'port2';
    gainField = 'median_gain_s21_db';
    corrField = 'median_corrected_s21_db';
    switchField = 'median_switch_s21_db';
    aboveField = 'above_zero_s21';
else
    portField = 'port1';
    gainField = 'median_gain_s12_db';
    corrField = 'median_corrected_s12_db';
    switchField = 'median_switch_s12_db';
    aboveField = 'above_zero_s12';
end

for portNum = 1:4
    hits = rows(arrayfun(@(r) r.(portField) == portNum, rows));
    if isempty(hits)
        continue;
    end
    key = sprintf('P%d', portNum);
    summary.(key) = struct( ...
        'num_pairs', numel(hits), ...
        'median_switch_db', median([hits.(switchField)], 'omitnan'), ...
        'median_corrected_db', median([hits.(corrField)], 'omitnan'), ...
        'median_gain_db', median([hits.(gainField)], 'omitnan'), ...
        'median_minus_t_db', median(extract_minus_t(hits, mode), 'omitnan'), ...
        'median_minus_den_db', median([hits.median_minus_den_db], 'omitnan'), ...
        'total_above_zero_count', sum([hits.(aboveField)]));
end
end

function values = extract_minus_t(rows, mode)
if strcmpi(mode, 's21')
    values = [rows.median_minus_t21_db];
else
    values = [rows.median_minus_t12_db];
end
end

function summary = summarize_by_geometry(rows)
summary = struct();
geomList = unique({rows.geometry});
for idxGeom = 1:numel(geomList)
    geomKey = geomList{idxGeom};
    hits = rows(strcmpi({rows.geometry}, geomKey));
    summary.(geomKey) = struct( ...
        'num_pairs', numel(hits), ...
        'median_switch_s21_db', median([hits.median_switch_s21_db], 'omitnan'), ...
        'median_corrected_s21_db', median([hits.median_corrected_s21_db], 'omitnan'), ...
        'median_gain_s21_db', median([hits.median_gain_s21_db], 'omitnan'), ...
        'median_switch_s12_db', median([hits.median_switch_s12_db], 'omitnan'), ...
        'median_corrected_s12_db', median([hits.median_corrected_s12_db], 'omitnan'), ...
        'median_gain_s12_db', median([hits.median_gain_s12_db], 'omitnan'), ...
        'median_minus_t21_db', median([hits.median_minus_t21_db], 'omitnan'), ...
        'median_minus_t12_db', median([hits.median_minus_t12_db], 'omitnan'), ...
        'median_minus_den_db', median([hits.median_minus_den_db], 'omitnan'));
end
end

function lines = build_band_note_block(bandAudit)
lines = {
    sprintf('Band %s', bandAudit.band)
    '-------------'
    'Pair-level summary:'
    };

for idxRow = 1:numel(bandAudit.rows)
    row = bandAudit.rows(idxRow);
    lines{end + 1} = sprintf(['  %s %s: switch/corrected S21 %.3f / %.3f dB, ' ...
        'gain %.3f dB, -|t21| %.3f dB, -|den| %.3f dB, switch/corrected S12 %.3f / %.3f dB'], ... %#ok<AGROW>
        row.pair, row.geometry, row.median_switch_s21_db, row.median_corrected_s21_db, ...
        row.median_gain_s21_db, row.median_minus_t21_db, row.median_minus_den_db, ...
        row.median_switch_s12_db, row.median_corrected_s12_db);
    lines{end + 1} = sprintf(['    ports [%d %d], return radius %.3f, S21>S0 count %d, S12>S0 count %d'], ... %#ok<AGROW>
        row.port1, row.port2, row.median_return_radius, row.above_zero_s21, row.above_zero_s12);
end

lines{end + 1} = 'Receive-port clustering (S21 uses port2, S12 uses port1):'; %#ok<AGROW>
lines = [lines; format_receive_summary(bandAudit.receive_port_summary_s21, 'S21'); format_receive_summary(bandAudit.receive_port_summary_s12, 'S12')]; %#ok<AGROW>

lines{end + 1} = 'Geometry clustering:'; %#ok<AGROW>
geomFields = fieldnames(bandAudit.geometry_summary);
for idxGeom = 1:numel(geomFields)
    item = bandAudit.geometry_summary.(geomFields{idxGeom});
    lines{end + 1} = sprintf(['  %s: median switch/corrected S21 %.3f / %.3f dB, gain %.3f dB; ' ...
        'median switch/corrected S12 %.3f / %.3f dB, gain %.3f dB; ' ...
        '-|t21| %.3f dB, -|t12| %.3f dB, -|den| %.3f dB'], ... %#ok<AGROW>
        geomFields{idxGeom}, ...
        item.median_switch_s21_db, item.median_corrected_s21_db, item.median_gain_s21_db, ...
        item.median_switch_s12_db, item.median_corrected_s12_db, item.median_gain_s12_db, ...
        item.median_minus_t21_db, item.median_minus_t12_db, item.median_minus_den_db);
end
end

function lines = format_receive_summary(summaryStruct, label)
portFields = fieldnames(summaryStruct);
lines = cell(numel(portFields), 1);
for idx = 1:numel(portFields)
    item = summaryStruct.(portFields{idx});
    lines{idx} = sprintf(['  %s receive %s: switch %.3f dB, corrected %.3f dB, gain %.3f dB, ' ...
        '-|t| %.3f dB, -|den| %.3f dB, above-0 count %d'], ...
        label, portFields{idx}, item.median_switch_db, item.median_corrected_db, ...
        item.median_gain_db, item.median_minus_t_db, item.median_minus_den_db, ...
        item.total_above_zero_count);
end
end

function fig = create_reciprocal_scale_figure(config, rows)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1700, 1100]);
tlo = tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

bandNames = config.bands;
colors = lines(numel(bandNames));
pairLabels = arrayfun(@(r) sprintf('%s %s', r.pair, r.geometry), rows, 'UniformOutput', false);
uniquePairs = unique(pairLabels, 'stable');

ax1 = nexttile(tlo);
hold(ax1, 'on');
grid(ax1, 'on');
title(ax1, 'Median switch and corrected |S21| by pair');
ylabel(ax1, 'Magnitude (dB)');

ax2 = nexttile(tlo);
hold(ax2, 'on');
grid(ax2, 'on');
title(ax2, 'Median S21 gain decomposition by pair');
ylabel(ax2, 'Contribution (dB)');

ax3 = nexttile(tlo);
hold(ax3, 'on');
grid(ax3, 'on');
title(ax3, 'Median switch and corrected |S12| by pair');
ylabel(ax3, 'Magnitude (dB)');

ax4 = nexttile(tlo);
hold(ax4, 'on');
grid(ax4, 'on');
title(ax4, 'Median S12 gain decomposition by pair');
ylabel(ax4, 'Contribution (dB)');

ax5 = nexttile(tlo);
hold(ax5, 'on');
grid(ax5, 'on');
title(ax5, 'Median S21 gain by receive port');
ylabel(ax5, 'Gain (dB)');

ax6 = nexttile(tlo);
hold(ax6, 'on');
grid(ax6, 'on');
title(ax6, 'Median S21 gain by geometry');
ylabel(ax6, 'Gain (dB)');

x = 1:numel(uniquePairs);
for idxBand = 1:numel(bandNames)
    bandName = bandNames{idxBand};
    hits = rows(strcmpi({rows.band}, bandName));
    [~, loc] = ismember(arrayfun(@(r) sprintf('%s %s', r.pair, r.geometry), hits, 'UniformOutput', false), uniquePairs);
    plot(ax1, loc - 0.12 + 0.12 * idxBand, [hits.median_switch_s21_db], 'o', 'Color', colors(idxBand, :), 'LineWidth', 1.2, 'DisplayName', sprintf('%s switch', bandName));
    plot(ax1, loc - 0.12 + 0.12 * idxBand, [hits.median_corrected_s21_db], 'x', 'Color', colors(idxBand, :), 'LineWidth', 1.8, 'HandleVisibility', 'off');
    plot(ax2, loc, [hits.median_minus_t21_db], '-', 'Color', colors(idxBand, :), 'LineWidth', 1.2, 'DisplayName', sprintf('%s -|t21|', bandName));
    plot(ax2, loc, [hits.median_minus_den_db], '--', 'Color', colors(idxBand, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');

    plot(ax3, loc - 0.12 + 0.12 * idxBand, [hits.median_switch_s12_db], 'o', 'Color', colors(idxBand, :), 'LineWidth', 1.2, 'DisplayName', sprintf('%s switch', bandName));
    plot(ax3, loc - 0.12 + 0.12 * idxBand, [hits.median_corrected_s12_db], 'x', 'Color', colors(idxBand, :), 'LineWidth', 1.8, 'HandleVisibility', 'off');
    plot(ax4, loc, [hits.median_minus_t12_db], '-', 'Color', colors(idxBand, :), 'LineWidth', 1.2, 'DisplayName', sprintf('%s -|t12|', bandName));
    plot(ax4, loc, [hits.median_minus_den_db], '--', 'Color', colors(idxBand, :), 'LineWidth', 1.0, 'HandleVisibility', 'off');
end

set([ax1, ax2, ax3, ax4], 'XTick', x, 'XTickLabel', uniquePairs, 'XTickLabelRotation', 35);
yline(ax1, 0, '--', 'Color', [0.6 0.6 0.6]);
yline(ax3, 0, '--', 'Color', [0.6 0.6 0.6]);
yline(ax2, 0, '--', 'Color', [0.6 0.6 0.6]);
yline(ax4, 0, '--', 'Color', [0.6 0.6 0.6]);
legend(ax1, 'Location', 'bestoutside');
legend(ax2, 'Location', 'bestoutside');
legend(ax3, 'Location', 'bestoutside');
legend(ax4, 'Location', 'bestoutside');

recvPorts = 1:4;
for idxBand = 1:numel(bandNames)
    bandName = bandNames{idxBand};
    hits = rows(strcmpi({rows.band}, bandName));
    recvVals = nan(size(recvPorts));
    geomVals = nan(1, 2);
    for idxPort = 1:numel(recvPorts)
        phits = hits([hits.port2] == recvPorts(idxPort));
        if ~isempty(phits)
            recvVals(idxPort) = median([phits.median_gain_s21_db], 'omitnan');
        end
    end
    straightHits = hits(strcmpi({hits.geometry}, 'STRAIGHT'));
    arcHits = hits(strcmpi({hits.geometry}, 'ARC'));
    if ~isempty(straightHits)
        geomVals(1) = median([straightHits.median_gain_s21_db], 'omitnan');
    end
    if ~isempty(arcHits)
        geomVals(2) = median([arcHits.median_gain_s21_db], 'omitnan');
    end

    plot(ax5, recvPorts, recvVals, '-o', 'Color', colors(idxBand, :), 'LineWidth', 1.4, 'DisplayName', bandName);
    plot(ax6, 1:2, geomVals, '-o', 'Color', colors(idxBand, :), 'LineWidth', 1.4, 'DisplayName', bandName);
end
set(ax5, 'XTick', recvPorts, 'XTickLabel', {'P1','P2','P3','P4'});
set(ax6, 'XTick', 1:2, 'XTickLabel', {'STRAIGHT','ARC'});
yline(ax5, 0, '--', 'Color', [0.6 0.6 0.6]);
yline(ax6, 0, '--', 'Color', [0.6 0.6 0.6]);
legend(ax5, 'Location', 'best');
legend(ax6, 'Location', 'best');

title(tlo, 'Phase II downstream-only reciprocal transmission-scale audit', 'FontWeight', 'bold');
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
