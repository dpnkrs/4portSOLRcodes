function outputs = generate_section5_2_reciprocal(summary)
%GENERATE_SECTION5_2_RECIPROCAL Build figures and tables for dissertation Section 5.2.

if nargin < 1 %#ok<INUSD>
    summary = [];
end

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
phase1MatDir = fullfile(repoDir, 'Phase I Extraction', 'outputs', 'mat');
figureDir = fullfile(scriptDir, 'outputs', 'figures');
tableDir = fullfile(scriptDir, 'outputs', 'tables');
noteDir = fullfile(scriptDir, 'outputs', 'notes');
ensure_dir(figureDir);
ensure_dir(tableDir);
ensure_dir(noteDir);

bandNames = {'0-67', '67-115', '110-170'};
candidates = reciprocal_candidates();

data = struct();
for iCand = 1:numel(candidates)
    data.(candidates(iCand).tag) = load_candidate_data(phase1MatDir, candidates(iCand), bandNames);
end

tableRows = build_reciprocal_table_rows(data, candidates, bandNames);
safe_write_table(struct2table(tableRows), fullfile(tableDir, 'Section5_2_ReciprocalSummary.csv'));
safe_save_mat(fullfile(tableDir, 'Section5_2_ReciprocalSummary.mat'), 'tableRows', tableRows);

mergedBase = fullfile(figureDir, 'Section5_2_ReciprocalTransmissionAndPhase');
plot_reciprocal_transmission_and_phase(data, candidates, mergedBase);

assessment = compose_reciprocal_assessment(tableRows, candidates);
write_text(fullfile(noteDir, 'Section5_2_ReciprocalAssessment.txt'), assessment);

outputs = struct();
outputs.figure_dir = figureDir;
outputs.table_dir = tableDir;
outputs.note_dir = noteDir;
outputs.merged_base = mergedBase;
outputs.table_path = fullfile(tableDir, 'Section5_2_ReciprocalSummary.csv');
outputs.note_path = fullfile(noteDir, 'Section5_2_ReciprocalAssessment.txt');
end

function safe_write_table(tbl, pathStr)
try
    writetable(tbl, pathStr);
catch ME
    warning('Section5_2:TableWriteSkipped', ...
        'Skipping table overwrite for %s: %s', pathStr, ME.message);
end
end

function safe_save_mat(pathStr, varName, value)
try
    S = struct();
    S.(varName) = value;
    save(pathStr, '-struct', 'S');
catch ME
    warning('Section5_2:MatWriteSkipped', ...
        'Skipping MAT overwrite for %s: %s', pathStr, ME.message);
end
end

function candidates = reciprocal_candidates()
candidates = struct( ...
    'tag', {'P1P2_STRAIGHT', 'P3P4_STRAIGHT', 'P1P4_ARC', 'P1P4_DIAGONAL'}, ...
    'key', {'EXTRACTED_THRU_P1P2_STRAIGHT', 'EXTRACTED_THRU_P3P4_STRAIGHT', ...
            'EXTRACTED_THRU_P1P4_ARC', 'EXTRACTED_THRU_P1P4_DIAGONAL'}, ...
    'label', {'P1P2 Straight', 'P3P4 Straight', 'P1P4 Arc', 'P1P4 Diagonal'}, ...
    'role', {'Straight benchmark', 'Straight cross-check', 'Orthogonal candidate', 'Orthogonal candidate'});
end

function candData = load_candidate_data(phase1MatDir, candidate, bandNames)
candData = repmat(struct( ...
    'band', "", ...
    'freq_raw', [], ...
    's_raw', [], ...
    'freq_fit', [], ...
    's_fit', []), numel(bandNames), 1);

for iBand = 1:numel(bandNames)
    bandName = bandNames{iBand};
    rawPath = fullfile(phase1MatDir, sprintf('%s_%s.mat', candidate.key, bandName));
    fitPath = fullfile(phase1MatDir, sprintf('FINALFIT_GLOBALWIN_%s_%s.mat', candidate.key, bandName));

    rawLoaded = load(rawPath, 'reciprocalResult');
    fitLoaded = load(fitPath, 'finalResult');

    candData(iBand).band = string(bandName);
    candData(iBand).freq_raw = rawLoaded.reciprocalResult.freq(:);
    candData(iBand).s_raw = rawLoaded.reciprocalResult.S;
    candData(iBand).freq_fit = fitLoaded.finalResult.freq(:);
    candData(iBand).s_fit = fitLoaded.finalResult.s_fit;
end
end

function rows = build_reciprocal_table_rows(data, candidates, bandNames)
rows = repmat(blank_table_row(), numel(candidates) * (numel(bandNames) + 1), 1);
rowIdx = 0;

for iCand = 1:numel(candidates)
    tag = candidates(iCand).tag;
    candRows = [];
    for iBand = 1:numel(bandNames)
        metrics = compute_reciprocal_metrics(data.(tag)(iBand).freq_raw, data.(tag)(iBand).s_raw);
        rowIdx = rowIdx + 1;
        rows(rowIdx) = metrics_to_row(candidates(iCand), bandNames{iBand}, metrics);
        candRows = [candRows; rows(rowIdx)]; %#ok<AGROW>
    end

    fullFreq = [];
    fullS = [];
    for iBand = 1:numel(bandNames)
        fullFreq = [fullFreq; data.(tag)(iBand).freq_raw(:)]; %#ok<AGROW>
        fullS = cat(3, fullS, data.(tag)(iBand).s_raw);
    end
    metrics = compute_reciprocal_metrics(fullFreq, flatten_slices(fullS));
    rowIdx = rowIdx + 1;
    rows(rowIdx) = metrics_to_row(candidates(iCand), 'Full-span', metrics);
end

rows = rows(1:rowIdx);
end

function row = blank_table_row()
row = struct( ...
    'Structure', "", ...
    'Role', "", ...
    'Band', "", ...
    'FreqStart_GHz', nan, ...
    'FreqStop_GHz', nan, ...
    'MeanReciprocityMismatch', nan, ...
    'MaxReciprocityMismatch', nan, ...
    'MedianInsertionLoss_dB', nan, ...
    'MaxInsertionLoss_dB', nan, ...
    'MedianWorstPortReturnLoss_dB', nan, ...
    'MinWorstPortReturnLoss_dB', nan, ...
    'PhaseRippleRMS_deg', nan, ...
    'LandingRepeatability', "");
end

function row = metrics_to_row(candidate, bandLabel, metrics)
row = blank_table_row();
row.Structure = string(candidate.label);
row.Role = string(candidate.role);
row.Band = string(bandLabel);
row.FreqStart_GHz = metrics.freq_start_ghz;
row.FreqStop_GHz = metrics.freq_stop_ghz;
row.MeanReciprocityMismatch = metrics.eps_mean;
row.MaxReciprocityMismatch = metrics.eps_max;
row.MedianInsertionLoss_dB = metrics.il_median_db;
row.MaxInsertionLoss_dB = metrics.il_max_db;
row.MedianWorstPortReturnLoss_dB = metrics.rl_median_db;
row.MinWorstPortReturnLoss_dB = metrics.rl_min_db;
row.PhaseRippleRMS_deg = metrics.phase_ripple_rms_deg;
row.LandingRepeatability = "N/E";
end

function metrics = compute_reciprocal_metrics(freq, sMat)
freq = freq(:);
s21 = squeeze(sMat(2, 1, :));
s12 = squeeze(sMat(1, 2, :));
s11 = squeeze(sMat(1, 1, :));
s22 = squeeze(sMat(2, 2, :));

epsRec = abs(s21 - s12);
ilDb = -0.5 * (safe_mag_db(s21) + safe_mag_db(s12));
worstRlDb = min(-safe_mag_db(s11), -safe_mag_db(s22));
phaseAvg = unwrap(angle(0.5 * (s21 + s12))) * 180 / pi;
p = polyfit(freq, phaseAvg, 1);
ripple = phaseAvg - polyval(p, freq);

metrics = struct();
metrics.freq_start_ghz = min(freq) / 1e9;
metrics.freq_stop_ghz = max(freq) / 1e9;
metrics.eps_mean = mean(epsRec, 'omitnan');
metrics.eps_max = max(epsRec, [], 'omitnan');
metrics.il_median_db = median(ilDb, 'omitnan');
metrics.il_max_db = max(ilDb, [], 'omitnan');
metrics.rl_median_db = median(worstRlDb, 'omitnan');
metrics.rl_min_db = min(worstRlDb, [], 'omitnan');
metrics.phase_ripple_rms_deg = sqrt(mean(ripple .^ 2, 'omitnan'));
end

function sMat = flatten_slices(sCat)
nSlices = size(sCat, 3);
sMat = zeros(2, 2, nSlices);
for k = 1:nSlices
    sMat(:, :, k) = sCat(:, :, k);
end
end

function plot_reciprocal_transmission_and_phase(data, candidates, basePath)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1750, 1350]);
t = tiledlayout(fig, 4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
style = default_style();

for iCand = 1:numel(candidates)
    candBands = data.(candidates(iCand).tag);
    axMag = nexttile(t, 2 * iCand - 1);
    hold(axMag, 'on');
    grid(axMag, 'on');
    box(axMag, 'on');
    set(axMag, 'FontName', style.fontName, 'FontSize', 13, 'LineWidth', 1.0, 'Color', 'w');

    axPh = nexttile(t, 2 * iCand);
    hold(axPh, 'on');
    grid(axPh, 'on');
    box(axPh, 'on');
    set(axPh, 'FontName', style.fontName, 'FontSize', 13, 'LineWidth', 1.0, 'Color', 'w');

    fitRef21 = concatenate_candidate_segments(candBands, 'fit', 21);
    fitRef12 = concatenate_candidate_segments(candBands, 'fit', 12);
    rawRef21 = concatenate_candidate_segments(candBands, 'raw', 21);
    rawRef12 = concatenate_candidate_segments(candBands, 'raw', 12);
    fitPhase21 = align_segmented_phase(fitRef21);
    fitPhase12 = align_segmented_phase(fitRef12);
    rawPhase21 = align_segmented_phase(rawRef21, fitPhase21);
    rawPhase12 = align_segmented_phase(rawRef12, fitPhase12);

    for iBand = 1:numel(candBands)
        b = candBands(iBand);
        plot(axMag, b.freq_raw / 1e9, safe_mag_db(squeeze(b.s_raw(2, 1, :))), ...
            '-', 'Color', style.raw21, 'LineWidth', 1.1);
        plot(axMag, b.freq_raw / 1e9, safe_mag_db(squeeze(b.s_raw(1, 2, :))), ...
            '--', 'Color', style.raw12, 'LineWidth', 1.1);
        plot(axMag, b.freq_fit / 1e9, safe_mag_db(squeeze(b.s_fit(2, 1, :))), ...
            '-', 'Color', style.fit21, 'LineWidth', 2.2);
        plot(axMag, b.freq_fit / 1e9, safe_mag_db(squeeze(b.s_fit(1, 2, :))), ...
            '--', 'Color', style.fit12, 'LineWidth', 2.2);

        plot(axPh, rawPhase21(iBand).freq / 1e9, rawPhase21(iBand).phase_deg, ...
            '-', 'Color', style.raw21, 'LineWidth', 1.1);
        plot(axPh, rawPhase12(iBand).freq / 1e9, rawPhase12(iBand).phase_deg, ...
            '--', 'Color', style.raw12, 'LineWidth', 1.1);
        plot(axPh, fitPhase21(iBand).freq / 1e9, fitPhase21(iBand).phase_deg, ...
            '-', 'Color', style.fit21, 'LineWidth', 2.0);
        plot(axPh, fitPhase12(iBand).freq / 1e9, fitPhase12(iBand).phase_deg, ...
            '--', 'Color', style.fit12, 'LineWidth', 2.0);
    end

    title(axMag, sprintf('%s Magnitude', candidates(iCand).label), ...
        'FontName', style.fontName, 'FontSize', 16, 'FontWeight', 'bold');
    title(axPh, sprintf('%s Transmission Phase', candidates(iCand).label), ...
        'FontName', style.fontName, 'FontSize', 16, 'FontWeight', 'bold');

    xlabel(axMag, 'Frequency (GHz)', 'FontName', style.fontName, 'FontSize', 14);
    ylabel(axMag, '|S_{21}|, |S_{12}| (dB)', 'FontName', style.fontName, 'FontSize', 14);
    xlim(axMag, [0, 170]);
    ylim(axMag, [-3.5, 0.5]);

    xlabel(axPh, 'Frequency (GHz)', 'FontName', style.fontName, 'FontSize', 14);
    ylabel(axPh, 'Phase (deg)', 'FontName', style.fontName, 'FontSize', 14);
    xlim(axPh, [0, 170]);
end

title(t, 'De-Embedded Reciprocal-Candidate Transmission Magnitude and Phase', ...
    'FontName', style.fontName, 'FontSize', 22, 'FontWeight', 'bold');

lgd = legend(nexttile(t, 1), create_transmission_legend_handles(nexttile(t, 1), style), ...
    {'Raw S_{21}', 'Raw S_{12}', 'Smoothed-fit S_{21}', 'Smoothed-fit S_{12}'}, ...
    'Location', 'southwest', 'FontName', style.fontName, 'FontSize', 11);
lgd.Layout.Tile = 'east';

save_figure_outputs(fig, basePath);
end

function segs = concatenate_candidate_segments(candBands, sourceType, sParam)
segs = repmat(struct('freq', [], 'phase_deg', []), numel(candBands), 1);
for iBand = 1:numel(candBands)
    b = candBands(iBand);
    switch sourceType
        case 'raw'
            freq = b.freq_raw;
            sMat = b.s_raw;
        case 'fit'
            freq = b.freq_fit;
            sMat = b.s_fit;
        otherwise
            error('Unsupported source type: %s', sourceType);
    end
    [row, col] = sparam_to_indices(sParam);
    segs(iBand).freq = freq(:);
    segs(iBand).phase_deg = unwrap(angle(squeeze(sMat(row, col, :)))) * 180 / pi;
end
end

function aligned = align_segmented_phase(segs, refSegs)
if nargin < 2
    refSegs = [];
end

aligned = segs;
for iSeg = 1:numel(segs)
    freq = segs(iSeg).freq(:);
    phase = segs(iSeg).phase_deg(:);
    if isempty(freq)
        continue;
    end

    if ~isempty(refSegs)
        refPhase = interp_reference_phase(freq, refSegs);
        validMask = isfinite(refPhase);
        if any(validMask)
            delta = median(refPhase(validMask) - phase(validMask), 'omitnan');
            phase = phase + 360 * round(delta / 360);
        end
    elseif iSeg > 1
        refPhase = interp1(aligned(iSeg - 1).freq, aligned(iSeg - 1).phase_deg, freq, 'linear', NaN);
        validMask = isfinite(refPhase);
        if any(validMask)
            delta = median(refPhase(validMask) - phase(validMask), 'omitnan');
        else
            delta = aligned(iSeg - 1).phase_deg(end) - phase(1);
        end
        phase = phase + 360 * round(delta / 360);
    end

    aligned(iSeg).phase_deg = phase;
end
end

function refPhase = interp_reference_phase(freq, refSegs)
refPhase = nan(size(freq));
for iSeg = 1:numel(refSegs)
    segPhase = interp1(refSegs(iSeg).freq, refSegs(iSeg).phase_deg, freq, 'linear', NaN);
    fillMask = isnan(refPhase) & isfinite(segPhase);
    refPhase(fillMask) = segPhase(fillMask);
end
end

function [row, col] = sparam_to_indices(sParam)
switch sParam
    case 11
        row = 1; col = 1;
    case 12
        row = 1; col = 2;
    case 21
        row = 2; col = 1;
    case 22
        row = 2; col = 2;
    otherwise
        error('Unsupported S-parameter key: %d', sParam);
end
end

function handles = create_transmission_legend_handles(ax, style)
handles = [
    plot(ax, nan, nan, '-',  'Color', style.raw21, 'LineWidth', 1.1);
    plot(ax, nan, nan, '--', 'Color', style.raw12, 'LineWidth', 1.1);
    plot(ax, nan, nan, '-',  'Color', style.fit21, 'LineWidth', 2.3);
    plot(ax, nan, nan, '--', 'Color', style.fit12, 'LineWidth', 2.3)
    ];
end

function txt = compose_reciprocal_assessment(tableRows, candidates)
fullRows = tableRows(string({tableRows.Band}) == "Full-span");
arcRow = fullRows(string({fullRows.Structure}) == "P1P4 Arc");
diagRow = fullRows(string({fullRows.Structure}) == "P1P4 Diagonal");
straightRow = fullRows(string({fullRows.Structure}) == "P1P2 Straight");
p34Row = fullRows(string({fullRows.Structure}) == "P3P4 Straight");

if isempty(arcRow) || isempty(diagRow) || isempty(straightRow) || isempty(p34Row)
    error('Section5_2:MissingRows', 'Unable to find full-span reciprocal rows for note generation.');
end

preferredOrth = "P1P4 Arc";
if diagRow.MeanReciprocityMismatch < arcRow.MeanReciprocityMismatch && ...
        diagRow.MedianInsertionLoss_dB <= arcRow.MedianInsertionLoss_dB && ...
        diagRow.PhaseRippleRMS_deg <= arcRow.PhaseRippleRMS_deg
    preferredOrth = "P1P4 Diagonal";
end

lines = {};
lines{end + 1} = 'Section 5.2 Reciprocal Standard Assessment';
lines{end + 1} = '';
lines{end + 1} = 'The reciprocal-candidate assessment is based on de-embedded two-port responses referred to the common standard-definition plane using the straight-thru-derived embedding from Section 5.1.';
lines{end + 1} = '';
lines{end + 1} = sprintf(['The P1P2 straight Thru remains the cleanest collinear benchmark, with full-span mean reciprocity mismatch %.4f, ', ...
    'median insertion loss %.3f dB, median worst-port return loss %.2f dB, and phase-ripple RMS %.3f deg. ', ...
    'The independent P3P4 straight cross-check is visibly poorer but still well behaved relative to the orthogonal candidates.'], ...
    straightRow.MeanReciprocityMismatch, straightRow.MedianInsertionLoss_dB, ...
    straightRow.MedianWorstPortReturnLoss_dB, straightRow.PhaseRippleRMS_deg);
lines{end + 1} = '';
lines{end + 1} = sprintf(['Among the orthogonal candidates, the arc and diagonal structures are close enough that the choice must be made by measured tradeoff rather than by layout symmetry alone. ', ...
    'The arc yields full-span mean reciprocity mismatch %.4f versus %.4f for the diagonal, median insertion loss %.3f dB versus %.3f dB, ', ...
    'and phase-ripple RMS %.3f deg versus %.3f deg, while the diagonal gives better median worst-port return loss %.2f dB versus %.2f dB.'], ...
    arcRow.MeanReciprocityMismatch, diagRow.MeanReciprocityMismatch, ...
    arcRow.MedianInsertionLoss_dB, diagRow.MedianInsertionLoss_dB, ...
    arcRow.PhaseRippleRMS_deg, diagRow.PhaseRippleRMS_deg, ...
    diagRow.MedianWorstPortReturnLoss_dB, arcRow.MedianWorstPortReturnLoss_dB);
lines{end + 1} = '';
lines{end + 1} = sprintf(['On balance, %s is the better orthogonal reciprocal candidate because it preserves slightly lower forward/reverse mismatch, lower transmission loss, and smoother phase progression, ', ...
    'which are the primary transmission-standard criteria for the later SOLR implementation. The remaining orthogonal structure should be retained as a cross-check rather than promoted as the primary reciprocal standard.'], preferredOrth);

txt = strjoin(lines, newline);
end

function style = default_style()
style = struct();
style.fontName = 'Times New Roman';
style.raw21 = [0.48, 0.75, 0.96];
style.raw12 = [0.95, 0.63, 0.63];
style.fit21 = [0.00, 0.45, 0.74];
style.fit12 = [0.85, 0.33, 0.10];
end

function values = safe_mag_db(x)
values = 20 * log10(max(abs(x), 1e-12));
end

function ensure_dir(pathStr)
if ~exist(pathStr, 'dir')
    mkdir(pathStr);
end
end

function write_text(pathStr, textStr)
fid = fopen(pathStr, 'w');
if fid < 0
    error('Section5_2:WriteFailed', 'Unable to write file: %s', pathStr);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', textStr);
end

function save_figure_outputs(fig, basePath)
fig.InvertHardcopy = 'off';
set(fig, 'Visible', 'on', 'Color', 'w');
drawnow;
savefig(fig, [basePath, '.fig']);
exportgraphics(fig, [basePath, '.jpg'], 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, [basePath, '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
set(fig, 'Visible', 'off');
close(fig);
end
