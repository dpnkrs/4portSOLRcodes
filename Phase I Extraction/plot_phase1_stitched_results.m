function plot_phase1_stitched_results(config, results)
%PLOT_PHASE1_STITCHED_RESULTS Plot stitched full-band SOL and thru responses.

figureDir = fullfile(fileparts(config.touchstone_output_dir), 'figures');
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

solMap = containers.Map();
thruMap = containers.Map();

for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    outputFields = fieldnames(bandResult.outputs);
    for idxField = 1:numel(outputFields)
        fieldName = outputFields{idxField};
        output = bandResult.outputs.(fieldName);
        baseKey = regexprep(fieldName, '_(0|67|110)\d*[-_]\d+$', '');
        baseKey = regexprep(baseKey, '_\d+_\d+$', '');

        if isfield(output, 'gamma')
            entry = struct('freq', output.freq, 'gamma', output.gamma, 'measured', output.measured_gamma);
            solMap = append_map_entry(solMap, baseKey, entry);
        elseif isfield(output, 'S')
            entry = struct('freq', output.freq, 's21', squeeze(output.S(2, 1, :)), 'measured', squeeze(output.measured_S(2, 1, :)));
            thruMap = append_map_entry(thruMap, baseKey, entry);
        end
    end
end

plot_stitched_sol(solMap, figureDir);
plot_stitched_thru(thruMap, figureDir);
end

function mapObj = append_map_entry(mapObj, key, entry)
if isKey(mapObj, key)
    entries = mapObj(key);
    entries{end + 1} = entry;
    mapObj(key) = entries;
else
    mapObj(key) = {entry};
end
end

function plot_stitched_sol(solMap, figureDir)
if solMap.Count == 0
    return;
end

keysList = sort(solMap.keys);

fig = figure('Visible', 'off', 'Color', 'w');
hold on;
for idxKey = 1:numel(keysList)
    stitched = stitch_gamma_entries(solMap(keysList{idxKey}));
    plot(stitched.freq / 1e9, 20 * log10(abs(stitched.measured)), '--', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('%s Measured', prettify_name(keysList{idxKey})));
    plot(stitched.freq / 1e9, 20 * log10(abs(stitched.gamma)), '-', 'LineWidth', 1.6, ...
        'DisplayName', sprintf('%s De-embedded', prettify_name(keysList{idxKey})));
end
title('Phase I Stitched SOL Magnitude - 10 MHz to 170 GHz');
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
grid on;
legend('Location', 'best');
apply_axes_style(gca);
exportgraphics(fig, fullfile(figureDir, 'PhaseI_Stitched_SOL_Magnitude.jpg'), 'Resolution', 300, 'BackgroundColor', 'white');
    save_visible_fig(fig, fullfile(figureDir, 'PhaseI_Stitched_SOL_Magnitude.fig'));
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
hold on;
for idxKey = 1:numel(keysList)
    stitched = stitch_gamma_entries(solMap(keysList{idxKey}));
    plot(stitched.freq / 1e9, rad2deg(unwrap(angle(stitched.measured))), '--', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('%s Measured', prettify_name(keysList{idxKey})));
    plot(stitched.freq / 1e9, rad2deg(unwrap(angle(stitched.gamma))), '-', 'LineWidth', 1.6, ...
        'DisplayName', sprintf('%s De-embedded', prettify_name(keysList{idxKey})));
end
title('Phase I Stitched SOL Phase - 10 MHz to 170 GHz');
xlabel('Frequency (GHz)');
ylabel('Phase (deg)');
grid on;
legend('Location', 'best');
apply_axes_style(gca);
exportgraphics(fig, fullfile(figureDir, 'PhaseI_Stitched_SOL_Phase.jpg'), 'Resolution', 300, 'BackgroundColor', 'white');
    save_visible_fig(fig, fullfile(figureDir, 'PhaseI_Stitched_SOL_Phase.fig'));
close(fig);
end

function plot_stitched_thru(thruMap, figureDir)
if thruMap.Count == 0
    return;
end

keysList = sort(thruMap.keys);

fig = figure('Visible', 'off', 'Color', 'w');
hold on;
for idxKey = 1:numel(keysList)
    stitched = stitch_s21_entries(thruMap(keysList{idxKey}));
    plot(stitched.freq / 1e9, 20 * log10(abs(stitched.measured)), '--', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('%s Measured', prettify_name(keysList{idxKey})));
    plot(stitched.freq / 1e9, 20 * log10(abs(stitched.s21)), '-', 'LineWidth', 1.6, ...
        'DisplayName', sprintf('%s De-embedded', prettify_name(keysList{idxKey})));
end
title('Phase I Stitched Thru S21 Magnitude - 10 MHz to 170 GHz');
xlabel('Frequency (GHz)');
ylabel('|S21| (dB)');
grid on;
legend('Location', 'best');
apply_axes_style(gca);
exportgraphics(fig, fullfile(figureDir, 'PhaseI_Stitched_THRU_S21_Magnitude.jpg'), 'Resolution', 300, 'BackgroundColor', 'white');
    save_visible_fig(fig, fullfile(figureDir, 'PhaseI_Stitched_THRU_S21_Magnitude.fig'));
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
hold on;
for idxKey = 1:numel(keysList)
    stitched = stitch_s21_entries(thruMap(keysList{idxKey}));
    plot(stitched.freq / 1e9, rad2deg(unwrap(angle(stitched.measured))), '--', 'LineWidth', 1.0, ...
        'DisplayName', sprintf('%s Measured', prettify_name(keysList{idxKey})));
    plot(stitched.freq / 1e9, rad2deg(unwrap(angle(stitched.s21))), '-', 'LineWidth', 1.6, ...
        'DisplayName', sprintf('%s De-embedded', prettify_name(keysList{idxKey})));
end
title('Phase I Stitched Thru S21 Phase - 10 MHz to 170 GHz');
xlabel('Frequency (GHz)');
ylabel('Phase (deg)');
grid on;
legend('Location', 'best');
apply_axes_style(gca);
exportgraphics(fig, fullfile(figureDir, 'PhaseI_Stitched_THRU_S21_Phase.jpg'), 'Resolution', 300, 'BackgroundColor', 'white');
    save_visible_fig(fig, fullfile(figureDir, 'PhaseI_Stitched_THRU_S21_Phase.fig'));
close(fig);
end

function stitched = stitch_gamma_entries(entries)
[freq, gamma, measured] = stitch_common_with_overlap(entries, 'gamma');
stitched = struct('freq', freq, 'gamma', gamma, 'measured', measured);
end

function stitched = stitch_s21_entries(entries)
[freq, s21, measured] = stitch_common_with_overlap(entries, 's21');
stitched = struct('freq', freq, 's21', s21, 'measured', measured);
end

function [freq, mainData, measuredData] = stitch_common_with_overlap(entries, fieldName)
overlaps = [63e9, 70e9; 112e9, 118e9];

startFreq = zeros(numel(entries), 1);
for idx = 1:numel(entries)
    startFreq(idx) = entries{idx}.freq(1);
end
[~, order] = sort(startFreq);
entries = entries(order);

freq = entries{1}.freq(:);
mainData = entries{1}.(fieldName)(:);
measuredData = entries{1}.measured(:);

for idx = 2:numel(entries)
    upperFreq = entries{idx}.freq(:);
    upperMain = entries{idx}.(fieldName)(:);

    overlap = find_overlap_window(freq(end), upperFreq(1), overlaps);
    [freq, mainData] = merge_with_overlap(freq, mainData, upperFreq, upperMain, overlap);
end

% Rebuild measured path in a second pass to avoid depending on merged freq state.
freq = entries{1}.freq(:);
measuredData = entries{1}.measured(:);
for idx = 2:numel(entries)
    upperFreq = entries{idx}.freq(:);
    upperMeasured = entries{idx}.measured(:);
    overlap = find_overlap_window(freq(end), upperFreq(1), overlaps);
    [freq, measuredData] = merge_with_overlap(freq, measuredData, upperFreq, upperMeasured, overlap);
end
mainData = mainData(:);
measuredData = measuredData(:);
end

function overlap = find_overlap_window(lowerEndFreq, upperStartFreq, overlaps)
overlap = [];
for idx = 1:size(overlaps, 1)
    if lowerEndFreq >= overlaps(idx, 1) && upperStartFreq <= overlaps(idx, 2)
        overlap = overlaps(idx, :);
        return;
    end
end
end

function [freqMerged, dataMerged] = merge_with_overlap(lowerFreq, lowerData, upperFreq, upperData, overlap)
if isempty(overlap)
    freqMerged = [lowerFreq; upperFreq];
    dataMerged = [lowerData; upperData];
    [freqMerged, uniqueIdx] = unique(freqMerged, 'stable');
    dataMerged = dataMerged(uniqueIdx);
    return;
end

lowerOverlapMask = lowerFreq >= overlap(1) & lowerFreq <= overlap(2);
upperOverlapMask = upperFreq >= overlap(1) & upperFreq <= overlap(2);

if ~any(lowerOverlapMask) || ~any(upperOverlapMask)
    freqMerged = [lowerFreq; upperFreq];
    dataMerged = [lowerData; upperData];
    [freqMerged, uniqueIdx] = unique(freqMerged, 'stable');
    dataMerged = dataMerged(uniqueIdx);
    return;
end

commonFreq = unique([lowerFreq(lowerOverlapMask); upperFreq(upperOverlapMask)]);
lowerInterp = interp1(lowerFreq, lowerData, commonFreq, 'linear');
upperInterp = interp1(upperFreq, upperData, commonFreq, 'linear');

validMask = ~isnan(lowerInterp) & ~isnan(upperInterp);
commonFreq = commonFreq(validMask);
lowerInterp = lowerInterp(validMask);
upperInterp = upperInterp(validMask);

if isempty(commonFreq)
    freqMerged = [lowerFreq; upperFreq];
    dataMerged = [lowerData; upperData];
    [freqMerged, uniqueIdx] = unique(freqMerged, 'stable');
    dataMerged = dataMerged(uniqueIdx);
    return;
end

alignFactor = sum(lowerInterp .* conj(upperInterp)) / max(sum(abs(upperInterp).^2), eps);
upperAligned = upperData .* alignFactor;
upperInterpAligned = interp1(upperFreq, upperAligned, commonFreq, 'linear');

w = (commonFreq - commonFreq(1)) / max(commonFreq(end) - commonFreq(1), eps);
blended = (1 - w) .* lowerInterp + w .* upperInterpAligned;

lowerKeepMask = lowerFreq < commonFreq(1);
upperKeepMask = upperFreq > commonFreq(end);

freqMerged = [lowerFreq(lowerKeepMask); commonFreq; upperFreq(upperKeepMask)];
dataMerged = [lowerData(lowerKeepMask); blended; upperAligned(upperKeepMask)];
end

function label = prettify_name(rawName)
label = strrep(rawName, '_', ' ');
end

function apply_axes_style(ax)
if nargin < 1 || isempty(ax)
    ax = gca;
end
set(ax, ...
    'Color', 'w', ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'GridColor', [0.82, 0.82, 0.82], ...
    'GridAlpha', 0.9, ...
    'LineWidth', 1.0, ...
    'FontSize', 11, ...
    'Box', 'on');
if ~isempty(ax.Title); ax.Title.Color = 'k'; end
if ~isempty(ax.XLabel); ax.XLabel.Color = 'k'; end
if ~isempty(ax.YLabel); ax.YLabel.Color = 'k'; end
end

function save_visible_fig(fig, filename)
prevVisible = get(fig, 'Visible');
set(fig, 'Visible', 'on');
drawnow;
savefig(fig, filename);
set(fig, 'Visible', prevVisible);
end

