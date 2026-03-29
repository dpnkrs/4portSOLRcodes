function outputs = generate_figure5_2(summary)
%GENERATE_FIGURE5_2 Create dissertation Figure 5.2 subfigures.

scriptDir = fileparts(mfilename('fullpath'));
figureDir = fullfile(scriptDir, 'outputs', 'figures');
noteDir = fullfile(scriptDir, 'outputs', 'notes');
ensure_dir(figureDir);
ensure_dir(noteDir);

stdList = {'OPEN', 'SHORT', 'LOAD'};
labelMap = struct('OPEN', 'Open', 'SHORT', 'Short', 'LOAD', 'Load');
colorMap = struct( ...
    'OPEN', [0.0000, 0.4470, 0.7410], ...
    'SHORT', [0.8500, 0.3250, 0.0980], ...
    'LOAD', [0.9290, 0.6940, 0.1250]);

figMag = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1200, 720]);
axMag = axes(figMag);
hold(axMag, 'on');
grid(axMag, 'on');
box(axMag, 'on');
set(axMag, 'FontName', 'Times New Roman', 'FontSize', 16, ...
    'LineWidth', 1.0, 'Color', 'w');

legendHandles = gobjects(0);
legendLabels = {};
for iStd = 1:numel(stdList)
    stdName = stdList{iStd};
    baseColor = colorMap.(stdName);
    acceptedBands = summary.accepted.(stdName).bands;
    acceptedHandle = gobjects(0);
    for iBand = 1:numel(acceptedBands)
        entry = acceptedBands(iBand);
        freqGHz = entry.freq / 1e9;
        magDb = safe_mag_db(entry.gamma);
        h = plot(axMag, freqGHz, magDb, ...
            'Color', tint_color(baseColor, 0.25), ...
            'LineWidth', 1.6, ...
            'LineStyle', '-');
        if iBand == 1
            acceptedHandle = h;
        end
    end

    stitched = summary.stitched.(stdName);
    stitchedHandle = plot(axMag, stitched.freq / 1e9, safe_mag_db(stitched.gamma_fit), ...
        'Color', baseColor, ...
        'LineWidth', 3.0, ...
        'LineStyle', '-');

    legendHandles(end + 1) = acceptedHandle; %#ok<AGROW>
    legendLabels{end + 1} = sprintf('%s accepted bandwise', labelMap.(stdName)); %#ok<AGROW>
    legendHandles(end + 1) = stitchedHandle; %#ok<AGROW>
    legendLabels{end + 1} = sprintf('%s stitched final', labelMap.(stdName)); %#ok<AGROW>
end

xlabel(axMag, 'Frequency (GHz)', 'FontName', 'Times New Roman', 'FontSize', 18);
ylabel(axMag, 'Magnitude (dB)', 'FontName', 'Times New Roman', 'FontSize', 18);
title(axMag, 'Extracted One-Port Standards: Accepted Bandwise and Stitched Full-Span Fits', ...
    'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
xlim(axMag, [0, 170]);
ylim(axMag, [-30, 2]);
legend(axMag, legendHandles, legendLabels, ...
    'Location', 'southwest', 'FontName', 'Times New Roman', 'FontSize', 12);

magBase = fullfile(figureDir, 'Figure5_2a_SOL_Magnitude_BandwiseAndStitched');
save_figure_outputs(figMag, magBase);

figSmith = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1000, 900]);
axSmith = axes(figSmith);
hold(axSmith, 'on');
axis(axSmith, 'equal');
axis(axSmith, [-1.1, 1.1, -1.1, 1.1]);
set(axSmith, 'Color', 'w', ...
    'FontName', 'Times New Roman', 'FontSize', 16, ...
    'LineWidth', 1.0, ...
    'XColor', 'none', 'YColor', 'none');
box(axSmith, 'off');
draw_smith_grid(axSmith);

legendHandles = gobjects(0);
legendLabels = {};
for iStd = 1:numel(stdList)
    stdName = stdList{iStd};
    baseColor = colorMap.(stdName);
    acceptedBands = summary.accepted.(stdName).bands;
    acceptedHandle = gobjects(0);
    for iBand = 1:numel(acceptedBands)
        entry = acceptedBands(iBand);
        h = plot(axSmith, real(entry.gamma), imag(entry.gamma), ...
            'Color', tint_color(baseColor, 0.25), ...
            'LineWidth', 1.6, ...
            'LineStyle', '-');
        if iBand == 1
            acceptedHandle = h;
        end
    end

    stitched = summary.stitched.(stdName);
    stitchedHandle = plot(axSmith, real(stitched.gamma_fit), imag(stitched.gamma_fit), ...
        'Color', baseColor, ...
        'LineWidth', 3.0, ...
        'LineStyle', '-');

    legendHandles(end + 1) = acceptedHandle; %#ok<AGROW>
    legendLabels{end + 1} = sprintf('%s accepted bandwise', labelMap.(stdName)); %#ok<AGROW>
    legendHandles(end + 1) = stitchedHandle; %#ok<AGROW>
    legendLabels{end + 1} = sprintf('%s stitched final', labelMap.(stdName)); %#ok<AGROW>
end

title(axSmith, 'Smith Chart of Extracted One-Port Standards', ...
    'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
legend(axSmith, legendHandles, legendLabels, ...
    'Location', 'northeastoutside', 'FontName', 'Times New Roman', 'FontSize', 12);

smithBase = fullfile(figureDir, 'Figure5_2b_SOL_Smith_BandwiseAndStitched');
save_figure_outputs(figSmith, smithBase);

noteLines = {
    'Figure 5.2 source note'
    ''
    'The stitched curves used in Figure 5.2 are the locked accepted-bandwise stitched final products:'
    '  FINALFIT_OPEN_STITCHED_FINAL.mat'
    '  FINALFIT_SHORT_STITCHED_FINAL.mat'
    '  FINALFIT_LOAD_STITCHED_FINAL.mat'
    ''
    'These stitched curves are regularized full-span representations of the accepted bandwise winners, not the alternate GLOBALWIN branch.'
    ''
    'Accepted bandwise winners used to build the stitched curves:'
    '  OPEN  : P2 (0-67), P2 (67-115), P1 (110-170)'
    '  SHORT : P1 (0-67), P2 (67-115), P2 (110-170)'
    '  LOAD  : P2 (0-67), P1 (67-115), P2 (110-170)'
    ''
    'Figure 5.2(a) plots accepted bandwise magnitude traces together with the stitched final curves.'
    'Figure 5.2(b) plots the accepted bandwise Smith-chart loci together with the stitched final curves.'
    };
write_text(fullfile(noteDir, 'Figure5_2_Sources.txt'), strjoin(noteLines, newline));

outputs = struct();
outputs.figure_dir = figureDir;
outputs.figure5_2a_base = magBase;
outputs.figure5_2b_base = smithBase;
end

function ensure_dir(pathStr)
if ~exist(pathStr, 'dir')
    mkdir(pathStr);
end
end

function values = safe_mag_db(gamma)
values = 20 * log10(max(abs(gamma), 1e-12));
end

function colorOut = tint_color(colorIn, mixFrac)
colorOut = colorIn + mixFrac * (1 - colorIn);
colorOut = min(max(colorOut, 0), 1);
end

function save_figure_outputs(fig, basePath)
fig.InvertHardcopy = 'off';
set(fig, 'Visible', 'on', 'Color', 'w');
apply_white_background(fig);
drawnow;
savefig(fig, [basePath, '.fig']);
exportgraphics(fig, [basePath, '.jpg'], 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, [basePath, '.pdf'], 'ContentType', 'vector', 'BackgroundColor', 'white');
set(fig, 'Visible', 'off');
close(fig);
end

function apply_white_background(fig)
set(fig, 'Color', 'w');

axesList = findall(fig, 'Type', 'Axes');
for k = 1:numel(axesList)
    try
        set(axesList(k), 'Color', 'w');
    catch
    end
end

panelList = findall(fig, 'Type', 'uipanel');
for k = 1:numel(panelList)
    try
        set(panelList(k), 'BackgroundColor', 'w');
    catch
    end
end

bgObjects = findall(fig, '-property', 'BackgroundColor');
for k = 1:numel(bgObjects)
    try
        set(bgObjects(k), 'BackgroundColor', 'w');
    catch
    end
end
end

function draw_smith_grid(ax)
gridColor = 0.72 * [1 1 1];
gridThin = 0.85 * [1 1 1];
theta = linspace(0, 2*pi, 1200);

plot(ax, cos(theta), sin(theta), 'Color', 0.55 * [1 1 1], 'LineWidth', 1.0);
plot(ax, [-1 1], [0 0], 'Color', 0.55 * [1 1 1], 'LineWidth', 1.0);

rValsMajor = [0.2 0.5 1 2 5];
rValsMinor = [0.1 0.3 0.7 1.5 3];
xValsMajor = [0.2 0.5 1 2 5];
xValsMinor = [0.1 0.3 0.7 1.5 3];

for r = rValsMinor
    draw_const_r(ax, r, gridThin, 0.7);
end
for r = rValsMajor
    draw_const_r(ax, r, gridColor, 0.9);
end
for x = xValsMinor
    draw_const_x(ax, x, gridThin, 0.7);
    draw_const_x(ax, -x, gridThin, 0.7);
end
for x = xValsMajor
    draw_const_x(ax, x, gridColor, 0.9);
    draw_const_x(ax, -x, gridColor, 0.9);
end

label_color = 0.15 * [1 1 1];
for r = rValsMajor
    gamma = (r - 1) / (r + 1);
    text(ax, gamma, 0.02, num2str(r), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'Color', label_color, ...
        'Rotation', 90);
end

xLabelVals = [0.2 0.5 1 2 5];
for x = xLabelVals
    gTop = (1i*x - 1) / (1i*x + 1);
    gBot = (-1i*x - 1) / (-1i*x + 1);
    text(ax, real(gTop)*1.03, imag(gTop)*1.03, sprintf('+j%g', x), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'Color', label_color);
    text(ax, real(gBot)*1.03, imag(gBot)*1.03, sprintf('-j%g', x), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontName', 'Times New Roman', 'FontSize', 12, 'Color', label_color);
end
text(ax, -1.02, 0, '0.0', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'FontName', 'Times New Roman', 'FontSize', 12, 'Color', label_color);
end

function draw_const_r(ax, r, colorVal, lineWidth)
cx = r / (1 + r);
rad = 1 / (1 + r);
t = linspace(0, 2*pi, 4000);
x = cx + rad * cos(t);
y = rad * sin(t);
mask = (x.^2 + y.^2) <= 1 + 1e-9;
plot(ax, x(mask), y(mask), 'Color', colorVal, 'LineWidth', lineWidth);
end

function draw_const_x(ax, xConst, colorVal, lineWidth)
cx = 1;
cy = 1 / xConst;
rad = abs(1 / xConst);
t = linspace(0, 2*pi, 4000);
x = cx + rad * cos(t);
y = cy + rad * sin(t);
mask = (x.^2 + y.^2) <= 1 + 1e-9 & x <= 1 + 1e-9;
plot(ax, x(mask), y(mask), 'Color', colorVal, 'LineWidth', lineWidth);
end

function write_text(pathStr, textStr)
fid = fopen(pathStr, 'w');
if fid < 0
    error('DissertationResults:WriteFailed', 'Unable to write file: %s', pathStr);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', textStr);
end
