function outputs = generate_figure5_3(summary)
%GENERATE_FIGURE5_3 Create dissertation Figure 5.3 for Open/Short extraction.

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
figureDir = fullfile(scriptDir, 'outputs', 'figures');
noteDir = fullfile(scriptDir, 'outputs', 'notes');
phase1MatDir = fullfile(repoDir, 'Phase I Extraction', 'outputs', 'mat');
ensure_dir(figureDir);
ensure_dir(noteDir);
cleanup_legacy_outputs(figureDir, {'Figure5_3a_*', 'Figure5_3b_*', 'Figure5_3c_*'});

loaded = load(fullfile(phase1MatDir, 'PHASE1_RESULTS_ALL.mat'), 'results');
results = loaded.results;
bandMap = build_band_map(results.bands);
z0 = results.config.z0;

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [100, 100, 1400, 980]);
t = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

style = default_plot_style();

ax1 = nexttile(t, 1);
plot_standard_panel(ax1, summary.accepted.OPEN.bands, summary.stitched.OPEN, bandMap, z0, ...
    'magnitude', 'Open', style);
title(ax1, 'Open Magnitude', 'FontName', style.fontName, 'FontSize', 18, 'FontWeight', 'bold');

ax2 = nexttile(t, 2);
plot_standard_panel(ax2, summary.accepted.OPEN.bands, summary.stitched.OPEN, bandMap, z0, ...
    'phase', 'Open', style);
title(ax2, 'Open Phase', 'FontName', style.fontName, 'FontSize', 18, 'FontWeight', 'bold');

ax3 = nexttile(t, 3);
plot_standard_panel(ax3, summary.accepted.SHORT.bands, summary.stitched.SHORT, bandMap, z0, ...
    'magnitude', 'Short', style);
title(ax3, 'Short Magnitude', 'FontName', style.fontName, 'FontSize', 18, 'FontWeight', 'bold');

ax4 = nexttile(t, 4);
plot_standard_panel(ax4, summary.accepted.SHORT.bands, summary.stitched.SHORT, bandMap, z0, ...
    'phase', 'Short', style);
title(ax4, 'Short Phase', 'FontName', style.fontName, 'FontSize', 18, 'FontWeight', 'bold');

title(t, 'Measured, Thru-Only Baseline, Final Extracted, Stitched, and Simulated Open/Short Responses', ...
    'FontName', style.fontName, 'FontSize', 22, 'FontWeight', 'bold');

lgd = legend(ax1, create_main_legend_handles(ax1, style, 'Open'), ...
    {'Measured at probe-tip plane', 'Thru-only baseline extraction', ...
    'Accepted bandwise extraction', 'Simulated de-embedded reference', ...
    'Stitched final reference'}, ...
    'Location', 'southwest', 'FontName', style.fontName, 'FontSize', 12);
lgd.Layout.Tile = 'east';

basePath = fullfile(figureDir, 'Figure5_3_OpenShort_ExtractionComparison');
save_figure_outputs(fig, basePath);

noteLines = {
    'Figure 5.3 source note'
    ''
    'Figure 5.3 replaces the placeholder metric-only figure with direct Open/Short response comparisons.'
    'Each panel overlays five quantities:'
    '  1. measured response at the probe-tip-calibrated plane (winner realization only)'
    '  2. thru-only baseline extraction obtained by de-embedding the measured response with embedding_baseline'
    '  3. final accepted bandwise extraction'
    '  4. stitched final full-span reference'
    '  5. simulated de-embedded trend reference'
    ''
    'The thru-only baseline is recomputed inside generate_figure5_3.m from PHASE1_RESULTS_ALL.mat'
    'using the saved embedding_baseline and the winner-port measured_gamma.'
    ''
    'The stitched curves are taken from FINALFIT_OPEN_STITCHED_FINAL.mat and FINALFIT_SHORT_STITCHED_FINAL.mat.'
    };
write_text(fullfile(noteDir, 'Figure5_3_Sources.txt'), strjoin(noteLines, newline));

outputs = struct();
outputs.figure_dir = figureDir;
outputs.figure5_3_base = basePath;
end

function plot_standard_panel(ax, bands, stitched, bandMap, z0, modeName, stdLabel, style)
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontName', style.fontName, 'FontSize', 15, ...
    'LineWidth', 1.0, 'Color', 'w');

stitchedVals = transform_trace(stitched.freq, stitched.gamma_fit, modeName, z0);

for iBand = 1:numel(bands)
    entry = bands(iBand);
    bandResult = bandMap(char(entry.band));
    baselineGamma = compute_thru_only_baseline(entry, bandResult.embedding_baseline.S);

    measuredVals = transform_trace(entry.freq, entry.measured_gamma, modeName, z0, stitched.freq, stitched.gamma_fit);
    baselineVals = transform_trace(entry.freq, baselineGamma, modeName, z0, stitched.freq, stitched.gamma_fit);
    finalVals = transform_trace(entry.freq, entry.gamma, modeName, z0, stitched.freq, stitched.gamma_fit);
    simVals = transform_trace(entry.freq, entry.sim_reference, modeName, z0, stitched.freq, stitched.gamma_fit);

    plot(ax, entry.freq / 1e9, measuredVals, ...
        'Color', style.measuredColor, 'LineWidth', 1.2, 'LineStyle', '--');
    plot(ax, entry.freq / 1e9, baselineVals, ...
        'Color', style.baselineColor, 'LineWidth', 1.4, 'LineStyle', ':');
    plot(ax, entry.freq / 1e9, finalVals, ...
        'Color', style.finalColor.(upper(stdLabel)), 'LineWidth', 2.0, 'LineStyle', '-');
    plot(ax, entry.freq / 1e9, simVals, ...
        'Color', style.simColor, 'LineWidth', 1.6, 'LineStyle', '-.');
end

plot(ax, stitched.freq / 1e9, stitchedVals, ...
    'Color', style.stitchedColor, 'LineWidth', 2.8, 'LineStyle', '-');

xlim(ax, [0, 170]);
xlabel(ax, 'Frequency (GHz)', 'FontName', style.fontName, 'FontSize', 17);
switch modeName
    case 'magnitude'
        ylabel(ax, '|\Gamma| (dB)', 'FontName', style.fontName, 'FontSize', 17);
    case 'phase'
        ylabel(ax, 'Phase (deg)', 'FontName', style.fontName, 'FontSize', 17);
        ylim(ax, [-180, 180]);
end
end

function baselineGamma = compute_thru_only_baseline(entry, embeddingBaselineS)
portNum = local_port_from_label(entry.winner_port);
gammaMeasured = entry.measured_gamma(:);
numFreq = numel(gammaMeasured);
baselineGamma = nan(numFreq, 1);
for idx = 1:numFreq
    s11 = embeddingBaselineS(1, 1, idx);
    s12 = embeddingBaselineS(1, 2, idx);
    s21 = embeddingBaselineS(2, 1, idx);
    s22 = embeddingBaselineS(2, 2, idx);
    if portNum == 1
        deltaGamma = gammaMeasured(idx) - s11;
        denominator = (s12 * s21) + (s22 * deltaGamma);
    else
        deltaGamma = gammaMeasured(idx) - s22;
        denominator = (s12 * s21) + (s11 * deltaGamma);
    end
    if abs(denominator) >= eps
        baselineGamma(idx) = deltaGamma / denominator;
    end
end
end

function portNum = local_port_from_label(portLabel)
if strcmpi(char(portLabel), 'P1') || strcmpi(char(portLabel), 'P3')
    portNum = 1;
else
    portNum = 2;
end
end

function values = transform_trace(freq, gamma, modeName, z0, refFreq, refGamma)
gamma = gamma(:);
switch modeName
    case 'magnitude'
        values = 20 * log10(max(abs(gamma), 1e-12));
    case 'phase'
        if nargin < 5
            refFreq = [];
            refGamma = [];
        end
        values = phase_trace_degrees(freq, gamma, refFreq, refGamma);
    case 'realz'
        z = z0 .* (1 + gamma) ./ max(1 - gamma, 1e-12);
        values = real(z);
    case 'imagz'
        z = z0 .* (1 + gamma) ./ max(1 - gamma, 1e-12);
        values = imag(z);
    otherwise
        error('Unsupported transform mode: %s', modeName);
end
end

function values = phase_trace_degrees(freq, gamma, refFreq, refGamma)
freq = freq(:);
gamma = gamma(:);
values = rad2deg(unwrap(angle(gamma)));

if nargin < 4 || isempty(refFreq) || isempty(refGamma)
    return;
end

refFreq = refFreq(:);
refVals = rad2deg(unwrap(angle(refGamma(:))));
refInterp = interp1(refFreq, refVals, freq, 'linear', NaN);
validMask = isfinite(values) & isfinite(refInterp);
if any(validMask)
    delta = median(refInterp(validMask) - values(validMask), 'omitnan');
    values = values + 360 * round(delta / 360);
end
end

function handles = create_main_legend_handles(ax, style, stdLabel)
stdColor = style.finalColor.(upper(stdLabel));
handles = [
    plot(ax, nan, nan, '--', 'Color', style.measuredColor, 'LineWidth', 1.2);
    plot(ax, nan, nan, ':',  'Color', style.baselineColor, 'LineWidth', 1.4);
    plot(ax, nan, nan, '-',  'Color', stdColor, 'LineWidth', 2.0);
    plot(ax, nan, nan, '-.', 'Color', style.simColor, 'LineWidth', 1.6);
    plot(ax, nan, nan, '-',  'Color', style.stitchedColor, 'LineWidth', 2.8)
    ];
end

function style = default_plot_style()
style = struct();
style.fontName = 'Times New Roman';
style.measuredColor = [0.45, 0.45, 0.45];
style.baselineColor = [0.49, 0.18, 0.56];
style.simColor = [0.20, 0.75, 0.20];
style.stitchedColor = [0.05, 0.05, 0.05];
style.finalColor = struct( ...
    'OPEN', [0.0000, 0.4470, 0.7410], ...
    'SHORT', [0.8500, 0.3250, 0.0980], ...
    'LOAD', [0.9290, 0.6940, 0.1250]);
end

function bandMap = build_band_map(bands)
keys = cell(numel(bands), 1);
vals = cell(numel(bands), 1);
for k = 1:numel(bands)
    keys{k} = char(bands{k}.band);
    vals{k} = bands{k};
end
bandMap = containers.Map(keys, vals);
end

function ensure_dir(pathStr)
if ~exist(pathStr, 'dir')
    mkdir(pathStr);
end
end

function cleanup_legacy_outputs(figureDir, patterns)
for iPat = 1:numel(patterns)
    files = dir(fullfile(figureDir, patterns{iPat}));
    for iFile = 1:numel(files)
        delete(fullfile(files(iFile).folder, files(iFile).name));
    end
end
end

function write_text(pathStr, textBody)
fid = fopen(pathStr, 'w');
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', textBody);
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
