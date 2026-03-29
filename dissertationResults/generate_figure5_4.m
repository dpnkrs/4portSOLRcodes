function outputs = generate_figure5_4(summary)
%GENERATE_FIGURE5_4 Create dissertation Figure 5.4 for Load extraction.

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);
figureDir = fullfile(scriptDir, 'outputs', 'figures');
noteDir = fullfile(scriptDir, 'outputs', 'notes');
phase1MatDir = fullfile(repoDir, 'Phase I Extraction', 'outputs', 'mat');
ensure_dir(figureDir);
ensure_dir(noteDir);

loaded = load(fullfile(phase1MatDir, 'PHASE1_RESULTS_ALL.mat'), 'results');
results = loaded.results;
bandMap = build_band_map(results.bands);
z0 = results.config.z0;
loadBands = summary.accepted.LOAD.bands;
stitched = summary.stitched.LOAD;

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Position', [120, 120, 1450, 980]);
t = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
style = default_plot_style();

ax1 = nexttile(t, 1);
plot_load_panel(ax1, loadBands, stitched, bandMap, z0, 'magnitude', style);
hold(ax1, 'on');
plot(ax1, [0, 170], [-15, -15], 'k--', 'LineWidth', 1.0);
title(ax1, 'Load Magnitude', 'FontName', style.fontName, 'FontSize', 18, 'FontWeight', 'bold');

ax2 = nexttile(t, 2);
plot_load_panel(ax2, loadBands, stitched, bandMap, z0, 'phase', style);
title(ax2, 'Load Phase', 'FontName', style.fontName, 'FontSize', 18, 'FontWeight', 'bold');

ax3 = nexttile(t, 3);
plot_load_impedance_panel(ax3, loadBands, stitched, z0, style);
title(ax3, 'Extracted Load Impedance', 'FontName', style.fontName, 'FontSize', 18, 'FontWeight', 'bold');

ax4 = nexttile(t, 4);
plot_load_uncertainty_panel(ax4, loadBands, bandMap, z0, style);
title(ax4, 'Load Repeatability and Held-Out Residual', 'FontName', style.fontName, 'FontSize', 18, 'FontWeight', 'bold');

title(t, 'Measured, Thru-Only Baseline, Final Extracted, Stitched, and Simulated Load Responses', ...
    'FontName', style.fontName, 'FontSize', 22, 'FontWeight', 'bold');

lgd = legend(ax1, create_main_legend_handles(ax1, style), ...
    {'Measured at probe-tip plane', 'Thru-only baseline extraction', ...
    'Accepted bandwise extraction', 'Simulated de-embedded reference', ...
    'Stitched final reference', '-15 dB guide'}, ...
    'Location', 'southwest', 'FontName', style.fontName, 'FontSize', 12);
lgd.Layout.Tile = 'east';

basePath = fullfile(figureDir, 'Figure5_4_Load_ExtractionValidation');
save_figure_outputs(fig, basePath);

noteLines = {
    'Figure 5.4 source note'
    ''
    'Figure 5.4 focuses on the Load as the limiting extracted one-port standard.'
    'Panels (a) and (b) overlay the measured winner realization, the thru-only baseline extraction,'
    'the final accepted bandwise extraction, the stitched final reference, and the simulated de-embedded reference.'
    ''
    'Panel (c) converts the accepted bandwise and stitched Load responses to impedance using'
    'Z = Z0 * (1 + Gamma) / (1 - Gamma), with the simulated de-embedded reference shown as a trend overlay.'
    ''
    'Panel (d) plots two practical uncertainty indicators for the Load:'
    '  1. |Delta Gamma| between accepted and rejected redundant extractions'
    '  2. held-out Load residual magnitude under the square one-port solve using the accepted winners'
    ''
    'This figure is intended to support the claim that the Load, not matrix conditioning, sets the practical limit.'
    };
write_text(fullfile(noteDir, 'Figure5_4_Sources.txt'), strjoin(noteLines, newline));

outputs = struct();
outputs.figure_dir = figureDir;
outputs.figure5_4_base = basePath;
end

function plot_load_panel(ax, bands, stitched, bandMap, z0, modeName, style)
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
        'Color', style.finalColor.LOAD, 'LineWidth', 2.0, 'LineStyle', '-');
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
        ylim(ax, [-30, 0]);
    case 'phase'
        ylabel(ax, 'Phase (deg)', 'FontName', style.fontName, 'FontSize', 17);
        ylim(ax, [-180, 180]);
end
end

function plot_load_impedance_panel(ax, bands, stitched, z0, style)
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontName', style.fontName, 'FontSize', 15, ...
    'LineWidth', 1.0, 'Color', 'w');

for iBand = 1:numel(bands)
    entry = bands(iBand);
    zFinal = gamma_to_z(entry.gamma, z0);
    zSim = gamma_to_z(entry.sim_reference, z0);
    plot(ax, entry.freq / 1e9, real(zFinal), '-', ...
        'Color', style.finalColor.LOAD, 'LineWidth', 2.0);
    plot(ax, entry.freq / 1e9, imag(zFinal), '--', ...
        'Color', style.finalColor.LOAD, 'LineWidth', 2.0);
    plot(ax, entry.freq / 1e9, real(zSim), '-.', ...
        'Color', style.simColor, 'LineWidth', 1.4);
    plot(ax, entry.freq / 1e9, imag(zSim), ':', ...
        'Color', style.simColor, 'LineWidth', 1.4);
end

zFit = gamma_to_z(stitched.gamma_fit, z0);
plot(ax, stitched.freq / 1e9, real(zFit), '-', ...
    'Color', style.stitchedColor, 'LineWidth', 2.8);
plot(ax, stitched.freq / 1e9, imag(zFit), '--', ...
    'Color', style.stitchedColor, 'LineWidth', 2.8);

xlim(ax, [0, 170]);
xlabel(ax, 'Frequency (GHz)', 'FontName', style.fontName, 'FontSize', 17);
ylabel(ax, 'Impedance (ohm)', 'FontName', style.fontName, 'FontSize', 17);
legend(ax, create_impedance_legend_handles(ax, style), ...
    {'Re\{Z_{load}\} accepted bandwise', 'Im\{Z_{load}\} accepted bandwise', ...
    'Re\{Z_{load}\} simulated', 'Im\{Z_{load}\} simulated', ...
    'Re\{Z_{load}\} stitched final', 'Im\{Z_{load}\} stitched final'}, ...
    'Location', 'best', 'FontName', style.fontName, 'FontSize', 11);
end

function plot_load_uncertainty_panel(ax, bands, bandMap, z0, style)
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'on');
set(ax, 'FontName', style.fontName, 'FontSize', 15, ...
    'LineWidth', 1.0, 'Color', 'w');

for iBand = 1:numel(bands)
    entry = bands(iBand);
    bandResult = bandMap(char(entry.band));
    [repeatability, heldoutResidual] = compute_load_validation_traces(bands(iBand), bandResult, z0);
    plot(ax, entry.freq / 1e9, repeatability, ...
        '-', 'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 1.8);
    plot(ax, entry.freq / 1e9, heldoutResidual, ...
        '--', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.8);
end

xlim(ax, [0, 170]);
xlabel(ax, 'Frequency (GHz)', 'FontName', style.fontName, 'FontSize', 17);
ylabel(ax, 'Magnitude', 'FontName', style.fontName, 'FontSize', 17);
legend(ax, create_uncertainty_legend_handles(ax), ...
    {'|\Delta\Gamma| winner vs rejected', 'Held-out |r_L|'}, ...
    'Location', 'best', 'FontName', style.fontName, 'FontSize', 12);
end

function [repeatability, heldoutResidual] = compute_load_validation_traces(loadEntry, bandResult, z0)
openEntry = bandResult.final_selection.OPEN;
shortEntry = bandResult.final_selection.SHORT;
loadSel = bandResult.final_selection.LOAD;

openOutput = bandResult.outputs.(matlab.lang.makeValidName(openEntry.winner_output_key));
shortOutput = bandResult.outputs.(matlab.lang.makeValidName(shortEntry.winner_output_key));
loadOutput = bandResult.outputs.(matlab.lang.makeValidName(loadSel.winner_output_key));

openGamma = openOutput.gamma(:);
shortGamma = shortOutput.gamma(:);
loadGamma = loadOutput.gamma(:);

openMeas = openOutput.measured_gamma(:);
shortMeas = shortOutput.measured_gamma(:);
loadMeas = loadOutput.measured_gamma(:);

if isfield(loadOutput, 'gamma') && isfield(loadEntry, 'gamma_loser') && any(isfinite(loadEntry.gamma_loser))
    repeatability = abs(loadEntry.gamma(:) - loadEntry.gamma_loser(:));
else
    repeatability = nan(size(loadEntry.freq));
end

heldoutResidual = nan(size(loadEntry.freq));
rejectedGamma = loadEntry.gamma_loser(:);
rejectedMeas = loadEntry.measured_gamma_loser(:);

for idx = 1:numel(loadEntry.freq)
    A = [1, openGamma(idx) * openMeas(idx), -openGamma(idx); ...
         1, shortGamma(idx) * shortMeas(idx), -shortGamma(idx); ...
         1, loadGamma(idx) * loadMeas(idx), -loadGamma(idx)];
    b = [openMeas(idx); shortMeas(idx); loadMeas(idx)];
    if any(~isfinite(A(:))) || any(~isfinite(b(:))) || ...
            ~isfinite(rejectedGamma(idx)) || ~isfinite(rejectedMeas(idx))
        continue;
    end
    x = A \ b;
    e00 = x(1);
    e11 = x(2);
    de = x(3);
    gHat = (e00 - de * rejectedGamma(idx)) / (1 - e11 * rejectedGamma(idx));
    heldoutResidual(idx) = abs(rejectedMeas(idx) - gHat);
end
end

function z = gamma_to_z(gamma, z0)
gamma = gamma(:);
den = 1 - gamma;
smallMask = abs(den) < 1e-12;
den(smallMask) = 1e-12 .* exp(1j * angle(den(smallMask) + eps));
z = z0 .* (1 + gamma) ./ den;
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

function handles = create_main_legend_handles(ax, style)
handles = [
    plot(ax, nan, nan, '--', 'Color', style.measuredColor, 'LineWidth', 1.2);
    plot(ax, nan, nan, ':',  'Color', style.baselineColor, 'LineWidth', 1.4);
    plot(ax, nan, nan, '-',  'Color', style.finalColor.LOAD, 'LineWidth', 2.0);
    plot(ax, nan, nan, '-.', 'Color', style.simColor, 'LineWidth', 1.6);
    plot(ax, nan, nan, '-',  'Color', style.stitchedColor, 'LineWidth', 2.8);
    plot(ax, nan, nan, '--', 'Color', [0 0 0], 'LineWidth', 1.0)
    ];
end

function handles = create_impedance_legend_handles(ax, style)
handles = [
    plot(ax, nan, nan, '-',  'Color', style.finalColor.LOAD, 'LineWidth', 2.0);
    plot(ax, nan, nan, '--', 'Color', style.finalColor.LOAD, 'LineWidth', 2.0);
    plot(ax, nan, nan, '-.', 'Color', style.simColor, 'LineWidth', 1.4);
    plot(ax, nan, nan, ':',  'Color', style.simColor, 'LineWidth', 1.4);
    plot(ax, nan, nan, '-',  'Color', style.stitchedColor, 'LineWidth', 2.8);
    plot(ax, nan, nan, '--', 'Color', style.stitchedColor, 'LineWidth', 2.8)
    ];
end

function handles = create_uncertainty_legend_handles(ax)
handles = [
    plot(ax, nan, nan, '-',  'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 1.8);
    plot(ax, nan, nan, '--', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 1.8)
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
