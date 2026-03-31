function plot_phase2_band_results(config, bandResult)
%PLOT_PHASE2_BAND_RESULTS Generate diagnostic figures for one Phase II band.

figureDir = config.figure_dir;
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

local_plot_sol_residuals(config, bandResult);
local_plot_reciprocal_validation(config, bandResult);
local_plot_dut_metrics(config, bandResult);
end

function local_plot_sol_residuals(config, bandResult)
fig = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
stdNames = {'Open', 'Short', 'Load'};
freqGHz = bandResult.freq / 1e9;
colorMap = lines(12);

for idxStd = 1:numel(stdNames)
    nexttile;
    residualMatrix = bandResult.sol_residual_summary.(stdNames{idxStd});
    hold on;
    for idxTrace = 1:size(residualMatrix, 1)
        plot(freqGHz, residualMatrix(idxTrace, :), 'LineWidth', 1.1, ...
            'Color', colorMap(mod(idxTrace - 1, size(colorMap, 1)) + 1, :));
    end
    grid on;
    xlabel('Frequency (GHz)');
    ylabel('|Residual|');
    title(sprintf('%s residual vs stitched Phase I reference', stdNames{idxStd}));
end

sgtitle(sprintf('Phase II SOL residuals - %s GHz', bandResult.band));
save_figure(fig, fullfile(config.figure_dir, sprintf('PhaseII_SOL_Residuals_%s', bandResult.band)));
end

function local_plot_reciprocal_validation(config, bandResult)
fig = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
freqGHz = bandResult.freq / 1e9;

nexttile;
hold on;
colors = lines(numel(bandResult.reference_pair_results) + numel(bandResult.holdout_results));
for idx = 1:numel(bandResult.reference_pair_results)
    pairResult = bandResult.reference_pair_results(idx);
    s21 = squeeze(pairResult.reference_corrected(2, 1, :));
    plot(freqGHz, 20 * log10(max(abs(s21), 1e-12)), 'LineWidth', 1.2, ...
        'Color', colors(idx, :), 'DisplayName', sprintf('%s %s used', pairResult.pair, pairResult.geometry));
end
for idx = 1:numel(bandResult.holdout_results)
    holdout = bandResult.holdout_results(idx);
    s21 = squeeze(holdout.corrected(2, 1, :));
    plot(freqGHz, 20 * log10(max(abs(s21), 1e-12)), '--', 'LineWidth', 1.2, ...
        'Color', colors(numel(bandResult.reference_pair_results) + idx, :), ...
        'DisplayName', sprintf('%s %s holdout', holdout.pair, holdout.geometry));
end
grid on;
xlabel('Frequency (GHz)');
ylabel('|S_{21}| (dB)');
title('Corrected reciprocal-thru transmission');
legend('Location', 'bestoutside');

nexttile;
hold on;
for idx = 1:numel(bandResult.holdout_results)
    holdout = bandResult.holdout_results(idx);
    plot(freqGHz, holdout.metrics.reciprocity_mismatch, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s %s |S21-S12|', holdout.pair, holdout.geometry));
end
grid on;
xlabel('Frequency (GHz)');
ylabel('Reciprocity mismatch');
title('Hold-out DIAG reciprocity residuals');
legend('Location', 'bestoutside');

sgtitle(sprintf('Phase II reciprocal validation - %s GHz', bandResult.band));
save_figure(fig, fullfile(config.figure_dir, sprintf('PhaseII_Reciprocal_Validation_%s', bandResult.band)));
end

function local_plot_dut_metrics(config, bandResult)
dut = bandResult.dut;
fig = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
freqGHz = dut.corrected.freq / 1e9;
colorMap = lines(8);

nexttile;
hold on;
for idxPort = 2:4
    sVal = squeeze(dut.corrected.S(idxPort, 1, :));
    plot(freqGHz, 20 * log10(max(abs(sVal), 1e-12)), 'LineWidth', 1.2, ...
        'Color', colorMap(idxPort - 1, :), ...
        'DisplayName', sprintf('S%d1', idxPort));
end
grid on;
xlabel('Frequency (GHz)');
ylabel('|S_{x1}| (dB)');
title('Corrected Lange coupler transmission from port 1');
legend('Location', 'best');

nexttile;
hold on;
for idxPort = 1:4
    sVal = squeeze(dut.corrected.S(idxPort, idxPort, :));
    plot(freqGHz, -20 * log10(max(abs(sVal), 1e-12)), 'LineWidth', 1.2, ...
        'Color', colorMap(idxPort, :), ...
        'DisplayName', sprintf('RL%d', idxPort));
end
grid on;
xlabel('Frequency (GHz)');
ylabel('Return loss (dB)');
title('Corrected Lange coupler return loss');
legend('Location', 'best');

sgtitle(sprintf('Phase II DUT diagnostics - %s GHz', bandResult.band));
save_figure(fig, fullfile(config.figure_dir, sprintf('PhaseII_DUT_Diagnostics_%s', bandResult.band)));
end

function save_figure(fig, basePath)
set(fig, 'Color', 'w');
ax = findall(fig, 'type', 'axes');
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
drawnow;
exportgraphics(fig, [basePath, '.jpg'], 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, [basePath, '.pdf'], 'BackgroundColor', 'white', 'ContentType', 'vector');
savefig(fig, [basePath, '.fig']);
close(fig);
end
