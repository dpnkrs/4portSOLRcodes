function plot_phase1_band_results(config, bandResult)
%PLOT_PHASE1_BAND_RESULTS Plot de-embedded standard magnitude and phase.

figureDir = fullfile(fileparts(config.touchstone_output_dir), 'figures');
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

if isfield(bandResult, 'alpha_fit') && ~isempty(bandResult.alpha_fit)
    fig = figure('Visible', 'on');
    alphaScores = squeeze(min(min(bandResult.alpha_fit.total_scores, [], 3), [], 2));
    plot(bandResult.alpha_fit.alpha_grid_np_per_m, alphaScores, 'LineWidth', 1.8);
    hold on;
    xline(bandResult.alpha_fit.best_alpha_np_per_m, '--r', 'LineWidth', 1.2);
    title(sprintf('Phase I Alpha-Fit Objective - %s GHz', bandResult.band));
    xlabel('\alpha (Np/m)');
    ylabel('Objective');
    grid on;
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_AlphaFit_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_AlphaFit_%s.fig', bandResult.band)));
    close(fig);

    fig = figure('Visible', 'on');
    epsScores = squeeze(min(min(bandResult.alpha_fit.total_scores, [], 3), [], 1));
    plot(bandResult.alpha_fit.eps_eff_grid, epsScores, 'LineWidth', 1.8);
    hold on;
    xline(bandResult.alpha_fit.best_eps_eff, '--r', 'LineWidth', 1.2);
    title(sprintf('Phase I EpsEff-Fit Objective - %s GHz', bandResult.band));
    xlabel('\epsilon_{eff}');
    ylabel('Objective');
    grid on;
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_EpsEffFit_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_EpsEffFit_%s.fig', bandResult.band)));
    close(fig);

    fig = figure('Visible', 'on');
    alphaEpsSurface = squeeze(min(bandResult.alpha_fit.total_scores, [], 3));
    imagesc(bandResult.alpha_fit.eps_eff_grid, bandResult.alpha_fit.alpha_grid_np_per_m, alphaEpsSurface);
    axis xy;
    colorbar;
    hold on;
    plot(bandResult.alpha_fit.best_eps_eff, bandResult.alpha_fit.best_alpha_np_per_m, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    title(sprintf('Phase I Alpha/EpsEff Objective Surface - %s GHz', bandResult.band));
    xlabel('\epsilon_{eff}');
    ylabel('\alpha (Np/m)');
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_AlphaEpsEff_Surface_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_AlphaEpsEff_Surface_%s.fig', bandResult.band)));
    close(fig);

    fig = figure('Visible', 'on');
    zcScores = squeeze(min(min(bandResult.alpha_fit.total_scores, [], 2), [], 1));
    plot(bandResult.alpha_fit.zc_grid, zcScores, 'LineWidth', 1.8);
    hold on;
    xline(bandResult.alpha_fit.best_zc, '--r', 'LineWidth', 1.2);
    if isfield(bandResult.alpha_fit, 'refined')
        refinedZcScores = squeeze(min(min(bandResult.alpha_fit.refined.total_scores, [], 2), [], 1));
        plot(bandResult.alpha_fit.refined.zc_grid, refinedZcScores, '--', 'LineWidth', 1.4);
    end
    title(sprintf('Phase I Zc-Fit Objective - %s GHz', bandResult.band));
    xlabel('Zc (ohm)');
    ylabel('Objective');
    grid on;
    legendEntries = {'Coarse', 'Selected'};
    if isfield(bandResult.alpha_fit, 'refined')
        legendEntries = {'Coarse', 'Selected', 'Refined'};
    end
    legend(legendEntries, 'Location', 'best');
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_ZcFit_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_ZcFit_%s.fig', bandResult.band)));
    close(fig);
end

if isfield(bandResult, 'embedding_diagnostics') && bandResult.embedding_diagnostics.summary.has_p34
    plot_embedding_comparison(figureDir, bandResult);
end

if isfield(bandResult, 'refine_diagnostics') && ~isempty(bandResult.refine_diagnostics)
    plot_refinement_diagnostics(figureDir, bandResult);
end

outputFields = fieldnames(bandResult.outputs);

solNames = {};
solFreq = {};
solGamma = {};
solMeasuredGamma = {};
thruNames = {};
thruFreq = {};
thruS21 = {};
thruMeasuredS21 = {};

for idx = 1:numel(outputFields)
    output = bandResult.outputs.(outputFields{idx});
    if isfield(output, 'gamma')
        solNames{end + 1} = outputFields{idx}; %#ok<AGROW>
        solFreq{end + 1} = output.freq; %#ok<AGROW>
        solGamma{end + 1} = output.gamma; %#ok<AGROW>
        solMeasuredGamma{end + 1} = output.measured_gamma; %#ok<AGROW>
    elseif isfield(output, 'S')
        thruNames{end + 1} = outputFields{idx}; %#ok<AGROW>
        thruFreq{end + 1} = output.freq; %#ok<AGROW>
        thruS21{end + 1} = squeeze(output.S(2, 1, :)); %#ok<AGROW>
        thruMeasuredS21{end + 1} = squeeze(output.measured_S(2, 1, :)); %#ok<AGROW>
    end
end

if ~isempty(solNames)
    fig = figure('Visible', 'on');
    hold on;
    for idx = 1:numel(solNames)
        baseLabel = prettify_name(solNames{idx});
        hasDistinctRaw = isfield(bandResult.outputs.(solNames{idx}), 'gamma_raw') && ...
            any(abs(bandResult.outputs.(solNames{idx}).gamma_raw - solGamma{idx}) > 1e-12);
        plot(solFreq{idx} / 1e9, 20 * log10(abs(solMeasuredGamma{idx})), '--', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s Measured', baseLabel));
        if hasDistinctRaw
            plot(solFreq{idx} / 1e9, 20 * log10(abs(bandResult.outputs.(solNames{idx}).gamma_raw)), ':', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('%s Raw De-embedded', baseLabel));
        end
        plot(solFreq{idx} / 1e9, 20 * log10(abs(solGamma{idx})), '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s De-embedded', baseLabel));
    end
    title(sprintf('Phase I De-Embedded SOL Magnitude - %s GHz', bandResult.band));
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    grid on;
    legend('Location', 'best');
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_SOL_Magnitude_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_SOL_Magnitude_%s.fig', bandResult.band)));
    close(fig);

    fig = figure('Visible', 'on');
    hold on;
    for idx = 1:numel(solNames)
        baseLabel = prettify_name(solNames{idx});
        hasDistinctRaw = isfield(bandResult.outputs.(solNames{idx}), 'gamma_raw') && ...
            any(abs(bandResult.outputs.(solNames{idx}).gamma_raw - solGamma{idx}) > 1e-12);
        plot(solFreq{idx} / 1e9, rad2deg(unwrap(angle(solMeasuredGamma{idx}))), '--', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s Measured', baseLabel));
        if hasDistinctRaw
            plot(solFreq{idx} / 1e9, rad2deg(unwrap(angle(bandResult.outputs.(solNames{idx}).gamma_raw))), ':', 'LineWidth', 1.2, ...
                'DisplayName', sprintf('%s Raw De-embedded', baseLabel));
        end
        plot(solFreq{idx} / 1e9, rad2deg(unwrap(angle(solGamma{idx}))), '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s De-embedded', baseLabel));
    end
    title(sprintf('Phase I De-Embedded SOL Phase - %s GHz', bandResult.band));
    xlabel('Frequency (GHz)');
    ylabel('Phase (deg)');
    grid on;
    legend('Location', 'best');
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_SOL_Phase_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_SOL_Phase_%s.fig', bandResult.band)));
    close(fig);
end

if ~isempty(thruNames)
    fig = figure('Visible', 'on');
    hold on;
    for idx = 1:numel(thruNames)
        baseLabel = prettify_name(thruNames{idx});
        plot(thruFreq{idx} / 1e9, 20 * log10(abs(thruMeasuredS21{idx})), '--', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s Measured', baseLabel));
        plot(thruFreq{idx} / 1e9, 20 * log10(abs(thruS21{idx})), '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s De-embedded', baseLabel));
    end
    title(sprintf('Phase I De-Embedded Thru S21 Magnitude - %s GHz', bandResult.band));
    xlabel('Frequency (GHz)');
    ylabel('|S21| (dB)');
    grid on;
    legend('Location', 'best');
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_THRU_S21_Magnitude_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_THRU_S21_Magnitude_%s.fig', bandResult.band)));
    close(fig);

    fig = figure('Visible', 'on');
    hold on;
    for idx = 1:numel(thruNames)
        baseLabel = prettify_name(thruNames{idx});
        plot(thruFreq{idx} / 1e9, rad2deg(unwrap(angle(thruMeasuredS21{idx}))), '--', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s Measured', baseLabel));
        plot(thruFreq{idx} / 1e9, rad2deg(unwrap(angle(thruS21{idx}))), '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s De-embedded', baseLabel));
    end
    title(sprintf('Phase I De-Embedded Thru S21 Phase - %s GHz', bandResult.band));
    xlabel('Frequency (GHz)');
    ylabel('Phase (deg)');
    grid on;
    legend('Location', 'best');
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_THRU_S21_Phase_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_THRU_S21_Phase_%s.fig', bandResult.band)));
    close(fig);
end

nonStraightMask = cellfun(@(name) contains(name, 'ARC') || contains(name, 'DIAGONAL'), thruNames);
if any(nonStraightMask)
    fig = figure('Visible', 'on');
    hold on;
    idxList = find(nonStraightMask);
    for k = 1:numel(idxList)
        idx = idxList(k);
        baseLabel = prettify_name(thruNames{idx});
        plot(thruFreq{idx} / 1e9, 20 * log10(abs(thruMeasuredS21{idx})), '--', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s Measured', baseLabel));
        plot(thruFreq{idx} / 1e9, 20 * log10(abs(thruS21{idx})), '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s De-embedded', baseLabel));
    end
    title(sprintf('Phase I De-Embedded Arc/Diagonal Thru S21 Magnitude - %s GHz', bandResult.band));
    xlabel('Frequency (GHz)');
    ylabel('|S21| (dB)');
    grid on;
    legend('Location', 'best');
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_THRU_S21_Magnitude_ArcDiag_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_THRU_S21_Magnitude_ArcDiag_%s.fig', bandResult.band)));
    close(fig);

    fig = figure('Visible', 'on');
    hold on;
    for k = 1:numel(idxList)
        idx = idxList(k);
        baseLabel = prettify_name(thruNames{idx});
        plot(thruFreq{idx} / 1e9, rad2deg(unwrap(angle(thruMeasuredS21{idx}))), '--', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s Measured', baseLabel));
        plot(thruFreq{idx} / 1e9, rad2deg(unwrap(angle(thruS21{idx}))), '-', 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s De-embedded', baseLabel));
    end
    title(sprintf('Phase I De-Embedded Arc/Diagonal Thru S21 Phase - %s GHz', bandResult.band));
    xlabel('Frequency (GHz)');
    ylabel('Phase (deg)');
    grid on;
    legend('Location', 'best');
    saveas(fig, fullfile(figureDir, sprintf('PhaseI_THRU_S21_Phase_ArcDiag_%s.jpg', bandResult.band)));
    savefig(fig, fullfile(figureDir, sprintf('PhaseI_THRU_S21_Phase_ArcDiag_%s.fig', bandResult.band)));
    close(fig);
end
end

function label = prettify_name(rawName)
label = strrep(rawName, '_', ' ');
end

function plot_embedding_comparison(figureDir, bandResult)
freqGHz = bandResult.embedding_diagnostics.freq / 1e9;
e12 = bandResult.embedding_diagnostics.embedding12.S;
e34 = bandResult.embedding_diagnostics.embedding34.S;

s11_12 = squeeze(e12(1, 1, :));
s21_12 = squeeze(e12(2, 1, :));
s22_12 = squeeze(e12(2, 2, :));
s11_34 = squeeze(e34(1, 1, :));
s21_34 = squeeze(e34(2, 1, :));
s22_34 = squeeze(e34(2, 2, :));

sigma12 = compute_sigma_max_local(e12);
sigma34 = compute_sigma_max_local(e34);

fig = figure('Visible', 'on');
subplot(2, 2, 1);
plot(freqGHz, 20 * log10(abs(s11_12)), 'LineWidth', 1.6, 'DisplayName', 'E12'); hold on;
plot(freqGHz, 20 * log10(abs(s11_34)), '--', 'LineWidth', 1.6, 'DisplayName', 'E34');
title('Embedding S11 Magnitude');
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
grid on;
legend('Location', 'best');

subplot(2, 2, 2);
plot(freqGHz, 20 * log10(abs(s21_12)), 'LineWidth', 1.6, 'DisplayName', 'E12'); hold on;
plot(freqGHz, 20 * log10(abs(s21_34)), '--', 'LineWidth', 1.6, 'DisplayName', 'E34');
title('Embedding S21 Magnitude');
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
grid on;
legend('Location', 'best');

subplot(2, 2, 3);
plot(freqGHz, 20 * log10(abs(s22_12)), 'LineWidth', 1.6, 'DisplayName', 'E12'); hold on;
plot(freqGHz, 20 * log10(abs(s22_34)), '--', 'LineWidth', 1.6, 'DisplayName', 'E34');
title('Embedding S22 Magnitude');
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
grid on;
legend('Location', 'best');

subplot(2, 2, 4);
plot(freqGHz, sigma12, 'LineWidth', 1.6, 'DisplayName', 'E12'); hold on;
plot(freqGHz, sigma34, '--', 'LineWidth', 1.6, 'DisplayName', 'E34');
title('Embedding \sigma_{max}');
xlabel('Frequency (GHz)');
ylabel('\sigma_{max}');
grid on;
legend('Location', 'best');

sgtitle(sprintf('Phase I Embedding Comparison - %s GHz', bandResult.band));
saveas(fig, fullfile(figureDir, sprintf('PhaseI_EmbeddingComparison_%s.jpg', bandResult.band)));
savefig(fig, fullfile(figureDir, sprintf('PhaseI_EmbeddingComparison_%s.fig', bandResult.band)));
close(fig);
end

function sigmaMax = compute_sigma_max_local(sParams)
numFreq = size(sParams, 3);
sigmaMax = zeros(numFreq, 1);
for idx = 1:numFreq
    sigmaMax(idx) = max(svd(sParams(:, :, idx)));
end
end

function plot_refinement_diagnostics(figureDir, bandResult)
freqGHz = bandResult.embedding_diagnostics.freq / 1e9;
diag = bandResult.refine_diagnostics;

fig = figure('Visible', 'on');
subplot(2, 2, 1);
plot(freqGHz, diag.abs_s11_shift, 'LineWidth', 1.6);
title('|\Delta S11|');
xlabel('Frequency (GHz)');
ylabel('Magnitude');
grid on;

subplot(2, 2, 2);
plot(freqGHz, diag.abs_s22_shift, 'LineWidth', 1.6);
title('|\Delta S22|');
xlabel('Frequency (GHz)');
ylabel('Magnitude');
grid on;

subplot(2, 2, 3);
plot(freqGHz, diag.short_mag_improvement, 'LineWidth', 1.6, 'DisplayName', 'Short'); hold on;
plot(freqGHz, diag.open_mag_improvement, 'LineWidth', 1.6, 'DisplayName', 'Open');
plot(freqGHz, diag.load_mag_improvement, 'LineWidth', 1.6, 'DisplayName', 'Load');
title('Improvement Metrics');
xlabel('Frequency (GHz)');
ylabel('Reduction');
grid on;
legend('Location', 'best');

subplot(2, 2, 4);
plot(freqGHz, diag.sigma_change, 'LineWidth', 1.6);
title('\Delta \sigma_{max}');
xlabel('Frequency (GHz)');
ylabel('Change');
grid on;

sgtitle(sprintf('Phase I Refinement Diagnostics - %s GHz', bandResult.band));
saveas(fig, fullfile(figureDir, sprintf('PhaseI_RefinementDiagnostics_%s.jpg', bandResult.band)));
savefig(fig, fullfile(figureDir, sprintf('PhaseI_RefinementDiagnostics_%s.fig', bandResult.band)));
close(fig);
end
