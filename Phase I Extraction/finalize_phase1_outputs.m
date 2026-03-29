function results = finalize_phase1_outputs(config, results)
%FINALIZE_PHASE1_OUTPUTS Final selection + global fitted outputs for Phase I.
% Keeps extraction results unchanged and adds:
%   1) Per-band winner selection between P1/P2 for OPEN/SHORT/LOAD
%   2) Per-band fitted standards exported from a global stitched fit
%   3) Stitched fit plots with bandwise traces and simulation overlays

figureDir = fullfile(fileparts(config.touchstone_output_dir), 'figures');
notesDir = fullfile(fileparts(config.touchstone_output_dir), 'notes');
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end
if ~exist(notesDir, 'dir')
    mkdir(notesDir);
end

stdList = {'OPEN', 'SHORT', 'LOAD'};
fitCfg = get_fit_config(config);

% Pass 1: winner selection without global-fit consistency metric
for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    selection = struct();

    for idxStd = 1:numel(stdList)
        stdName = stdList{idxStd};
        p1Key = find_output_key(bandResult.outputs, stdName, 'P1');
        p2Key = find_output_key(bandResult.outputs, stdName, 'P2');
        if isempty(p1Key) || isempty(p2Key)
            continue;
        end

        p1 = bandResult.outputs.(p1Key);
        p2 = bandResult.outputs.(p2Key);
        simRef = get_sim_reference(config, stdName, p1.freq);

        p1Metrics = score_candidate_metrics(stdName, p1, bandResult.embedding, simRef, 1, []);
        p2Metrics = score_candidate_metrics(stdName, p2, bandResult.embedding, simRef, 2, []);

        [winnerPort, reason] = select_winner(p1Metrics, p2Metrics);
        if winnerPort == 1
            winner = p1;
            winnerKey = p1Key;
        else
            winner = p2;
            winnerKey = p2Key;
        end

        s = struct();
        s.freq = winner.freq(:);
        s.gamma_winner = winner.gamma(:);
        s.gamma_fit = winner.gamma(:); % replaced after global fit
        s.winner_port = sprintf('P%d', winnerPort);
        s.winner_output_key = winnerKey;
        s.metrics_p1 = p1Metrics;
        s.metrics_p2 = p2Metrics;
        s.selection_reason = reason;
        s.sim_reference = simRef;
        s.candidate_p1_key = p1Key;
        s.candidate_p2_key = p2Key;
        selection.(stdName) = s;
    end

    results.bands{idxBand}.final_selection = selection;
end

% Build provisional global fits and rerank with global-fit consistency.
globalFits = build_global_fits(config, results, stdList, fitCfg, 'provisional');

for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    selection = bandResult.final_selection;
    stdFields = fieldnames(selection);

    for idxStd = 1:numel(stdFields)
        stdName = stdFields{idxStd};
        s = selection.(stdName);
        p1 = bandResult.outputs.(s.candidate_p1_key);
        p2 = bandResult.outputs.(s.candidate_p2_key);
        simRef = get_sim_reference(config, stdName, p1.freq);

        if isfield(globalFits, stdName) && ~isempty(globalFits.(stdName).freq)
            globalRef = interp1(globalFits.(stdName).freq, globalFits.(stdName).gamma_fit, p1.freq(:), 'linear', 'extrap');
        else
            globalRef = [];
        end

        p1Metrics = score_candidate_metrics(stdName, p1, bandResult.embedding, simRef, 1, globalRef);
        p2Metrics = score_candidate_metrics(stdName, p2, bandResult.embedding, simRef, 2, globalRef);
        [winnerPort, reason] = select_winner(p1Metrics, p2Metrics);

        if winnerPort == 1
            winner = p1;
            winnerKey = s.candidate_p1_key;
        else
            winner = p2;
            winnerKey = s.candidate_p2_key;
        end

        s.freq = winner.freq(:);
        s.gamma_winner = winner.gamma(:);
        s.gamma_fit = winner.gamma(:); % replaced after final global fit
        s.winner_port = sprintf('P%d', winnerPort);
        s.winner_output_key = winnerKey;
        s.metrics_p1 = p1Metrics;
        s.metrics_p2 = p2Metrics;
        s.selection_reason = reason;
        selection.(stdName) = s;
    end

    results.bands{idxBand}.final_selection = selection;
end

% Final global fits based on reranked winners.
globalFits = build_global_fits(config, results, stdList, fitCfg, 'final');

% Export per-band standards from global fit and write notes.
for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    selection = bandResult.final_selection;
    stdFields = fieldnames(selection);

    for idxStd = 1:numel(stdFields)
        stdName = stdFields{idxStd};
        s = selection.(stdName);

        if isfield(globalFits, stdName) && ~isempty(globalFits.(stdName).freq)
            gammaFit = interp1(globalFits.(stdName).freq, globalFits.(stdName).gamma_fit, s.freq(:), 'linear', 'extrap');
        else
            gammaFit = s.gamma_winner;
        end

        outputBase = sprintf('FINALFIT_%s_%s', stdName, bandResult.band);
        write_touchstone_1port(fullfile(config.touchstone_output_dir, [outputBase '.s1p']), s.freq, gammaFit, config.z0);

        finalResult = s;
        finalResult.gamma_fit = gammaFit;
        save(fullfile(config.mat_output_dir, [outputBase '.mat']), 'finalResult');
        selection.(stdName) = finalResult;
    end

    results.bands{idxBand}.final_selection = selection;
    write_selection_notes(notesDir, bandResult.band, selection);
end

plot_stitched_final_fits(config, results, globalFits, figureDir);
plot_stitched_thru_fits(results, figureDir, fitCfg);

% Parallel path: global winner (single port choice across all bands) for each standard.
[results, globalWinner, globalWinnerFits] = run_global_winner_path(config, results, stdList, fitCfg, notesDir);
plot_stitched_final_fits_global(config, results, globalWinner, globalWinnerFits, fitCfg, figureDir);
results = run_global_winner_thru_path(config, results, fitCfg, notesDir, figureDir);

fprintf('Phase I finalization complete (winner selection + global-fitted exports + stitched fit plots).\n');
end

function fitCfg = get_fit_config(config)
fitCfg = struct();
fitCfg.overlap_windows_hz = [63e9, 70e9; 112e9, 118e9];
fitCfg.lambda = 1800;
fitCfg.thru_lambda = 1400;
fitCfg.max_delay_s = 8e-12;
fitCfg.max_gain_db = 2.0;
fitCfg.min_overlap_points = 12;
fitCfg.robust_trim_frac = 0.15;
fitCfg.quality_phase_scale = 0.7;
fitCfg.quality_mag_scale = 0.45;

if isfield(config, 'finalfit_overlap_windows_hz')
    fitCfg.overlap_windows_hz = config.finalfit_overlap_windows_hz;
end
if isfield(config, 'finalfit_lambda')
    fitCfg.lambda = config.finalfit_lambda;
end
if isfield(config, 'finalfit_thru_lambda')
    fitCfg.thru_lambda = config.finalfit_thru_lambda;
end
if isfield(config, 'finalfit_max_delay_s')
    fitCfg.max_delay_s = config.finalfit_max_delay_s;
end
if isfield(config, 'finalfit_max_gain_db')
    fitCfg.max_gain_db = config.finalfit_max_gain_db;
end
if isfield(config, 'finalfit_min_overlap_points')
    fitCfg.min_overlap_points = config.finalfit_min_overlap_points;
end
if isfield(config, 'finalfit_robust_trim_frac')
    fitCfg.robust_trim_frac = config.finalfit_robust_trim_frac;
end
if isfield(config, 'finalfit_quality_phase_scale')
    fitCfg.quality_phase_scale = config.finalfit_quality_phase_scale;
end
if isfield(config, 'finalfit_quality_mag_scale')
    fitCfg.quality_mag_scale = config.finalfit_quality_mag_scale;
end
end

function key = find_output_key(outputs, stdName, portLabel)
fields = fieldnames(outputs);
key = '';
for idx = 1:numel(fields)
    f = upper(fields{idx});
    if contains(f, ['EXTRACTED_' upper(stdName) '_' upper(portLabel) '_'])
        key = fields{idx};
        return;
    end
end
end

function simRef = get_sim_reference(config, stdName, freq)
simRef = [];
if ~isfield(config, 'simulation_priors') || ~isfield(config.simulation_priors, 'enabled') || ~config.simulation_priors.enabled
    return;
end
stdKey = lower(stdName);
if ~isfield(config.simulation_priors.ref, stdKey)
    return;
end
ref = config.simulation_priors.ref.(stdKey);
simRef = interp1(ref.freq(:), ref.gamma(:), freq(:), 'linear', 'extrap');
end

function metrics = score_candidate_metrics(stdName, out, embedding, simRef, portNum, globalRef)
gamma = out.gamma(:);
freq = out.freq(:);
measured = out.measured_gamma(:);
sEmbed = embedding.S;

metrics = struct();
metrics.passivity = passivity_metric(stdName, gamma);
metrics.forward_error = forward_error_metric(gamma, measured, sEmbed, portNum);
metrics.global_fit = global_fit_metric(gamma, globalRef);
metrics.smoothness = smoothness_metric(gamma, freq);
metrics.fit_consistency = fit_consistency_metric(freq, gamma);
metrics.sim_closeness = simulation_metric(gamma, simRef);
end

function val = passivity_metric(stdName, gamma)
switch upper(stdName)
    case {'OPEN', 'SHORT'}
        val = mean(abs(abs(gamma) - 1), 'omitnan');
    otherwise
        val = mean(abs(gamma), 'omitnan');
end
end

function val = forward_error_metric(gammaActual, gammaMeasured, sEmbed, portNum)
numFreq = numel(gammaActual);
errs = zeros(numFreq, 1);
for idx = 1:numFreq
    s11 = sEmbed(1, 1, idx);
    s12 = sEmbed(1, 2, idx);
    s21 = sEmbed(2, 1, idx);
    s22 = sEmbed(2, 2, idx);
    if portNum == 1
        pred = s11 + (s12 * s21 * gammaActual(idx)) / (1 - s22 * gammaActual(idx));
    else
        pred = s22 + (s12 * s21 * gammaActual(idx)) / (1 - s11 * gammaActual(idx));
    end
    errs(idx) = abs(pred - gammaMeasured(idx));
end
val = mean(errs, 'omitnan');
end

function val = global_fit_metric(gamma, globalRef)
if isempty(globalRef)
    val = inf;
    return;
end
val = mean(abs(gamma(:) - globalRef(:)), 'omitnan');
end

function val = smoothness_metric(gamma, freq)
if numel(gamma) < 3
    val = 0;
    return;
end
f = freq(:);
g = gamma(:);
d1 = diff(g) ./ max(diff(f), eps);
d2 = diff(d1) ./ max(diff(f(2:end)), eps);
val = mean(abs(d2), 'omitnan');
end

function val = simulation_metric(gamma, simRef)
if isempty(simRef)
    val = inf;
    return;
end
val = mean(abs(gamma(:) - simRef(:)), 'omitnan');
end

function [winnerPort, reason] = select_winner(m1, m2)
% Priority order:
% passivity -> forward_error -> global_fit -> smoothness -> fit_consistency -> sim_closeness
metrics = {'passivity', 'forward_error', 'global_fit', 'smoothness', 'fit_consistency', 'sim_closeness'};
tols = [0.03, 0.05, 0.08, 0.10, 0.10, 0.10]; % relative tie thresholds

winnerPort = 1;
reason = 'P1 default tie-break';

for idx = 1:numel(metrics)
    k = metrics{idx};
    a = m1.(k);
    b = m2.(k);
    if ~isfinite(a) && ~isfinite(b)
        continue;
    end
    [better, isTie] = compare_metric(a, b, tols(idx));
    if ~isTie
        winnerPort = better;
        reason = sprintf('Selected by %s (P%d better)', k, better);
        return;
    end
end
end

function val = fit_consistency_metric(freq, gamma)
gammaFit = smooth_complex_pchip(freq, gamma);
val = mean(abs(gamma - gammaFit), 'omitnan');
end

function [betterPort, isTie] = compare_metric(a, b, relTol)
if ~isfinite(a) && isfinite(b)
    betterPort = 2;
    isTie = false;
    return;
end
if isfinite(a) && ~isfinite(b)
    betterPort = 1;
    isTie = false;
    return;
end
if ~isfinite(a) && ~isfinite(b)
    betterPort = 1;
    isTie = true;
    return;
end

scale = max(max(abs([a, b])), eps);
delta = (a - b) / scale;
if abs(delta) <= relTol
    isTie = true;
    betterPort = 1;
    return;
end
isTie = false;
if a < b
    betterPort = 1;
else
    betterPort = 2;
end
end

function gammaFit = smooth_complex_pchip(freq, gamma)
f = freq(:);
g = gamma(:);
numPts = numel(f);
if numPts < 8
    gammaFit = g;
    return;
end

knotCount = max(12, round(numPts / 12));
knotIdx = unique(round(linspace(1, numPts, knotCount)));
kf = f(knotIdx);
kr = real(g(knotIdx));
ki = imag(g(knotIdx));

fitReal = interp1(kf, kr, f, 'pchip', 'extrap');
fitImag = interp1(kf, ki, f, 'pchip', 'extrap');
gammaFit = fitReal + 1i * fitImag;
end

function write_selection_notes(notesDir, bandName, selection)
filename = fullfile(notesDir, sprintf('PhaseI_FinalSelection_%s.txt', bandName));
fid = fopen(filename, 'w');
if fid == -1
    warning('Unable to write final selection note: %s', filename);
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Phase I Final Selection Summary - %s GHz band\n\n', bandName);
stdList = fieldnames(selection);
for idx = 1:numel(stdList)
    key = stdList{idx};
    s = selection.(key);
    fprintf(fid, '%s:\n', key);
    fprintf(fid, '  Winner: %s\n', s.winner_port);
    fprintf(fid, '  Reason: %s\n', s.selection_reason);
    fprintf(fid, ['  P1 metrics: passivity=%.6g, forward=%.6g, globalfit=%.6g, ' ...
        'smooth=%.6g, fit=%.6g, sim=%.6g\n'], ...
        s.metrics_p1.passivity, s.metrics_p1.forward_error, s.metrics_p1.global_fit, ...
        s.metrics_p1.smoothness, s.metrics_p1.fit_consistency, s.metrics_p1.sim_closeness);
    fprintf(fid, ['  P2 metrics: passivity=%.6g, forward=%.6g, globalfit=%.6g, ' ...
        'smooth=%.6g, fit=%.6g, sim=%.6g\n\n'], ...
        s.metrics_p2.passivity, s.metrics_p2.forward_error, s.metrics_p2.global_fit, ...
        s.metrics_p2.smoothness, s.metrics_p2.fit_consistency, s.metrics_p2.sim_closeness);
end
end

function globalFits = build_global_fits(config, results, stdList, fitCfg, tag)
globalFits = struct();
for idxStd = 1:numel(stdList)
    stdName = stdList{idxStd};
    entries = collect_std_entries(results, stdName);
    if isempty(entries)
        continue;
    end

    alignOpts = struct('use_gain', true, 'lambda', fitCfg.lambda);
    fitOut = build_global_stitched_fit(entries, 'gamma_winner', fitCfg, alignOpts);
    globalFits.(stdName) = fitOut;

    if isfield(config, 'mat_output_dir') && ~isempty(config.mat_output_dir)
        outFile = fullfile(config.mat_output_dir, sprintf('FINALFIT_%s_STITCHED_%s.mat', stdName, upper(tag)));
        save(outFile, 'fitOut');
    end
end
end

function plot_stitched_final_fits(config, results, globalFits, figureDir)
stdList = {'OPEN', 'SHORT', 'LOAD'};

for idxStd = 1:numel(stdList)
    stdName = stdList{idxStd};
    entries = collect_std_entries(results, stdName);
    if isempty(entries)
        continue;
    end
    if ~isfield(globalFits, stdName)
        continue;
    end

    fitOut = globalFits.(stdName);
    sim = collect_sim_overlay(config, stdName, fitOut.freq);

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on;
    for idx = 1:numel(entries)
        plot(entries{idx}.freq / 1e9, 20 * log10(abs(entries{idx}.gamma_winner)), '-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s %s winner', entries{idx}.band, entries{idx}.winner_port));
    end
    plot(fitOut.freq / 1e9, 20 * log10(abs(fitOut.gamma_fit)), 'k-', 'LineWidth', 2.4, ...
        'DisplayName', sprintf('%s stitched smooth fit', stdName));
    if ~isempty(sim)
        plot(sim.freq / 1e9, 20 * log10(abs(sim.gamma)), '--', 'LineWidth', 2.0, ...
            'DisplayName', sprintf('%s simulated de-embedded', stdName));
    end
    title(sprintf('Phase I Stitched %s Magnitude (Bandwise + Smooth Fit)', stdName));
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    grid on;
    legend('Location', 'best');
    apply_axes_style(gca);
    save_figure(fig, fullfile(figureDir, sprintf('PhaseI_StitchedFit_%s_Magnitude', stdName)));
    close(fig);

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on;
    for idx = 1:numel(entries)
        plot(entries{idx}.freq / 1e9, rad2deg(unwrap(angle(entries{idx}.gamma_winner))), '-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s %s winner', entries{idx}.band, entries{idx}.winner_port));
    end
    plot(fitOut.freq / 1e9, rad2deg(unwrap(angle(fitOut.gamma_fit))), 'k-', 'LineWidth', 2.4, ...
        'DisplayName', sprintf('%s stitched smooth fit', stdName));
    if ~isempty(sim)
        plot(sim.freq / 1e9, rad2deg(unwrap(angle(sim.gamma))), '--', 'LineWidth', 2.0, ...
            'DisplayName', sprintf('%s simulated de-embedded', stdName));
    end
    title(sprintf('Phase I Stitched %s Phase (Bandwise + Smooth Fit)', stdName));
    xlabel('Frequency (GHz)');
    ylabel('Phase (deg)');
    grid on;
    legend('Location', 'best');
    apply_axes_style(gca);
    save_figure(fig, fullfile(figureDir, sprintf('PhaseI_StitchedFit_%s_Phase', stdName)));
    close(fig);
end
end

function plot_stitched_thru_fits(results, figureDir, fitCfg)
thruMap = containers.Map();
for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    outFields = fieldnames(bandResult.outputs);
    for idx = 1:numel(outFields)
        out = bandResult.outputs.(outFields{idx});
        if ~isfield(out, 'S')
            continue;
        end
        key = regexprep(outFields{idx}, '_\d+_\d+$', '');
        entry = struct();
        entry.band = bandResult.band;
        entry.freq = out.freq(:);
        entry.s21 = squeeze(out.S(2, 1, :));
        if isKey(thruMap, key)
            tmp = thruMap(key);
            tmp{end + 1} = entry;
            thruMap(key) = tmp;
        else
            thruMap(key) = {entry};
        end
    end
end

keysList = sort(thruMap.keys);
for idxKey = 1:numel(keysList)
    k = keysList{idxKey};
    entries = thruMap(k);
    alignOpts = struct('use_gain', false, 'lambda', fitCfg.thru_lambda);
    fitOut = build_global_stitched_fit(entries, 's21', fitCfg, alignOpts);

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on;
    for ii = 1:numel(entries)
        plot(entries{ii}.freq / 1e9, 20 * log10(abs(entries{ii}.s21)), '-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s raw', entries{ii}.band));
    end
    plot(fitOut.freq / 1e9, 20 * log10(abs(fitOut.gamma_fit)), 'k-', 'LineWidth', 2.4, ...
        'DisplayName', 'stitched smooth fit');
    title(sprintf('Phase I Stitched THRU %s S21 Magnitude (Bandwise + Smooth Fit)', strrep(k, '_', ' ')));
    xlabel('Frequency (GHz)');
    ylabel('|S21| (dB)');
    grid on;
    legend('Location', 'best');
    apply_axes_style(gca);
    save_figure(fig, fullfile(figureDir, sprintf('PhaseI_StitchedFit_THRU_%s_Magnitude', k)));
    close(fig);

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on;
    for ii = 1:numel(entries)
        plot(entries{ii}.freq / 1e9, rad2deg(unwrap(angle(entries{ii}.s21))), '-', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%s raw', entries{ii}.band));
    end
    plot(fitOut.freq / 1e9, rad2deg(unwrap(angle(fitOut.gamma_fit))), 'k-', 'LineWidth', 2.4, ...
        'DisplayName', 'stitched smooth fit');
    title(sprintf('Phase I Stitched THRU %s S21 Phase (Bandwise + Smooth Fit)', strrep(k, '_', ' ')));
    xlabel('Frequency (GHz)');
    ylabel('Phase (deg)');
    grid on;
    legend('Location', 'best');
    apply_axes_style(gca);
    save_figure(fig, fullfile(figureDir, sprintf('PhaseI_StitchedFit_THRU_%s_Phase', k)));
    close(fig);
end
end

function entries = collect_std_entries(results, stdName)
entries = {};
for idxBand = 1:numel(results.bands)
    b = results.bands{idxBand};
    if ~isfield(b, 'final_selection') || ~isfield(b.final_selection, stdName)
        continue;
    end
    fs = b.final_selection.(stdName);
    entry = struct();
    entry.band = b.band;
    entry.freq = fs.freq(:);
    entry.gamma_winner = fs.gamma_winner(:);
    entry.winner_port = fs.winner_port;
    entries{end + 1} = entry; %#ok<AGROW>
end
end

function sim = collect_sim_overlay(config, stdName, freqTarget)
sim = [];
if ~isfield(config, 'simulation_priors') || ~isfield(config.simulation_priors, 'enabled') || ~config.simulation_priors.enabled
    return;
end
stdKey = lower(stdName);
if ~isfield(config.simulation_priors.ref, stdKey)
    return;
end
ref = config.simulation_priors.ref.(stdKey);
sim = struct();
sim.freq = freqTarget(:);
sim.gamma = interp1(ref.freq(:), ref.gamma(:), freqTarget(:), 'linear', 'extrap');
end

function fitOut = build_global_stitched_fit(entries, fieldName, fitCfg, alignOpts)
% Build a true global fit:
% 1) overlap-based robust band alignment (gain + phi + tau)
% 2) balanced weighting per band
% 3) complex-domain smooth solve over the full union frequency grid

if nargin < 4 || isempty(alignOpts)
    alignOpts = struct('use_gain', true, 'lambda', fitCfg.lambda);
end

startFreq = zeros(numel(entries), 1);
for idx = 1:numel(entries)
    startFreq(idx) = entries{idx}.freq(1);
end
[~, order] = sort(startFreq);
entries = entries(order);

allFreq = [];
allData = [];
allW = [];
alignInfo = struct('band', {}, 'gain', {}, 'phi_rad', {}, 'tau_s', {});

assembledFreq = [];
assembledData = [];

for idx = 1:numel(entries)
    f = entries{idx}.freq(:);
    g = entries{idx}.(fieldName)(:);

    [gain, phi, tau, quality] = estimate_alignment(f, g, assembledFreq, assembledData, fitCfg, alignOpts);
    gAligned = gain .* g .* exp(1i * (phi + 2 * pi * (f - mean(f)) * tau));

    alignInfo(idx).band = entries{idx}.band; %#ok<AGROW>
    alignInfo(idx).gain = gain; %#ok<AGROW>
    alignInfo(idx).phi_rad = phi; %#ok<AGROW>
    alignInfo(idx).tau_s = tau; %#ok<AGROW>
    alignInfo(idx).quality = quality; %#ok<AGROW>

    bandWeight = 1 / max(numel(entries), 1);
    pointWeight = bandWeight / max(numel(f), 1);

    allFreq = [allFreq; f]; %#ok<AGROW>
    allData = [allData; gAligned]; %#ok<AGROW>
    allW = [allW; pointWeight * ones(numel(f), 1)]; %#ok<AGROW>

    assembledFreq = [assembledFreq; f]; %#ok<AGROW>
    assembledData = [assembledData; gAligned]; %#ok<AGROW>
    [assembledFreq, sortIdx] = sort(assembledFreq);
    assembledData = assembledData(sortIdx);
end

[freqU, ~, ic] = unique(allFreq);
wU = accumarray(ic, allW);
yU = accumarray(ic, allW .* allData) ./ max(wU, eps);

gammaFit = smooth_complex_global(yU, wU, alignOpts.lambda);

fitOut = struct();
fitOut.freq = freqU(:);
fitOut.gamma_fit = gammaFit(:);
fitOut.alignment = alignInfo;
end

function [gain, phi, tau, quality] = estimate_alignment(f, g, refF, refG, fitCfg, alignOpts)
gain = 1;
phi = 0;
tau = 0;
quality = 0;
if isempty(refF) || isempty(refG)
    return;
end

[f, ~, fIdx] = unique(f(:), 'sorted');
if numel(f) < 2
    return;
end
g = accumarray(fIdx, g(:), [], @(x) mean(x, 'omitnan'));

% interp1 requires strictly unique sample points.
[refF, ~, refIdx] = unique(refF(:), 'sorted');
if numel(refF) < 2
    return;
end
refG = accumarray(refIdx, refG(:), [], @(x) mean(x, 'omitnan'));

% Find overlap region between current band and already assembled data.
fMin = min(f);
fMax = max(f);
win = fitCfg.overlap_windows_hz;
haveWindow = false;
ovLo = -inf;
ovHi = inf;
for idx = 1:size(win, 1)
    lo = max([fMin, min(refF), win(idx, 1)]);
    hi = min([fMax, max(refF), win(idx, 2)]);
    if hi > lo
        ovLo = lo;
        ovHi = hi;
        haveWindow = true;
        break;
    end
end

if ~haveWindow
    % Fallback: constant phase anchor at nearest frequency.
    [~, iRef] = min(abs(refF - f(1)));
    if abs(g(1)) > 0
        r = refG(iRef) / g(1);
        gain = abs(r);
        phi = angle(r);
    else
        gain = 1;
        phi = 0;
    end
    if ~alignOpts.use_gain
        gain = 1;
    end
    gain = clamp_gain(gain, fitCfg.max_gain_db);
    tau = 0;
    quality = 0.25;
    return;
end

grid = linspace(ovLo, ovHi, 200).';
gRef = interp1(refF, refG, grid, 'linear', NaN);
gCur = interp1(f, g, grid, 'linear', NaN);
valid = ~isnan(gRef) & ~isnan(gCur) & (abs(gCur) > 1e-5) & (abs(gRef) > 1e-5);
if nnz(valid) < fitCfg.min_overlap_points
    [~, iRef] = min(abs(refF - f(1)));
    if abs(g(1)) > 0
        r = refG(iRef) / g(1);
        gain = abs(r);
        phi = angle(r);
    else
        gain = 1;
        phi = 0;
    end
    if ~alignOpts.use_gain
        gain = 1;
    end
    gain = clamp_gain(gain, fitCfg.max_gain_db);
    tau = 0;
    quality = 0.25;
    return;
end

grid = grid(valid);
ratio = gRef(valid) ./ gCur(valid);
magRatio = abs(ratio);
phaseRatio = angle(ratio);
finiteMask = isfinite(magRatio) & isfinite(phaseRatio) & magRatio > 1e-6;
magRatio = magRatio(finiteMask);
phaseRatio = phaseRatio(finiteMask);
grid = grid(finiteMask);

if isempty(magRatio)
    gain = 1;
    phi = 0;
    tau = 0;
    quality = 0;
    return;
end

logMag = log(magRatio);
trimN = floor(numel(logMag) * fitCfg.robust_trim_frac);
if trimN > 0 && (2 * trimN) < numel(logMag)
    sorted = sort(logMag);
    lo = sorted(trimN + 1);
    hi = sorted(end - trimN);
    keep = logMag >= lo & logMag <= hi;
    logMag = logMag(keep);
    phaseRatio = phaseRatio(keep);
    grid = grid(keep);
end

if isempty(logMag)
    gainRaw = 1;
else
    gainRaw = exp(median(logMag, 'omitnan'));
end

phaseDiff = unwrap(phaseRatio);
x = 2 * pi * (grid - mean(grid));
if numel(x) >= 2
    p = polyfit(x, phaseDiff, 1);
    tauRaw = p(1);
    phiRaw = p(2);
else
    tauRaw = 0;
    phiRaw = median(phaseDiff, 'omitnan');
end

resid = wrap_to_pi(phaseDiff - (tauRaw * x + phiRaw));
phaseSpread = std(resid, 'omitnan');
magSpread = std(logMag, 'omitnan');
coverage = numel(logMag) / max(numel(valid), 1);
quality = coverage * exp(-phaseSpread / max(fitCfg.quality_phase_scale, eps)) * ...
    exp(-magSpread / max(fitCfg.quality_mag_scale, eps));
quality = min(max(quality, 0), 1);

if alignOpts.use_gain
    gain = 1 + quality * (gainRaw - 1);
else
    gain = 1;
end
phi = quality * phiRaw;
tau = quality * tauRaw;

gain = clamp_gain(gain, fitCfg.max_gain_db);
tau = min(max(tau, -fitCfg.max_delay_s), fitCfg.max_delay_s);
end

function [results, globalWinner, globalWinnerFits] = run_global_winner_path(config, results, stdList, fitCfg, notesDir)
globalWinner = struct();
globalWinnerFits = struct();

for idxStd = 1:numel(stdList)
    stdName = stdList{idxStd};
    [bestPort, summary] = pick_global_winner_port(config, results, stdName);
    if isempty(bestPort)
        continue;
    end

    entries = collect_std_entries_fixed_port(results, stdName, bestPort);
    if isempty(entries)
        continue;
    end

    alignOpts = struct('use_gain', true, 'lambda', fitCfg.lambda);
    fitOut = build_global_stitched_fit(entries, 'gamma_winner', fitCfg, alignOpts);
    globalWinnerFits.(stdName) = fitOut;

    globalWinner.(stdName) = struct();
    globalWinner.(stdName).port = bestPort;
    globalWinner.(stdName).summary = summary;

    for idxBand = 1:numel(results.bands)
        bandResult = results.bands{idxBand};
        outKey = find_output_key(bandResult.outputs, stdName, bestPort);
        if isempty(outKey)
            continue;
        end
        out = bandResult.outputs.(outKey);
        gammaFit = interp1(fitOut.freq, fitOut.gamma_fit, out.freq(:), 'linear', 'extrap');

        outputBase = sprintf('FINALFIT_GLOBALWIN_%s_%s', stdName, bandResult.band);
        write_touchstone_1port(fullfile(config.touchstone_output_dir, [outputBase '.s1p']), out.freq(:), gammaFit, config.z0);

        finalResult = struct();
        finalResult.freq = out.freq(:);
        finalResult.gamma_winner = out.gamma(:);
        finalResult.gamma_fit = gammaFit(:);
        finalResult.winner_port = bestPort;
        finalResult.winner_output_key = outKey;
        finalResult.global_summary = summary;
        save(fullfile(config.mat_output_dir, [outputBase '.mat']), 'finalResult');

        if ~isfield(results.bands{idxBand}, 'final_selection_global')
            results.bands{idxBand}.final_selection_global = struct();
        end
        results.bands{idxBand}.final_selection_global.(stdName) = finalResult;
    end
end

write_global_winner_notes(notesDir, globalWinner);
end

function [bestPort, summary] = pick_global_winner_port(config, results, stdName)
bestPort = '';
summary = struct();

metricNames = {'passivity', 'forward_error', 'smoothness', 'fit_consistency', 'sim_closeness'};
metricWeights = [4.0, 3.0, 1.2, 1.0, 0.6];

scoreP1 = 0;
scoreP2 = 0;
usedBands = {};
bandDetail = struct('band', {}, 'p1', {}, 'p2', {}, 'p1_norm', {}, 'p2_norm', {});

for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    p1Key = find_output_key(bandResult.outputs, stdName, 'P1');
    p2Key = find_output_key(bandResult.outputs, stdName, 'P2');
    if isempty(p1Key) || isempty(p2Key)
        continue;
    end

    p1 = bandResult.outputs.(p1Key);
    p2 = bandResult.outputs.(p2Key);
    simRef = get_sim_reference(config, stdName, p1.freq(:));

    m1 = score_candidate_metrics(stdName, p1, bandResult.embedding, simRef, 1, []);
    m2 = score_candidate_metrics(stdName, p2, bandResult.embedding, simRef, 2, []);

    p1Norm = 0;
    p2Norm = 0;
    for idxMetric = 1:numel(metricNames)
        key = metricNames{idxMetric};
        a = m1.(key);
        b = m2.(key);
        if ~isfinite(a) && ~isfinite(b)
            continue;
        end
        if ~isfinite(a)
            a = 10 * max(abs(b), 1);
        end
        if ~isfinite(b)
            b = 10 * max(abs(a), 1);
        end
        denom = max(min([a, b]), eps);
        p1Norm = p1Norm + metricWeights(idxMetric) * (a / denom);
        p2Norm = p2Norm + metricWeights(idxMetric) * (b / denom);
    end

    scoreP1 = scoreP1 + p1Norm;
    scoreP2 = scoreP2 + p2Norm;
    usedBands{end + 1} = bandResult.band; %#ok<AGROW>

    bd = struct();
    bd.band = bandResult.band;
    bd.p1 = m1;
    bd.p2 = m2;
    bd.p1_norm = p1Norm;
    bd.p2_norm = p2Norm;
    bandDetail(end + 1) = bd; %#ok<AGROW>
end

if isempty(usedBands)
    return;
end

if scoreP1 <= scoreP2
    bestPort = 'P1';
else
    bestPort = 'P2';
end

summary.metric_names = metricNames;
summary.metric_weights = metricWeights;
summary.used_bands = usedBands;
summary.total_score_p1 = scoreP1;
summary.total_score_p2 = scoreP2;
summary.band_detail = bandDetail;
summary.reason = sprintf('Global winner by aggregated normalized weighted score (%s selected)', bestPort);
end

function entries = collect_std_entries_fixed_port(results, stdName, portLabel)
entries = {};
for idxBand = 1:numel(results.bands)
    b = results.bands{idxBand};
    outKey = find_output_key(b.outputs, stdName, portLabel);
    if isempty(outKey)
        continue;
    end
    out = b.outputs.(outKey);
    entry = struct();
    entry.band = b.band;
    entry.freq = out.freq(:);
    entry.gamma_winner = out.gamma(:);
    entry.winner_port = portLabel;
    entries{end + 1} = entry; %#ok<AGROW>
end
end

function write_global_winner_notes(notesDir, globalWinner)
filename = fullfile(notesDir, 'PhaseI_GlobalWinnerSelection.txt');
fid = fopen(filename, 'w');
if fid == -1
    warning('Unable to write global winner note: %s', filename);
    return;
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'Phase I Global Winner Selection (single port across all bands)\n\n');
stdFields = fieldnames(globalWinner);
for idx = 1:numel(stdFields)
    stdName = stdFields{idx};
    g = globalWinner.(stdName);
    fprintf(fid, '%s:\n', stdName);
    fprintf(fid, '  Winner: %s\n', g.port);
    fprintf(fid, '  Score P1: %.6g\n', g.summary.total_score_p1);
    fprintf(fid, '  Score P2: %.6g\n', g.summary.total_score_p2);
    fprintf(fid, '  Reason: %s\n', g.summary.reason);
    fprintf(fid, '  Bands used: %s\n', strjoin(g.summary.used_bands, ', '));
    for ii = 1:numel(g.summary.band_detail)
        bd = g.summary.band_detail(ii);
        fprintf(fid, '    %s: p1_norm=%.6g, p2_norm=%.6g\n', bd.band, bd.p1_norm, bd.p2_norm);
    end
    fprintf(fid, '\n');
end
end

function plot_stitched_final_fits_global(config, results, globalWinner, globalWinnerFits, fitCfg, figureDir)
stdList = {'OPEN', 'SHORT', 'LOAD'};
for idxStd = 1:numel(stdList)
    stdName = stdList{idxStd};
    if ~isfield(globalWinner, stdName) || ~isfield(globalWinnerFits, stdName)
        continue;
    end
    winnerPort = globalWinner.(stdName).port;
    loserPort = opposite_port_label(winnerPort);
    winnerEntries = collect_std_entries_fixed_port(results, stdName, winnerPort);
    loserEntries = collect_std_entries_fixed_port(results, stdName, loserPort);
    if isempty(winnerEntries)
        continue;
    end
    fitOut = globalWinnerFits.(stdName);
    loserFit = [];
    if ~isempty(loserEntries)
        alignOpts = struct('use_gain', true, 'lambda', fitCfg.lambda);
        loserFit = build_global_stitched_fit(loserEntries, 'gamma_winner', fitCfg, alignOpts);
    end
    sim = collect_sim_overlay(config, stdName, fitOut.freq);

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on;
    for idx = 1:numel(winnerEntries)
        plot(winnerEntries{idx}.freq / 1e9, 20 * log10(abs(winnerEntries{idx}.gamma_winner)), '-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('%s %s raw', winnerEntries{idx}.band, winnerPort));
    end
    for idx = 1:numel(loserEntries)
        plot(loserEntries{idx}.freq / 1e9, 20 * log10(abs(loserEntries{idx}.gamma_winner)), '--', 'LineWidth', 1.0, ...
            'DisplayName', sprintf('%s %s rejected', loserEntries{idx}.band, loserPort));
    end
    plot(fitOut.freq / 1e9, 20 * log10(abs(fitOut.gamma_fit)), 'k-', 'LineWidth', 2.4, ...
        'DisplayName', sprintf('%s global-winner smooth fit', stdName));
    if ~isempty(loserFit)
        plot(loserFit.freq / 1e9, 20 * log10(abs(loserFit.gamma_fit)), '-', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s rejected smooth fit', stdName));
    end
    if ~isempty(sim)
        plot(sim.freq / 1e9, 20 * log10(abs(sim.gamma)), '--', 'LineWidth', 2.0, ...
            'DisplayName', sprintf('%s simulated de-embedded', stdName));
    end
    title(sprintf('Phase I Stitched %s Magnitude (Global Winner %s)', stdName, winnerPort));
    xlabel('Frequency (GHz)');
    ylabel('Magnitude (dB)');
    grid on;
    legend('Location', 'best');
    apply_axes_style(gca);
    save_figure(fig, fullfile(figureDir, sprintf('PhaseI_StitchedFit_GLOBALWIN_%s_Magnitude', stdName)));
    close(fig);

    fig = figure('Visible', 'off', 'Color', 'w');
    hold on;
    for idx = 1:numel(winnerEntries)
        plot(winnerEntries{idx}.freq / 1e9, rad2deg(unwrap(angle(winnerEntries{idx}.gamma_winner))), '-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('%s %s raw', winnerEntries{idx}.band, winnerPort));
    end
    for idx = 1:numel(loserEntries)
        plot(loserEntries{idx}.freq / 1e9, rad2deg(unwrap(angle(loserEntries{idx}.gamma_winner))), '--', 'LineWidth', 1.0, ...
            'DisplayName', sprintf('%s %s rejected', loserEntries{idx}.band, loserPort));
    end
    plot(fitOut.freq / 1e9, rad2deg(unwrap(angle(fitOut.gamma_fit))), 'k-', 'LineWidth', 2.4, ...
        'DisplayName', sprintf('%s global-winner smooth fit', stdName));
    if ~isempty(loserFit)
        plot(loserFit.freq / 1e9, rad2deg(unwrap(angle(loserFit.gamma_fit))), '-', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.8, ...
            'DisplayName', sprintf('%s rejected smooth fit', stdName));
    end
    if ~isempty(sim)
        plot(sim.freq / 1e9, rad2deg(unwrap(angle(sim.gamma))), '--', 'LineWidth', 2.0, ...
            'DisplayName', sprintf('%s simulated de-embedded', stdName));
    end
    title(sprintf('Phase I Stitched %s Phase (Global Winner %s)', stdName, winnerPort));
    xlabel('Frequency (GHz)');
    ylabel('Phase (deg)');
    grid on;
    legend('Location', 'best');
    apply_axes_style(gca);
    save_figure(fig, fullfile(figureDir, sprintf('PhaseI_StitchedFit_GLOBALWIN_%s_Phase', stdName)));
    close(fig);
end
end

function portLabel = opposite_port_label(portLabelIn)
if strcmpi(portLabelIn, 'P1')
    portLabel = 'P2';
else
    portLabel = 'P1';
end
end

function results = run_global_winner_thru_path(config, results, fitCfg, notesDir, figureDir)
thruMap = collect_thru_map(results);
if thruMap.Count == 0
    return;
end

keysList = sort(thruMap.keys);
noteFile = fullfile(notesDir, 'PhaseI_GlobalWinnerThruSelection.txt');
fid = fopen(noteFile, 'w');
if fid ~= -1
    c = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, 'Phase I GLOBALWIN THRU smoothing summary\n\n');
end

for idxKey = 1:numel(keysList)
    key = keysList{idxKey};
    entries = thruMap(key);

    alignTx = struct('use_gain', false, 'lambda', fitCfg.thru_lambda);
    fit21 = build_global_stitched_fit(entries, 's21', fitCfg, alignTx);
    fit12 = build_global_stitched_fit(entries, 's12', fitCfg, alignTx);

    plot_globalwin_thru(figureDir, key, entries, fit21);

    for idxBand = 1:numel(results.bands)
        bandResult = results.bands{idxBand};
        outFields = fieldnames(bandResult.outputs);
        for idxField = 1:numel(outFields)
            fieldName = outFields{idxField};
            baseKey = regexprep(fieldName, '_\d+_\d+$', '');
            if ~strcmpi(baseKey, key)
                continue;
            end
            out = bandResult.outputs.(fieldName);
            if ~isfield(out, 'S')
                continue;
            end

            freq = out.freq(:);
            sFit = out.S;
            s21Fit = interp1(fit21.freq, fit21.gamma_fit, freq, 'linear', 'extrap');
            s12Fit = interp1(fit12.freq, fit12.gamma_fit, freq, 'linear', 'extrap');
            sFit(2, 1, :) = reshape(s21Fit, 1, 1, []);
            sFit(1, 2, :) = reshape(s12Fit, 1, 1, []);

            outputBase = sprintf('FINALFIT_GLOBALWIN_%s_%s', key, bandResult.band);
            write_touchstone_2port(fullfile(config.touchstone_output_dir, [outputBase '.s2p']), freq, sFit, config.z0);
            finalResult = struct();
            finalResult.freq = freq;
            finalResult.s_fit = sFit;
            finalResult.s_raw = out.S;
            finalResult.fit21 = fit21;
            finalResult.fit12 = fit12;
            finalResult.key = key;
            save(fullfile(config.mat_output_dir, [outputBase '.mat']), 'finalResult');

            if ~isfield(results.bands{idxBand}, 'final_selection_global_thru')
                results.bands{idxBand}.final_selection_global_thru = struct();
            end
            results.bands{idxBand}.final_selection_global_thru.(fieldName) = finalResult;
        end
    end

    if fid ~= -1
        fprintf(fid, '%s:\n', key);
        fprintf(fid, '  fit21 points: %d\n', numel(fit21.freq));
        fprintf(fid, '  fit12 points: %d\n', numel(fit12.freq));
        fprintf(fid, '\n');
    end
end
end

function thruMap = collect_thru_map(results)
thruMap = containers.Map();
for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    outFields = fieldnames(bandResult.outputs);
    for idx = 1:numel(outFields)
        fieldName = outFields{idx};
        out = bandResult.outputs.(fieldName);
        if ~isfield(out, 'S')
            continue;
        end
        key = regexprep(fieldName, '_\d+_\d+$', '');
        entry = struct();
        entry.band = bandResult.band;
        entry.freq = out.freq(:);
        entry.s21 = squeeze(out.S(2, 1, :));
        entry.s12 = squeeze(out.S(1, 2, :));
        if isKey(thruMap, key)
            tmp = thruMap(key);
            tmp{end + 1} = entry;
            thruMap(key) = tmp;
        else
            thruMap(key) = {entry};
        end
    end
end
end

function plot_globalwin_thru(figureDir, key, entries, fit21)
fig = figure('Visible', 'off', 'Color', 'w');
hold on;
for ii = 1:numel(entries)
    plot(entries{ii}.freq / 1e9, 20 * log10(abs(entries{ii}.s21)), '-', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s raw', entries{ii}.band));
end
plot(fit21.freq / 1e9, 20 * log10(abs(fit21.gamma_fit)), 'k-', 'LineWidth', 2.4, ...
    'DisplayName', 'GLOBALWIN smooth fit');
title(sprintf('Phase I GLOBALWIN THRU %s S21 Magnitude', strrep(key, '_', ' ')));
xlabel('Frequency (GHz)');
ylabel('|S21| (dB)');
grid on;
legend('Location', 'best');
apply_axes_style(gca);
save_figure(fig, fullfile(figureDir, sprintf('PhaseI_StitchedFit_GLOBALWIN_THRU_%s_Magnitude', key)));
close(fig);

fig = figure('Visible', 'off', 'Color', 'w');
hold on;
for ii = 1:numel(entries)
    plot(entries{ii}.freq / 1e9, rad2deg(unwrap(angle(entries{ii}.s21))), '-', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s raw', entries{ii}.band));
end
plot(fit21.freq / 1e9, rad2deg(unwrap(angle(fit21.gamma_fit))), 'k-', 'LineWidth', 2.4, ...
    'DisplayName', 'GLOBALWIN smooth fit');
title(sprintf('Phase I GLOBALWIN THRU %s S21 Phase', strrep(key, '_', ' ')));
xlabel('Frequency (GHz)');
ylabel('Phase (deg)');
grid on;
legend('Location', 'best');
apply_axes_style(gca);
save_figure(fig, fullfile(figureDir, sprintf('PhaseI_StitchedFit_GLOBALWIN_THRU_%s_Phase', key)));
close(fig);
end

function g = clamp_gain(g, maxGainDb)
if ~isfinite(g) || g <= 0
    g = 1;
    return;
end
maxLin = 10^(abs(maxGainDb) / 20);
g = min(max(g, 1 / maxLin), maxLin);
end

function y = wrap_to_pi(x)
y = mod(x + pi, 2 * pi) - pi;
end

function x = smooth_complex_global(y, w, lambda)
n = numel(y);
if n <= 4
    x = y;
    return;
end

e = ones(n, 1);
D2 = spdiags([e, -2 * e, e], 0:2, n - 2, n);
W = spdiags(w(:), 0, n, n);
L = D2' * D2;
A = W + lambda * L + 1e-12 * speye(n);

xr = A \ (w(:) .* real(y(:)));
xi = A \ (w(:) .* imag(y(:)));
x = xr + 1i * xi;
end

function save_figure(fig, basePathNoExt)
exportgraphics(fig, [basePathNoExt '.jpg'], 'Resolution', 300, 'BackgroundColor', 'white');
save_visible_fig(fig, [basePathNoExt '.fig']);
end

function save_visible_fig(fig, filename)
prevVisible = get(fig, 'Visible');
set(fig, 'Visible', 'on');
drawnow;
savefig(fig, filename);
set(fig, 'Visible', prevVisible);
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

lgd = legend(ax);
if ~isempty(lgd) && isvalid(lgd)
    lgd.TextColor = 'k';
    lgd.Color = 'w';
end
end

