function bandResult = extract_phase1_band(config, bandName, bandFiles)
%EXTRACT_PHASE1_BAND Run the Phase I extraction for one band.

bandResult = struct();
bandResult.band = bandName;
bandResult.files = bandFiles;
bandResult.outputs = struct();
bandResult.missing = {};

straightKey = 'THRU_P1P2_STRAIGHT';
if ~isfield(bandFiles, straightKey)
    error('Band %s is missing the required straight thru file %s.', bandName, straightKey);
end

straightData = read_touchstone_file(bandFiles.(straightKey));
p34StraightData = [];
if isfield(bandFiles, 'THRU_P3P4_STRAIGHT')
    p34StraightData = read_touchstone_file(bandFiles.THRU_P3P4_STRAIGHT);
end
alphaFit = fit_line_alpha_for_band(config, bandFiles, straightData);
fprintf('  Best alpha for %s GHz band: %.3f Np/m (%.5f dB/um)\n', ...
    bandName, alphaFit.best_alpha_np_per_m, alphaFit.best_alpha_db_per_um);
fprintf('  Best eps_eff for %s GHz band: %.5f\n', bandName, alphaFit.best_eps_eff);
fprintf('  Best Zc for %s GHz band: %.3f ohm\n', bandName, alphaFit.best_zc);
print_alpha_fit_breakdown(alphaFit.best_detail);
lineModel = build_uniform_line_2port(straightData.freq, config.line_length_m, alphaFit.best_eps_eff, alphaFit.best_zc, alphaFit.best_alpha_np_per_m, config.z0);
embeddingExtractOptions = struct( ...
    'jump_threshold', config.embedding_root_jump_threshold, ...
    'blend_floor', config.embedding_root_blend_floor, ...
    'enable_jump_blend', config.embedding_root_enable_jump_blend, ...
    'use_sqrt_candidate', false);
embeddingBaseline = extract_embedding_from_thru(straightData.S, lineModel.S, embeddingExtractOptions);
embedding = embeddingBaseline;
if config.embedding_refine_with_sol
    [embedding, refineDiagnostics] = refine_embedding_with_sol( ...
        embedding, bandFiles, config, straightData, p34StraightData, lineModel, bandName);
else
    refineDiagnostics = [];
end

bandResult.embedding = embedding;
bandResult.embedding_baseline = embeddingBaseline;
bandResult.line_model = lineModel;
bandResult.alpha_fit = alphaFit;
bandResult.straight_p12 = straightData;
bandResult.straight_p34 = p34StraightData;
bandResult.embedding_diagnostics = build_embedding_diagnostics(straightData, p34StraightData, lineModel, config);
bandResult.refine_diagnostics = refineDiagnostics;

embeddingFilename = sprintf('EXTRACTED_EMBEDDING_%s.mat', bandName);
save(fullfile(config.mat_output_dir, embeddingFilename), 'embedding', 'lineModel', 'straightData', 'p34StraightData', 'alphaFit');

solFields = fieldnames(bandFiles);
for idxField = 1:numel(solFields)
    fieldName = solFields{idxField};
    if startsWith(fieldName, 'SHORT_') || startsWith(fieldName, 'OPEN_') || startsWith(fieldName, 'LOAD_')
        solData = read_touchstone_file(bandFiles.(fieldName));
        measuredGamma = select_driven_reflection(solData.S, fieldName);
        parts = split(fieldName, '_');
        portNumber = sscanf(parts{2}, 'P%d');
        localPort = 1;
        if ismember(portNumber, [2, 4])
            localPort = 2;
        end
        gammaDeembeddedRaw = deembed_oneport_standard(measuredGamma, embedding.S, localPort);
        enforcement = [];
        gammaDeembedded = gammaDeembeddedRaw;

        standardName = upper(parts{1});
        portLabel = upper(parts{2});

        outputBase = sprintf('EXTRACTED_%s_%s_%s', standardName, portLabel, bandName);
        write_touchstone_1port(fullfile(config.touchstone_output_dir, [outputBase '.s1p']), solData.freq, gammaDeembedded, config.z0);

        standardResult = struct();
        standardResult.freq = solData.freq;
        standardResult.gamma = gammaDeembedded;
        standardResult.gamma_raw = gammaDeembeddedRaw;
        standardResult.source_file = bandFiles.(fieldName);
        standardResult.measured_gamma = measuredGamma;
        standardResult.passivity_enforcement = enforcement;

        save(fullfile(config.mat_output_dir, [outputBase '.mat']), 'standardResult');
        bandResult.outputs.(matlab.lang.makeValidName(outputBase)) = standardResult;
    end
end

for idxField = 1:numel(solFields)
    fieldName = solFields{idxField};
    if startsWith(fieldName, 'THRU_')
        thruData = read_touchstone_file(bandFiles.(fieldName));
        sDeembedded = deembed_twport_standard(thruData.S, embedding.S);

        outputBase = sprintf('EXTRACTED_%s_%s', fieldName, bandName);
        write_touchstone_2port(fullfile(config.touchstone_output_dir, [outputBase '.s2p']), thruData.freq, sDeembedded, config.z0);

        reciprocalResult = struct();
        reciprocalResult.freq = thruData.freq;
        reciprocalResult.S = sDeembedded;
        reciprocalResult.source_file = bandFiles.(fieldName);
        reciprocalResult.measured_S = thruData.S;

        save(fullfile(config.mat_output_dir, [outputBase '.mat']), 'reciprocalResult');
        bandResult.outputs.(matlab.lang.makeValidName(outputBase)) = reciprocalResult;
    end
end

plot_phase1_band_results(config, bandResult);

bandSummary = struct();
bandSummary.band = bandName;
bandSummary.available_files = fieldnames(bandFiles);
bandSummary.touchstone_output_dir = config.touchstone_output_dir;
bandSummary.mat_output_dir = config.mat_output_dir;
bandSummary.alpha_fit = alphaFit;
bandSummary.embedding_diagnostics = bandResult.embedding_diagnostics.summary;
bandSummary.refine_diagnostics = summarize_refine_diagnostics(refineDiagnostics);
save(fullfile(config.mat_output_dir, sprintf('PHASE1_SUMMARY_%s.mat', bandName)), 'bandSummary');

write_alpha_fit_breakdown(config, bandName, alphaFit);
write_embedding_diagnostic_notes(config, bandName, bandResult.embedding_diagnostics);
write_refinement_notes(config, bandName, refineDiagnostics);
end

function gammaMeasured = select_driven_reflection(sParams, fieldName)
parts = split(fieldName, '_');
portNumber = sscanf(parts{2}, 'P%d');

if ismember(portNumber, [1, 3])
    gammaMeasured = squeeze(sParams(1, 1, :));
elseif ismember(portNumber, [2, 4])
    gammaMeasured = squeeze(sParams(2, 2, :));
else
    error('Unsupported port label in field name: %s', fieldName);
end
end

function alphaFit = fit_line_alpha_for_band(config, bandFiles, straightData)
p12Data = straightData;

p34Data = [];
if isfield(bandFiles, 'THRU_P3P4_STRAIGHT')
    p34Data = read_touchstone_file(bandFiles.THRU_P3P4_STRAIGHT);
end

shortData = [];
if isfield(bandFiles, 'SHORT_P1')
    shortData = read_touchstone_file(bandFiles.SHORT_P1);
end

openData = [];
if isfield(bandFiles, 'OPEN_P1')
    openData = read_touchstone_file(bandFiles.OPEN_P1);
end

loadData = [];
if isfield(bandFiles, 'LOAD_P1')
    loadData = read_touchstone_file(bandFiles.LOAD_P1);
end

weights = config.alpha_fit_weights;
coarse = evaluate_search_grid( ...
    config.alpha_search_np_per_m, ...
    config.eps_eff_search, ...
    config.zc_search, ...
    config, weights, p12Data, p34Data, shortData, openData, loadData);

bestCoarse = coarse.best_detail;
alphaRefine = linspace( ...
    max(0, bestCoarse.alpha_np_per_m - config.alpha_refine_halfspan_np_per_m), ...
    bestCoarse.alpha_np_per_m + config.alpha_refine_halfspan_np_per_m, ...
    config.alpha_refine_points);
epsRefine = linspace( ...
    bestCoarse.eps_eff - config.eps_eff_refine_halfspan, ...
    bestCoarse.eps_eff + config.eps_eff_refine_halfspan, ...
    config.eps_eff_refine_points);
zcRefine = linspace( ...
    bestCoarse.zc - config.zc_refine_halfspan, ...
    bestCoarse.zc + config.zc_refine_halfspan, ...
    config.zc_refine_points);

refined = evaluate_search_grid( ...
    alphaRefine, ...
    epsRefine, ...
    zcRefine, ...
    config, weights, p12Data, p34Data, shortData, openData, loadData);

alphaFit = struct();
alphaFit.best_alpha_np_per_m = refined.best_detail.alpha_np_per_m;
alphaFit.best_alpha_db_per_um = alphaFit.best_alpha_np_per_m * 8.686 / 1e6;
alphaFit.best_eps_eff = refined.best_detail.eps_eff;
alphaFit.best_zc = refined.best_detail.zc;
alphaFit.alpha_grid_np_per_m = coarse.alpha_grid_np_per_m;
alphaFit.eps_eff_grid = coarse.eps_eff_grid;
alphaFit.zc_grid = coarse.zc_grid;
alphaFit.total_scores = coarse.total_scores;
alphaFit.details = coarse.details;
alphaFit.best_detail = refined.best_detail;
alphaFit.coarse = coarse;
alphaFit.refined = refined;
end

function penalty = compute_embedding_passivity_penalty(sEmbedding)
numFreq = size(sEmbedding, 3);
excessTerms = zeros(numFreq, 1);
for idx = 1:numFreq
    sigmaMax = max(svd(sEmbedding(:, :, idx)));
    excessTerms(idx) = max(sigmaMax - 1, 0);
end
penalty = mean(excessTerms .^ 2, 'omitnan') + max(excessTerms, [], 'omitnan')^2;
end

function searchResult = evaluate_search_grid(alphaGrid, epsGrid, zcGrid, config, weights, p12Data, p34Data, shortData, openData, loadData)
numAlpha = numel(alphaGrid);
numEps = numel(epsGrid);
numZc = numel(zcGrid);
numCandidates = numAlpha * numEps * numZc;

candidateList = zeros(numCandidates, 3);
idx = 1;
for iAlpha = 1:numAlpha
    for iEps = 1:numEps
        for iZc = 1:numZc
            candidateList(idx, :) = [alphaGrid(iAlpha), epsGrid(iEps), zcGrid(iZc)];
            idx = idx + 1;
        end
    end
end

scoreList = zeros(numCandidates, 1);
detailList = repmat(empty_detail_struct(), numCandidates, 1);

useParallel = config.use_parallel && license('test', 'Distrib_Computing_Toolbox');
if useParallel
    parfor idxCandidate = 1:numCandidates
        [scoreList(idxCandidate), detailList(idxCandidate)] = evaluate_candidate(candidateList(idxCandidate, :), config, weights, p12Data, p34Data, shortData, openData, loadData);
    end
else
    for idxCandidate = 1:numCandidates
        [scoreList(idxCandidate), detailList(idxCandidate)] = evaluate_candidate(candidateList(idxCandidate, :), config, weights, p12Data, p34Data, shortData, openData, loadData);
    end
end

totalScores = reshape(scoreList, [numAlpha, numEps, numZc]);
details = reshape(detailList, [numAlpha, numEps, numZc]);

[~, bestLinearIdx] = min(scoreList);

searchResult = struct();
searchResult.alpha_grid_np_per_m = alphaGrid;
searchResult.eps_eff_grid = epsGrid;
searchResult.zc_grid = zcGrid;
searchResult.total_scores = totalScores;
searchResult.details = details;
searchResult.best_detail = detailList(bestLinearIdx);
end

function [totalScore, detail] = evaluate_candidate(candidate, config, weights, p12Data, p34Data, shortData, openData, loadData)
alphaCandidate = candidate(1);
epsCandidate = candidate(2);
zcCandidate = candidate(3);

lineCandidate = build_uniform_line_2port(p12Data.freq, config.line_length_m, epsCandidate, zcCandidate, alphaCandidate, config.z0);
            embeddingCandidate = extract_embedding_from_thru(p12Data.S, lineCandidate.S);

passivityPenalty = 0;
targetPenalty = 0;
if ~isempty(shortData)
    gammaShort = deembed_oneport_standard(select_driven_reflection(shortData.S, 'SHORT_P1'), embeddingCandidate.S);
    passivityPenalty = passivityPenalty + mean(max(abs(gammaShort) - 1, 0).^2, 'omitnan');
    targetPenalty = targetPenalty + mean((abs(gammaShort) - 1).^2, 'omitnan');
end
if ~isempty(openData)
    gammaOpen = deembed_oneport_standard(select_driven_reflection(openData.S, 'OPEN_P1'), embeddingCandidate.S);
    passivityPenalty = passivityPenalty + mean(max(abs(gammaOpen) - 1, 0).^2, 'omitnan');
    targetPenalty = targetPenalty + mean((abs(gammaOpen) - 1).^2, 'omitnan');
end

loadPenalty = 0;
if ~isempty(loadData)
    gammaLoad = deembed_oneport_standard(select_driven_reflection(loadData.S, 'LOAD_P1'), embeddingCandidate.S);
    loadPenalty = mean(abs(gammaLoad).^2, 'omitnan');
end

straightLineMatchPenalty = 0;
straightReturnPenalty = 0;
embeddingPassivityPenalty = compute_embedding_passivity_penalty(embeddingCandidate.S);
lineS21 = squeeze(lineCandidate.S(2, 1, :));

if ~isempty(p34Data)
    p34Deembedded = deembed_twport_standard(p34Data.S, embeddingCandidate.S);
    p34S21 = squeeze(p34Deembedded(2, 1, :));
    p34S11 = squeeze(p34Deembedded(1, 1, :));
    p34S22 = squeeze(p34Deembedded(2, 2, :));
    straightLineMatchPenalty = mean(abs(p34S21 - lineS21).^2, 'omitnan');
    straightReturnPenalty = mean(abs(p34S11).^2 + abs(p34S22).^2, 'omitnan');
end

alphaRegularization = ((alphaCandidate - config.alpha_regularization_center_np_per_m) / config.alpha_regularization_scale_np_per_m)^2;
epsRegularization = ((epsCandidate - config.eps_eff_regularization_center) / config.eps_eff_regularization_scale)^2;
zcRegularization = ((zcCandidate - config.zc_regularization_center) / config.zc_regularization_scale)^2;

totalScore = ...
    weights.open_short_passivity * passivityPenalty + ...
    weights.open_short_target * targetPenalty + ...
    weights.load_magnitude * loadPenalty + ...
    weights.straight_line_match * straightLineMatchPenalty + ...
    weights.straight_return_loss * straightReturnPenalty + ...
    weights.embedding_passivity * embeddingPassivityPenalty + ...
    weights.alpha_regularization * alphaRegularization + ...
    weights.eps_regularization * epsRegularization + ...
    weights.zc_regularization * zcRegularization;

detail = empty_detail_struct();
detail.alpha_np_per_m = alphaCandidate;
detail.alpha_db_per_um = alphaCandidate * 8.686 / 1e6;
detail.eps_eff = epsCandidate;
detail.zc = zcCandidate;
detail.open_short_passivity = passivityPenalty;
detail.open_short_target = targetPenalty;
detail.load_magnitude = loadPenalty;
detail.straight_line_match = straightLineMatchPenalty;
detail.straight_return_loss = straightReturnPenalty;
detail.embedding_passivity = embeddingPassivityPenalty;
detail.alpha_regularization = alphaRegularization;
detail.eps_regularization = epsRegularization;
detail.zc_regularization = zcRegularization;
detail.total_score = totalScore;
end

function detail = empty_detail_struct()
detail = struct( ...
    'alpha_np_per_m', 0, ...
    'alpha_db_per_um', 0, ...
    'eps_eff', 0, ...
    'zc', 0, ...
    'open_short_passivity', NaN, ...
    'open_short_target', NaN, ...
    'load_magnitude', NaN, ...
    'straight_line_match', NaN, ...
    'straight_return_loss', NaN, ...
    'embedding_passivity', NaN, ...
    'alpha_regularization', NaN, ...
    'eps_regularization', NaN, ...
    'zc_regularization', NaN, ...
    'total_score', NaN);
end

function print_alpha_fit_breakdown(bestDetail)
fprintf('  Objective breakdown:\n');
fprintf('    open_short_passivity: %.6g\n', bestDetail.open_short_passivity);
fprintf('    open_short_target: %.6g\n', bestDetail.open_short_target);
fprintf('    load_magnitude: %.6g\n', bestDetail.load_magnitude);
fprintf('    straight_line_match: %.6g\n', bestDetail.straight_line_match);
fprintf('    straight_return_loss: %.6g\n', bestDetail.straight_return_loss);
fprintf('    embedding_passivity: %.6g\n', bestDetail.embedding_passivity);
fprintf('    alpha_regularization: %.6g\n', bestDetail.alpha_regularization);
fprintf('    eps_regularization: %.6g\n', bestDetail.eps_regularization);
fprintf('    zc_regularization: %.6g\n', bestDetail.zc_regularization);
fprintf('    total_score: %.6g\n', bestDetail.total_score);
end

function write_alpha_fit_breakdown(config, bandName, alphaFit)
notesDir = fullfile(fileparts(config.touchstone_output_dir), 'notes');
if ~exist(notesDir, 'dir')
    mkdir(notesDir);
end

filename = fullfile(notesDir, sprintf('PhaseI_ObjectiveBreakdown_%s.txt', bandName));
fid = fopen(filename, 'w');
if fid == -1
    warning('Unable to write objective breakdown note: %s', filename);
    return;
end

cleanup = onCleanup(@() fclose(fid));
best = alphaFit.best_detail;

fprintf(fid, 'Phase I Objective Breakdown - %s GHz band\n\n', bandName);
fprintf(fid, 'Best alpha: %.9g Np/m\n', alphaFit.best_alpha_np_per_m);
fprintf(fid, 'Best alpha: %.9g dB/um\n', alphaFit.best_alpha_db_per_um);
fprintf(fid, 'Best eps_eff: %.9g\n', alphaFit.best_eps_eff);
fprintf(fid, 'Best Zc: %.9g ohm\n\n', alphaFit.best_zc);
fprintf(fid, 'Objective terms at optimum:\n');
fprintf(fid, 'open_short_passivity = %.9g\n', best.open_short_passivity);
fprintf(fid, 'open_short_target = %.9g\n', best.open_short_target);
fprintf(fid, 'load_magnitude = %.9g\n', best.load_magnitude);
fprintf(fid, 'straight_line_match = %.9g\n', best.straight_line_match);
fprintf(fid, 'straight_return_loss = %.9g\n', best.straight_return_loss);
fprintf(fid, 'embedding_passivity = %.9g\n', best.embedding_passivity);
fprintf(fid, 'alpha_regularization = %.9g\n', best.alpha_regularization);
fprintf(fid, 'eps_regularization = %.9g\n', best.eps_regularization);
fprintf(fid, 'zc_regularization = %.9g\n', best.zc_regularization);
fprintf(fid, 'total_score = %.9g\n', best.total_score);
end

function diagnostics = build_embedding_diagnostics(straightData, p34StraightData, lineModel, config)
embeddingExtractOptions = struct( ...
    'jump_threshold', config.embedding_root_jump_threshold, ...
    'blend_floor', config.embedding_root_blend_floor, ...
    'enable_jump_blend', config.embedding_root_enable_jump_blend, ...
    'use_sqrt_candidate', false);
embedding12 = extract_embedding_from_thru(straightData.S, lineModel.S, embeddingExtractOptions);

diagnostics = struct();
diagnostics.freq = straightData.freq;
diagnostics.embedding12 = embedding12;
diagnostics.embedding34 = [];
diagnostics.summary = struct();

diagnostics.summary.has_p34 = ~isempty(p34StraightData);
diagnostics.summary.mean_s11_diff = NaN;
diagnostics.summary.mean_s21_diff = NaN;
diagnostics.summary.mean_s22_diff = NaN;
diagnostics.summary.mean_sigma_diff = NaN;

if isempty(p34StraightData)
    return;
end

embedding34 = extract_embedding_from_thru(p34StraightData.S, lineModel.S, embeddingExtractOptions);
diagnostics.embedding34 = embedding34;

s11Diff = abs(squeeze(embedding12.S(1, 1, :)) - squeeze(embedding34.S(1, 1, :)));
s21Diff = abs(squeeze(embedding12.S(2, 1, :)) - squeeze(embedding34.S(2, 1, :)));
s22Diff = abs(squeeze(embedding12.S(2, 2, :)) - squeeze(embedding34.S(2, 2, :)));

sigma12 = compute_sigma_max_vector(embedding12.S);
sigma34 = compute_sigma_max_vector(embedding34.S);
sigmaDiff = abs(sigma12 - sigma34);

diagnostics.summary.mean_s11_diff = mean(s11Diff, 'omitnan');
diagnostics.summary.mean_s21_diff = mean(s21Diff, 'omitnan');
diagnostics.summary.mean_s22_diff = mean(s22Diff, 'omitnan');
diagnostics.summary.mean_sigma_diff = mean(sigmaDiff, 'omitnan');
end

function sigmaMax = compute_sigma_max_vector(sParams)
numFreq = size(sParams, 3);
sigmaMax = zeros(numFreq, 1);
for idx = 1:numFreq
    sigmaMax(idx) = max(svd(sParams(:, :, idx)));
end
end

function write_embedding_diagnostic_notes(config, bandName, diagnostics)
notesDir = fullfile(fileparts(config.touchstone_output_dir), 'notes');
if ~exist(notesDir, 'dir')
    mkdir(notesDir);
end

filename = fullfile(notesDir, sprintf('PhaseI_EmbeddingComparison_%s.txt', bandName));
fid = fopen(filename, 'w');
if fid == -1
    warning('Unable to write embedding diagnostic note: %s', filename);
    return;
end

cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Phase I Embedding Comparison - %s GHz band\n\n', bandName);
fprintf(fid, 'Has P3P4 straight-thru comparison: %d\n', diagnostics.summary.has_p34);
fprintf(fid, 'Mean |S11(E12) - S11(E34)| = %.9g\n', diagnostics.summary.mean_s11_diff);
fprintf(fid, 'Mean |S21(E12) - S21(E34)| = %.9g\n', diagnostics.summary.mean_s21_diff);
fprintf(fid, 'Mean |S22(E12) - S22(E34)| = %.9g\n', diagnostics.summary.mean_s22_diff);
fprintf(fid, 'Mean |sigma_max(E12) - sigma_max(E34)| = %.9g\n', diagnostics.summary.mean_sigma_diff);
end

function summary = summarize_refine_diagnostics(refineDiagnostics)
if isempty(refineDiagnostics)
    summary = [];
    return;
end

summary = struct();
summary.mean_s11_shift = mean(refineDiagnostics.abs_s11_shift, 'omitnan');
summary.mean_s12_shift = mean(refineDiagnostics.abs_s12_shift, 'omitnan');
summary.mean_s21_shift = mean(refineDiagnostics.abs_s21_shift, 'omitnan');
summary.mean_s22_shift = mean(refineDiagnostics.abs_s22_shift, 'omitnan');
summary.mean_short_improvement = mean(refineDiagnostics.short_mag_improvement, 'omitnan');
summary.mean_open_improvement = mean(refineDiagnostics.open_mag_improvement, 'omitnan');
summary.mean_load_improvement = mean(refineDiagnostics.load_mag_improvement, 'omitnan');
summary.mean_sigma_change = mean(refineDiagnostics.sigma_change, 'omitnan');
summary.mean_thru_reconstruction_error = mean(refineDiagnostics.thru_reconstruction_error, 'omitnan');
summary.mean_p34_line_match_error = mean(refineDiagnostics.p34_line_match_error, 'omitnan');
summary.mean_p34_return_loss_error = mean(refineDiagnostics.p34_return_loss_error, 'omitnan');
summary.mean_simple_short_discrepancy = mean(refineDiagnostics.simple_short_discrepancy, 'omitnan');
summary.mean_simple_open_discrepancy = mean(refineDiagnostics.simple_open_discrepancy, 'omitnan');
summary.mean_simple_load_discrepancy = mean(refineDiagnostics.simple_load_discrepancy, 'omitnan');
if isfield(refineDiagnostics, 'sim_ref_error')
    summary.mean_sim_ref_error = mean(refineDiagnostics.sim_ref_error, 'omitnan');
    summary.mean_sim_raw_error = mean(refineDiagnostics.sim_raw_error, 'omitnan');
end
end

function write_refinement_notes(config, bandName, refineDiagnostics)
if isempty(refineDiagnostics)
    return;
end

notesDir = fullfile(fileparts(config.touchstone_output_dir), 'notes');
if ~exist(notesDir, 'dir')
    mkdir(notesDir);
end

filename = fullfile(notesDir, sprintf('PhaseI_RefinementSummary_%s.txt', bandName));
fid = fopen(filename, 'w');
if fid == -1
    warning('Unable to write refinement note: %s', filename);
    return;
end

cleanup = onCleanup(@() fclose(fid));

fprintf(fid, 'Phase I Reflective Refinement Summary - %s GHz band\n\n', bandName);
fprintf(fid, 'Mean |delta S11| = %.9g\n', mean(refineDiagnostics.abs_s11_shift, 'omitnan'));
fprintf(fid, 'Mean |delta S12| = %.9g\n', mean(refineDiagnostics.abs_s12_shift, 'omitnan'));
fprintf(fid, 'Mean |delta S21| = %.9g\n', mean(refineDiagnostics.abs_s21_shift, 'omitnan'));
fprintf(fid, 'Mean |delta S22| = %.9g\n', mean(refineDiagnostics.abs_s22_shift, 'omitnan'));
fprintf(fid, 'Mean short |Gamma|-1 reduction = %.9g\n', mean(refineDiagnostics.short_mag_improvement, 'omitnan'));
fprintf(fid, 'Mean open |Gamma|-1 reduction = %.9g\n', mean(refineDiagnostics.open_mag_improvement, 'omitnan'));
fprintf(fid, 'Mean load |Gamma| reduction = %.9g\n', mean(refineDiagnostics.load_mag_improvement, 'omitnan'));
fprintf(fid, 'Mean sigma_max change = %.9g\n', mean(refineDiagnostics.sigma_change, 'omitnan'));
fprintf(fid, 'Mean straight-thru reconstruction error = %.9g\n', mean(refineDiagnostics.thru_reconstruction_error, 'omitnan'));
fprintf(fid, 'Mean P3P4 line-match error = %.9g\n', mean(refineDiagnostics.p34_line_match_error, 'omitnan'));
fprintf(fid, 'Mean P3P4 return-loss error = %.9g\n', mean(refineDiagnostics.p34_return_loss_error, 'omitnan'));
fprintf(fid, 'Mean simple short discrepancy = %.9g\n', mean(refineDiagnostics.simple_short_discrepancy, 'omitnan'));
fprintf(fid, 'Mean simple open discrepancy = %.9g\n', mean(refineDiagnostics.simple_open_discrepancy, 'omitnan'));
fprintf(fid, 'Mean simple load discrepancy = %.9g\n', mean(refineDiagnostics.simple_load_discrepancy, 'omitnan'));
if isfield(refineDiagnostics, 'sim_ref_error')
    fprintf(fid, 'Mean simulation ref error = %.9g\n', mean(refineDiagnostics.sim_ref_error, 'omitnan'));
    fprintf(fid, 'Mean simulation raw error = %.9g\n', mean(refineDiagnostics.sim_raw_error, 'omitnan'));
end
end

