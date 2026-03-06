function [embeddingOut, diagnostics] = refine_embedding_with_sol(embeddingIn, bandFiles, config, straightData, p34StraightData, lineModel)
%REFINE_EMBEDDING_WITH_SOL Weighted Stage A embedding refinement.
% Refines only S11/S22 per frequency while keeping S12/S21 fixed. Adds
% stronger weighting near band edges and overlap windows.

requiredFields = {'SHORT_P1', 'OPEN_P1', 'LOAD_P1', 'SHORT_P2', 'OPEN_P2', 'LOAD_P2'};
for idx = 1:numel(requiredFields)
    if ~isfield(bandFiles, requiredFields{idx})
        embeddingOut = embeddingIn;
        diagnostics = [];
        return;
    end
end

shortData = read_touchstone_file(bandFiles.SHORT_P1);
openData = read_touchstone_file(bandFiles.OPEN_P1);
loadData = read_touchstone_file(bandFiles.LOAD_P1);
shortDataP2 = read_touchstone_file(bandFiles.SHORT_P2);
openDataP2 = read_touchstone_file(bandFiles.OPEN_P2);
loadDataP2 = read_touchstone_file(bandFiles.LOAD_P2);

gammaShortMeasuredP1 = squeeze(shortData.S(1, 1, :));
gammaOpenMeasuredP1 = squeeze(openData.S(1, 1, :));
gammaLoadMeasuredP1 = squeeze(loadData.S(1, 1, :));
gammaShortMeasuredP2 = squeeze(shortDataP2.S(2, 2, :));
gammaOpenMeasuredP2 = squeeze(openDataP2.S(2, 2, :));
gammaLoadMeasuredP2 = squeeze(loadDataP2.S(2, 2, :));

sParams = embeddingIn.S;
numFreq = size(sParams, 3);
weights = config.embedding_refine_weights;
diagnostics = initialize_diagnostics(numFreq);

freqVector = straightData.freq(:);
freqRange = [freqVector(1), freqVector(end)];

options = optimset( ...
    'Display', 'off', ...
    'MaxIter', config.embedding_refine_max_iter, ...
    'MaxFunEvals', config.embedding_refine_max_fun_evals, ...
    'TolX', 1e-6, ...
    'TolFun', 1e-6);

for idx = 1:numFreq
    s11Base = sParams(1, 1, idx);
    s12 = sParams(1, 2, idx);
    s21 = sParams(2, 1, idx);
    s22Base = sParams(2, 2, idx);

    if idx == 1
        x0 = [real(s11Base), imag(s11Base), real(s22Base), imag(s22Base)];
        prevS11 = [];
        prevS22 = [];
    else
        x0 = [real(sParams(1, 1, idx - 1)), imag(sParams(1, 1, idx - 1)), real(sParams(2, 2, idx - 1)), imag(sParams(2, 2, idx - 1))];
        prevS11 = sParams(1, 1, idx - 1);
        prevS22 = sParams(2, 2, idx - 1);
    end

    localWeight = frequency_weight(freqVector(idx), freqRange, config);
    objective = @(x) local_refinement_objective( ...
        x, ...
        s11Base, s12, s21, s22Base, ...
        prevS11, prevS22, ...
        gammaShortMeasuredP1(idx), gammaOpenMeasuredP1(idx), gammaLoadMeasuredP1(idx), ...
        gammaShortMeasuredP2(idx), gammaOpenMeasuredP2(idx), gammaLoadMeasuredP2(idx), ...
        squeeze(straightData.S(:, :, idx)), ...
        squeeze(lineModel.S(:, :, idx)), ...
        get_optional_slice(p34StraightData, idx), ...
        weights, localWeight, config.z0);

    xOpt = fminsearch(objective, x0, options);
    s11Refined = xOpt(1) + 1i * xOpt(2);
    s22Refined = xOpt(3) + 1i * xOpt(4);

    sParams(1, 1, idx) = s11Refined;
    sParams(2, 2, idx) = s22Refined;

    diagnostics = update_diagnostics( ...
        diagnostics, idx, ...
        s11Base, s12, s21, s22Base, ...
        s11Refined, s22Refined, ...
        gammaShortMeasuredP1(idx), gammaOpenMeasuredP1(idx), gammaLoadMeasuredP1(idx), ...
        gammaShortMeasuredP2(idx), gammaOpenMeasuredP2(idx), gammaLoadMeasuredP2(idx), ...
        squeeze(straightData.S(:, :, idx)), ...
        squeeze(lineModel.S(:, :, idx)), ...
        get_optional_slice(p34StraightData, idx), ...
        config.z0);
end

tParams = zeros(2, 2, numFreq);
for idx = 1:numFreq
    tParams(:, :, idx) = s_to_abcd_local(sParams(:, :, idx));
end

embeddingOut = struct();
embeddingOut.S = sParams;
embeddingOut.T = tParams;
end

function diagnostics = initialize_diagnostics(numFreq)
diagnostics = struct();
diagnostics.abs_s11_shift = zeros(numFreq, 1);
diagnostics.abs_s12_shift = zeros(numFreq, 1);
diagnostics.abs_s21_shift = zeros(numFreq, 1);
diagnostics.abs_s22_shift = zeros(numFreq, 1);
diagnostics.short_mag_improvement = zeros(numFreq, 1);
diagnostics.open_mag_improvement = zeros(numFreq, 1);
diagnostics.load_mag_improvement = zeros(numFreq, 1);
diagnostics.sigma_change = zeros(numFreq, 1);
diagnostics.thru_reconstruction_error = zeros(numFreq, 1);
diagnostics.p34_line_match_error = NaN(numFreq, 1);
diagnostics.p34_return_loss_error = NaN(numFreq, 1);
diagnostics.simple_short_discrepancy = zeros(numFreq, 1);
diagnostics.simple_open_discrepancy = zeros(numFreq, 1);
diagnostics.simple_load_discrepancy = zeros(numFreq, 1);
end

function w = frequency_weight(freqHz, freqRange, config)
fLo = freqRange(1);
fHi = freqRange(2);
span = max(config.embedding_refine_edge_span_hz, 1);
distEdge = min(freqHz - fLo, fHi - freqHz);
edgeTaper = max(0, 1 - distEdge / span);
edgeWeight = config.embedding_refine_edge_boost * edgeTaper;

overlapWeight = 0;
windows = config.embedding_refine_overlap_windows_hz;
for idx = 1:size(windows, 1)
    if freqHz >= windows(idx, 1) && freqHz <= windows(idx, 2)
        overlapWeight = overlapWeight + config.embedding_refine_overlap_boost;
    end
end
w = 1 + edgeWeight + overlapWeight;
end

function score = local_refinement_objective(x, s11Base, s12, s21, s22Base, prevS11, prevS22, gammaShortMeasuredP1, gammaOpenMeasuredP1, gammaLoadMeasuredP1, gammaShortMeasuredP2, gammaOpenMeasuredP2, gammaLoadMeasuredP2, p12Measured, lineS, p34Measured, weights, localWeight, z0)
s11 = x(1) + 1i * x(2);
s22 = x(3) + 1i * x(4);

gammaShortP1 = oneport_deembed_formula_port1(gammaShortMeasuredP1, s11, s12, s21, s22);
gammaOpenP1 = oneport_deembed_formula_port1(gammaOpenMeasuredP1, s11, s12, s21, s22);
gammaLoadP1 = oneport_deembed_formula_port1(gammaLoadMeasuredP1, s11, s12, s21, s22);

gammaShortP2 = oneport_deembed_formula_port2(gammaShortMeasuredP2, s11, s12, s21, s22);
gammaOpenP2 = oneport_deembed_formula_port2(gammaOpenMeasuredP2, s11, s12, s21, s22);
gammaLoadP2 = oneport_deembed_formula_port2(gammaLoadMeasuredP2, s11, s12, s21, s22);

openShortPassivityP1 = max(abs(gammaShortP1) - 1, 0)^2 + max(abs(gammaOpenP1) - 1, 0)^2;
openShortTargetP1 = (abs(gammaShortP1) - 1)^2 + (abs(gammaOpenP1) - 1)^2;
loadMagnitudeP1 = abs(gammaLoadP1)^2;

openShortPassivityP2 = max(abs(gammaShortP2) - 1, 0)^2 + max(abs(gammaOpenP2) - 1, 0)^2;
openShortTargetP2 = (abs(gammaShortP2) - 1)^2 + (abs(gammaOpenP2) - 1)^2;
loadMagnitudeP2 = abs(gammaLoadP2)^2;

baselinePenalty = ...
    weights.baseline_s11 * abs(s11 - s11Base)^2 + ...
    weights.baseline_s22 * abs(s22 - s22Base)^2;

smoothPenalty = 0;
if ~isempty(prevS11)
    smoothPenalty = ...
        weights.smoothness_s11 * abs(s11 - prevS11)^2 + ...
        weights.smoothness_s22 * abs(s22 - prevS22)^2;
end

sMatrix = [s11, s12; s21, s22];
sigmaMax = max(svd(sMatrix));
sigmaExcess = max(sigmaMax - 1, 0);
embeddingPassivity = ...
    weights.embedding_passivity_mean * sigmaExcess^2 + ...
    weights.embedding_passivity_local * sigmaExcess^4;

thruRecon = reconstruct_thru_from_embedding(sMatrix, lineS, z0);
thruPenalty = norm(thruRecon - p12Measured, 'fro')^2;

p34Penalty = 0;
if ~isempty(p34Measured)
    p34Actual = deembed_twport_standard_single(p34Measured, sMatrix, z0);
    p34Penalty = ...
        weights.p34_line_match * norm(p34Actual - lineS, 'fro')^2 + ...
        weights.p34_return_loss * (abs(p34Actual(1, 1))^2 + abs(p34Actual(2, 2))^2);
end

score = ...
    localWeight * weights.open_short_passivity_p1 * openShortPassivityP1 + ...
    localWeight * weights.open_short_target_p1 * openShortTargetP1 + ...
    localWeight * weights.load_magnitude_p1 * loadMagnitudeP1 + ...
    localWeight * weights.open_short_passivity_p2 * openShortPassivityP2 + ...
    localWeight * weights.open_short_target_p2 * openShortTargetP2 + ...
    localWeight * weights.load_magnitude_p2 * loadMagnitudeP2 + ...
    localWeight * weights.thru_reconstruct * thruPenalty + ...
    localWeight * p34Penalty + ...
    baselinePenalty + ...
    smoothPenalty + ...
    localWeight * embeddingPassivity;
end

function diagnostics = update_diagnostics(diagnostics, idx, s11Base, s12, s21, s22Base, s11Refined, s22Refined, gammaShortMeasuredP1, gammaOpenMeasuredP1, gammaLoadMeasuredP1, gammaShortMeasuredP2, gammaOpenMeasuredP2, gammaLoadMeasuredP2, p12Measured, lineS, p34Measured, z0)
gammaShortBaseP1 = oneport_deembed_formula_port1(gammaShortMeasuredP1, s11Base, s12, s21, s22Base);
gammaOpenBaseP1 = oneport_deembed_formula_port1(gammaOpenMeasuredP1, s11Base, s12, s21, s22Base);
gammaLoadBaseP1 = oneport_deembed_formula_port1(gammaLoadMeasuredP1, s11Base, s12, s21, s22Base);
gammaShortBaseP2 = oneport_deembed_formula_port2(gammaShortMeasuredP2, s11Base, s12, s21, s22Base);
gammaOpenBaseP2 = oneport_deembed_formula_port2(gammaOpenMeasuredP2, s11Base, s12, s21, s22Base);
gammaLoadBaseP2 = oneport_deembed_formula_port2(gammaLoadMeasuredP2, s11Base, s12, s21, s22Base);

gammaShortRefinedP1 = oneport_deembed_formula_port1(gammaShortMeasuredP1, s11Refined, s12, s21, s22Refined);
gammaOpenRefinedP1 = oneport_deembed_formula_port1(gammaOpenMeasuredP1, s11Refined, s12, s21, s22Refined);
gammaLoadRefinedP1 = oneport_deembed_formula_port1(gammaLoadMeasuredP1, s11Refined, s12, s21, s22Refined);
gammaShortRefinedP2 = oneport_deembed_formula_port2(gammaShortMeasuredP2, s11Refined, s12, s21, s22Refined);
gammaOpenRefinedP2 = oneport_deembed_formula_port2(gammaOpenMeasuredP2, s11Refined, s12, s21, s22Refined);
gammaLoadRefinedP2 = oneport_deembed_formula_port2(gammaLoadMeasuredP2, s11Refined, s12, s21, s22Refined);

gammaShortSimpleP1 = oneport_simple_transmission_only(gammaShortMeasuredP1, s12, s21);
gammaOpenSimpleP1 = oneport_simple_transmission_only(gammaOpenMeasuredP1, s12, s21);
gammaLoadSimpleP1 = oneport_simple_transmission_only(gammaLoadMeasuredP1, s12, s21);

sigmaBase = max(svd([s11Base, s12; s21, s22Base]));
sigmaRefined = max(svd([s11Refined, s12; s21, s22Refined]));

sEmbedRefined = [s11Refined, s12; s21, s22Refined];
sThruRecon = reconstruct_thru_from_embedding(sEmbedRefined, lineS, z0);

diagnostics.abs_s11_shift(idx) = abs(s11Refined - s11Base);
diagnostics.abs_s12_shift(idx) = 0;
diagnostics.abs_s21_shift(idx) = 0;
diagnostics.abs_s22_shift(idx) = abs(s22Refined - s22Base);
diagnostics.short_mag_improvement(idx) = 0.5 * ( ...
    abs(abs(gammaShortBaseP1) - 1) - abs(abs(gammaShortRefinedP1) - 1) + ...
    abs(abs(gammaShortBaseP2) - 1) - abs(abs(gammaShortRefinedP2) - 1));
diagnostics.open_mag_improvement(idx) = 0.5 * ( ...
    abs(abs(gammaOpenBaseP1) - 1) - abs(abs(gammaOpenRefinedP1) - 1) + ...
    abs(abs(gammaOpenBaseP2) - 1) - abs(abs(gammaOpenRefinedP2) - 1));
diagnostics.load_mag_improvement(idx) = 0.5 * ( ...
    abs(gammaLoadBaseP1) - abs(gammaLoadRefinedP1) + ...
    abs(gammaLoadBaseP2) - abs(gammaLoadRefinedP2));
diagnostics.sigma_change(idx) = sigmaRefined - sigmaBase;
diagnostics.thru_reconstruction_error(idx) = norm(sThruRecon - p12Measured, 'fro')^2;

if ~isempty(p34Measured)
    p34Actual = deembed_twport_standard_single(p34Measured, sEmbedRefined, z0);
    diagnostics.p34_line_match_error(idx) = norm(p34Actual - lineS, 'fro')^2;
    diagnostics.p34_return_loss_error(idx) = abs(p34Actual(1, 1))^2 + abs(p34Actual(2, 2))^2;
end

diagnostics.simple_short_discrepancy(idx) = abs(gammaShortRefinedP1 - gammaShortSimpleP1);
diagnostics.simple_open_discrepancy(idx) = abs(gammaOpenRefinedP1 - gammaOpenSimpleP1);
diagnostics.simple_load_discrepancy(idx) = abs(gammaLoadRefinedP1 - gammaLoadSimpleP1);
end

function gammaActual = oneport_deembed_formula_port1(gammaMeasured, s11, s12, s21, s22)
deltaGamma = gammaMeasured - s11;
denominator = (s12 * s21) + (s22 * deltaGamma);
if abs(denominator) < eps
    gammaActual = NaN;
else
    gammaActual = deltaGamma / denominator;
end
end

function gammaActual = oneport_deembed_formula_port2(gammaMeasured, s11, s12, s21, s22)
deltaGamma = gammaMeasured - s22;
denominator = (s12 * s21) + (s11 * deltaGamma);
if abs(denominator) < eps
    gammaActual = NaN;
else
    gammaActual = deltaGamma / denominator;
end
end

function gammaSimple = oneport_simple_transmission_only(gammaMeasured, s12, s21)
denominator = s12 * s21;
if abs(denominator) < eps
    gammaSimple = NaN;
else
    gammaSimple = gammaMeasured / denominator;
end
end

function measuredSlice = get_optional_slice(dataStruct, idx)
if isempty(dataStruct)
    measuredSlice = [];
else
    measuredSlice = squeeze(dataStruct.S(:, :, idx));
end
end

function sThru = reconstruct_thru_from_embedding(sEmbedding, sLine, z0)
tEmbedding = s_to_abcd_local(sEmbedding, z0);
tLine = s_to_abcd_local(sLine, z0);
tThru = tEmbedding * tLine * tEmbedding;
sThru = abcd_to_s(tThru, z0);
end

function sActual = deembed_twport_standard_single(sMeasured, sEmbedding, z0)
tMeasured = s_to_abcd_local(sMeasured, z0);
tEmbedding = s_to_abcd_local(sEmbedding, z0);
tActual = tEmbedding \ tMeasured / tEmbedding;
sActual = abcd_to_s(tActual, z0);
end
