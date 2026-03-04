function [embeddingOut, diagnostics] = refine_embedding_with_sol(embeddingIn, bandFiles, config)
%REFINE_EMBEDDING_WITH_SOL Softly refine reflective terms of E(f) using SOL.

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
diagnostics = struct();
diagnostics.abs_s11_shift = zeros(numFreq, 1);
diagnostics.abs_s22_shift = zeros(numFreq, 1);
diagnostics.short_mag_improvement = zeros(numFreq, 1);
diagnostics.open_mag_improvement = zeros(numFreq, 1);
diagnostics.load_mag_improvement = zeros(numFreq, 1);
diagnostics.sigma_change = zeros(numFreq, 1);

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

    objective = @(x) local_refinement_objective( ...
        x, ...
        s11Base, s12, s21, s22Base, ...
        prevS11, prevS22, ...
        gammaShortMeasuredP1(idx), gammaOpenMeasuredP1(idx), gammaLoadMeasuredP1(idx), ...
        gammaShortMeasuredP2(idx), gammaOpenMeasuredP2(idx), gammaLoadMeasuredP2(idx), ...
        weights);

    xOpt = fminsearch(objective, x0, options);
    s11Refined = xOpt(1) + 1i * xOpt(2);
    s22Refined = xOpt(3) + 1i * xOpt(4);

    sParams(1, 1, idx) = s11Refined;
    sParams(2, 2, idx) = s22Refined;

    gammaShortBaseP1 = oneport_deembed_formula_port1(gammaShortMeasuredP1(idx), s11Base, s12, s21, s22Base);
    gammaOpenBaseP1 = oneport_deembed_formula_port1(gammaOpenMeasuredP1(idx), s11Base, s12, s21, s22Base);
    gammaLoadBaseP1 = oneport_deembed_formula_port1(gammaLoadMeasuredP1(idx), s11Base, s12, s21, s22Base);
    gammaShortBaseP2 = oneport_deembed_formula_port2(gammaShortMeasuredP2(idx), s11Base, s12, s21, s22Base);
    gammaOpenBaseP2 = oneport_deembed_formula_port2(gammaOpenMeasuredP2(idx), s11Base, s12, s21, s22Base);
    gammaLoadBaseP2 = oneport_deembed_formula_port2(gammaLoadMeasuredP2(idx), s11Base, s12, s21, s22Base);

    gammaShortRefinedP1 = oneport_deembed_formula_port1(gammaShortMeasuredP1(idx), s11Refined, s12, s21, s22Refined);
    gammaOpenRefinedP1 = oneport_deembed_formula_port1(gammaOpenMeasuredP1(idx), s11Refined, s12, s21, s22Refined);
    gammaLoadRefinedP1 = oneport_deembed_formula_port1(gammaLoadMeasuredP1(idx), s11Refined, s12, s21, s22Refined);
    gammaShortRefinedP2 = oneport_deembed_formula_port2(gammaShortMeasuredP2(idx), s11Refined, s12, s21, s22Refined);
    gammaOpenRefinedP2 = oneport_deembed_formula_port2(gammaOpenMeasuredP2(idx), s11Refined, s12, s21, s22Refined);
    gammaLoadRefinedP2 = oneport_deembed_formula_port2(gammaLoadMeasuredP2(idx), s11Refined, s12, s21, s22Refined);

    sigmaBase = max(svd([s11Base, s12; s21, s22Base]));
    sigmaRefined = max(svd([s11Refined, s12; s21, s22Refined]));

    diagnostics.abs_s11_shift(idx) = abs(s11Refined - s11Base);
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
end

tParams = zeros(2, 2, numFreq);
for idx = 1:numFreq
    tParams(:, :, idx) = s_to_abcd_local(sParams(:, :, idx));
end

embeddingOut = struct();
embeddingOut.S = sParams;
embeddingOut.T = tParams;
end

function score = local_refinement_objective(x, s11Base, s12, s21, s22Base, prevS11, prevS22, gammaShortMeasuredP1, gammaOpenMeasuredP1, gammaLoadMeasuredP1, gammaShortMeasuredP2, gammaOpenMeasuredP2, gammaLoadMeasuredP2, weights)
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
embeddingPassivity = max(sigmaMax - 1, 0)^2;

score = ...
    weights.open_short_passivity_p1 * openShortPassivityP1 + ...
    weights.open_short_target_p1 * openShortTargetP1 + ...
    weights.load_magnitude_p1 * loadMagnitudeP1 + ...
    weights.open_short_passivity_p2 * openShortPassivityP2 + ...
    weights.open_short_target_p2 * openShortTargetP2 + ...
    weights.load_magnitude_p2 * loadMagnitudeP2 + ...
    baselinePenalty + ...
    smoothPenalty + ...
    weights.embedding_passivity * embeddingPassivity;
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
