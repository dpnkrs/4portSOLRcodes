function [embeddingOut, diagnostics] = refine_embedding_with_sol(embeddingIn, bandFiles, config, straightData, p34StraightData, lineModel, bandName)
%REFINE_EMBEDDING_WITH_SOL Weighted Stage A embedding refinement.
% Refines full S11/S12/S21/S22 per frequency with strong transmission regularization.
% Includes optional soft simulation priors (raw + virtual-mTRL).

if nargin < 7
    bandName = '';
end

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
baselineGamma = compute_baseline_deembedded_gammas(sParams, gammaShortMeasuredP1, gammaOpenMeasuredP1, gammaLoadMeasuredP1, gammaShortMeasuredP2, gammaOpenMeasuredP2, gammaLoadMeasuredP2, config);
simBand = build_sim_priors_on_band(config, freqVector, baselineGamma);

options = optimset( ...
    'Display', 'off', ...
    'MaxIter', config.embedding_refine_max_iter, ...
    'MaxFunEvals', config.embedding_refine_max_fun_evals, ...
    'TolX', 1e-6, ...
    'TolFun', 1e-6);

for idx = 1:numFreq
    s11Base = sParams(1, 1, idx);
    s12Base = sParams(1, 2, idx);
    s21Base = sParams(2, 1, idx);
    s22Base = sParams(2, 2, idx);

    if idx == 1
        x0 = [real(s11Base), imag(s11Base), real(s12Base), imag(s12Base), real(s21Base), imag(s21Base), real(s22Base), imag(s22Base)];
        prevS11 = [];
        prevS12 = [];
        prevS21 = [];
        prevS22 = [];
        prevPrevS11 = [];
        prevPrevS22 = [];
    else
        x0 = [real(sParams(1, 1, idx - 1)), imag(sParams(1, 1, idx - 1)), real(sParams(1, 2, idx - 1)), imag(sParams(1, 2, idx - 1)), real(sParams(2, 1, idx - 1)), imag(sParams(2, 1, idx - 1)), real(sParams(2, 2, idx - 1)), imag(sParams(2, 2, idx - 1))];
        prevS11 = sParams(1, 1, idx - 1);
        prevS12 = sParams(1, 2, idx - 1);
        prevS21 = sParams(2, 1, idx - 1);
        prevS22 = sParams(2, 2, idx - 1);
        if idx > 2
            prevPrevS11 = sParams(1, 1, idx - 2);
            prevPrevS22 = sParams(2, 2, idx - 2);
        else
            prevPrevS11 = [];
            prevPrevS22 = [];
        end
    end

    localWeight = frequency_weight(freqVector(idx), freqRange, config);
    edgeGate = edge_refinement_gate(freqVector(idx), freqRange, config);
    localWeight = localWeight * max(edgeGate, 0.2);

    objective = @(x) local_refinement_objective( ...
        x, ...
        s11Base, s12Base, s21Base, s22Base, ...
        prevS11, prevS12, prevS21, prevS22, prevPrevS11, prevPrevS22, ...
        gammaShortMeasuredP1(idx), gammaOpenMeasuredP1(idx), gammaLoadMeasuredP1(idx), ...
        gammaShortMeasuredP2(idx), gammaOpenMeasuredP2(idx), gammaLoadMeasuredP2(idx), ...
        squeeze(straightData.S(:, :, idx)), ...
        squeeze(lineModel.S(:, :, idx)), ...
        get_optional_slice(p34StraightData, idx), ...
        simBand, idx, freqVector(idx), ...
        weights, localWeight, config.z0, ...
        config.embedding_refine_denominator_floor, ...
        config.embedding_refine_s12s21_trust_radius, ...
        config.embedding_refine_highband_start_hz, ...
        config.embedding_refine_highband_damp_start_hz, ...
        config.embedding_refine_highband_baseline_s22_scale, ...
        config.embedding_refine_highband_strong_s22_start_hz, ...
        config.embedding_refine_highband_strong_s22_scale, ...
        config.embedding_refine_highband_open_short_p2_scale, ...
        config.embedding_refine_highband_load_p2_scale, ...
        config.embedding_refine_transition_windows_hz, ...
        config.embedding_refine_transition_p34_scale, ...
        config.embedding_refine_transition_trust_scale, ...
        config.embedding_refine_transition_baseline_s11_scale, ...
        config.embedding_refine_transition_baseline_s22_scale, ...
        config.embedding_refine_transition_soft_hz, ...
        config.embedding_refine_0_67_seconddiff_window_hz, ...
        config.embedding_refine_0_67_seconddiff_weight_s11, ...
        config.embedding_refine_0_67_seconddiff_weight_s22, ...
        config.embedding_refine_0_67_s22_focus_window_hz, ...
        config.embedding_refine_0_67_s22_focus_weight, ...
        config.embedding_refine_0_67_disable_sim_window_hz, ...
        bandName);

    xOpt = fminsearch(objective, x0, options);
    s11Refined = xOpt(1) + 1i * xOpt(2);
    s12Refined = xOpt(3) + 1i * xOpt(4);
    s21Refined = xOpt(5) + 1i * xOpt(6);
    s22Refined = xOpt(7) + 1i * xOpt(8);

    sParams(1, 1, idx) = s11Refined;
    sParams(1, 2, idx) = s12Refined;
    sParams(2, 1, idx) = s21Refined;
    sParams(2, 2, idx) = s22Refined;

    diagnostics = update_diagnostics( ...
        diagnostics, idx, ...
        s11Base, s12Base, s21Base, s22Base, ...
        s11Refined, s12Refined, s21Refined, s22Refined, ...
        gammaShortMeasuredP1(idx), gammaOpenMeasuredP1(idx), gammaLoadMeasuredP1(idx), ...
        gammaShortMeasuredP2(idx), gammaOpenMeasuredP2(idx), gammaLoadMeasuredP2(idx), ...
        squeeze(straightData.S(:, :, idx)), ...
        squeeze(lineModel.S(:, :, idx)), ...
        get_optional_slice(p34StraightData, idx), ...
        simBand, idx, ...
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
diagnostics.sim_ref_error = NaN(numFreq, 1);
diagnostics.sim_raw_error = NaN(numFreq, 1);
end

function baselineGamma = compute_baseline_deembedded_gammas(sParams, gammaShortMeasuredP1, gammaOpenMeasuredP1, gammaLoadMeasuredP1, gammaShortMeasuredP2, gammaOpenMeasuredP2, gammaLoadMeasuredP2, config)
numFreq = size(sParams, 3);
baselineGamma = struct();
baselineGamma.p1.short = zeros(numFreq, 1);
baselineGamma.p1.open = zeros(numFreq, 1);
baselineGamma.p1.load = zeros(numFreq, 1);
baselineGamma.p2.short = zeros(numFreq, 1);
baselineGamma.p2.open = zeros(numFreq, 1);
baselineGamma.p2.load = zeros(numFreq, 1);

for ii = 1:numFreq
    s11 = sParams(1, 1, ii);
    s12 = sParams(1, 2, ii);
    s21 = sParams(2, 1, ii);
    s22 = sParams(2, 2, ii);

    baselineGamma.p1.short(ii) = oneport_deembed_formula_port1(gammaShortMeasuredP1(ii), s11, s12, s21, s22);
    baselineGamma.p1.open(ii) = oneport_deembed_formula_port1(gammaOpenMeasuredP1(ii), s11, s12, s21, s22);
    baselineGamma.p1.load(ii) = oneport_deembed_formula_port1(gammaLoadMeasuredP1(ii), s11, s12, s21, s22);

    baselineGamma.p2.short(ii) = oneport_deembed_formula_port2(gammaShortMeasuredP2(ii), s11, s12, s21, s22);
    baselineGamma.p2.open(ii) = oneport_deembed_formula_port2(gammaOpenMeasuredP2(ii), s11, s12, s21, s22);
    baselineGamma.p2.load(ii) = oneport_deembed_formula_port2(gammaLoadMeasuredP2(ii), s11, s12, s21, s22);
end

if ~isfield(config, 'simulation_prior_use_load') || ~config.simulation_prior_use_load
    baselineGamma.p1 = rmfield(baselineGamma.p1, 'load');
    baselineGamma.p2 = rmfield(baselineGamma.p2, 'load');
end
end
function simBand = build_sim_priors_on_band(config, freqVector, baselineGamma)
simBand = struct('enabled', false);
if ~isfield(config, 'simulation_priors') || ~isfield(config.simulation_priors, 'enabled') || ~config.simulation_priors.enabled
    return;
end
sim = config.simulation_priors;

fields = {'short', 'open', 'load'};
if ~isfield(config, 'simulation_prior_use_load') || ~config.simulation_prior_use_load
    fields = {'short', 'open'};
end

for i = 1:numel(fields)
    key = fields{i};
    simRef = interp1(sim.ref.(key).freq, sim.ref.(key).gamma, freqVector, 'linear', 'extrap');
    simRaw = interp1(sim.raw.(key).freq, sim.raw.(key).gamma, freqVector, 'linear', 'extrap');

    if isfield(config, 'simulation_prior_alignment') && config.simulation_prior_alignment
        robustAlign = isfield(config, 'simulation_prior_alignment_robust') && config.simulation_prior_alignment_robust;
        [a1, b1] = fit_complex_affine(simRef, baselineGamma.p1.(key), robustAlign);
        [a2, b2] = fit_complex_affine(simRef, baselineGamma.p2.(key), robustAlign);
        simBand.ref_p1.(key) = a1 * simRef + b1;
        simBand.ref_p2.(key) = a2 * simRef + b2;
    else
        simBand.ref_p1.(key) = simRef;
        simBand.ref_p2.(key) = simRef;
    end

    simBand.raw.(key) = simRaw;
end

simBand.fields = fields;
simBand.enabled = true;
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

function gate = edge_refinement_gate(freqHz, freqRange, config)
fLo = freqRange(1);
fHi = freqRange(2);
distEdge = min(freqHz - fLo, fHi - freqHz);
freezeHz = max(config.embedding_refine_edge_freeze_hz, 0);
taperHz = max(config.embedding_refine_edge_taper_hz, eps);

if distEdge <= freezeHz
    gate = 0.2;
elseif distEdge >= (freezeHz + taperHz)
    gate = 1;
else
    gate = (distEdge - freezeHz) / taperHz;
end
end

function score = local_refinement_objective(x, s11Base, s12Base, s21Base, s22Base, prevS11, prevS12, prevS21, prevS22, prevPrevS11, prevPrevS22, gammaShortMeasuredP1, gammaOpenMeasuredP1, gammaLoadMeasuredP1, gammaShortMeasuredP2, gammaOpenMeasuredP2, gammaLoadMeasuredP2, p12Measured, lineS, p34Measured, simBand, idx, freqHz, weights, localWeight, z0, denominatorFloor, trustRadius, highbandStartHz, highbandDampStartHz, highbandBaselineS22Scale, highbandStrongS22StartHz, highbandStrongS22Scale, highbandOpenShortP2Scale, highbandLoadP2Scale, transitionWindowsHz, transitionP34Scale, transitionTrustScale, transitionBaselineS11Scale, transitionBaselineS22Scale, transitionSoftHz, secondDiffWindowHz, secondDiffWeightS11, secondDiffWeightS22, secondDiffS22FocusWindowHz, secondDiffS22FocusWeight, disableSimWindowHz, bandName)
s11 = x(1) + 1i * x(2);
s12 = x(3) + 1i * x(4);
s21 = x(5) + 1i * x(6);
s22 = x(7) + 1i * x(8);

gammaShortP1 = oneport_deembed_formula_port1(gammaShortMeasuredP1, s11, s12, s21, s22);
gammaOpenP1 = oneport_deembed_formula_port1(gammaOpenMeasuredP1, s11, s12, s21, s22);
gammaLoadP1 = oneport_deembed_formula_port1(gammaLoadMeasuredP1, s11, s12, s21, s22);

gammaShortP2 = oneport_deembed_formula_port2(gammaShortMeasuredP2, s11, s12, s21, s22);
gammaOpenP2 = oneport_deembed_formula_port2(gammaOpenMeasuredP2, s11, s12, s21, s22);
gammaLoadP2 = oneport_deembed_formula_port2(gammaLoadMeasuredP2, s11, s12, s21, s22);

shortPassivityP1 = max(abs(gammaShortP1) - 1, 0)^2;
openPassivityP1 = max(abs(gammaOpenP1) - 1, 0)^2;
shortTargetP1 = (abs(gammaShortP1) - 1)^2;
openTargetP1 = (abs(gammaOpenP1) - 1)^2;
loadMagnitudeP1 = abs(gammaLoadP1)^2;

shortPassivityP2 = max(abs(gammaShortP2) - 1, 0)^2;
openPassivityP2 = max(abs(gammaOpenP2) - 1, 0)^2;
shortTargetP2 = (abs(gammaShortP2) - 1)^2;
openTargetP2 = (abs(gammaOpenP2) - 1)^2;
loadMagnitudeP2 = abs(gammaLoadP2)^2;

openBoost = 1;
if freqHz >= highbandStartHz
    openBoost = max(1, weights.open_highband_boost);
end

openShortPassivityP1 = shortPassivityP1 + openBoost * openPassivityP1;
openShortTargetP1 = shortTargetP1 + openBoost * openTargetP1;
openShortPassivityP2 = shortPassivityP2 + openBoost * openPassivityP2;
openShortTargetP2 = shortTargetP2 + openBoost * openTargetP2;

% Local transition damping around known 0-67 sub-band boundaries.
transitionWeight = 0;
for kk = 1:size(transitionWindowsHz, 1)
    transitionWeight = max(transitionWeight, smooth_window_weight(freqHz, transitionWindowsHz(kk, :), transitionSoftHz));
end

baselineS11Scale = 1;
baselineS22Scale = 1;
p34Scale = 1;
trustScale = 1;
baselineS11Scale = baselineS11Scale + transitionWeight * (transitionBaselineS11Scale - 1);
baselineS22Scale = baselineS22Scale + transitionWeight * (transitionBaselineS22Scale - 1);
p34Scale = p34Scale + transitionWeight * (transitionP34Scale - 1);
trustScale = trustScale + transitionWeight * (transitionTrustScale - 1);

% Above the high-band damping threshold, keep S22 closer to baseline and
% reduce aggressive P2 reflective forcing that was distorting open.
p2OpenShortScale = 1;
p2LoadScale = 1;
if freqHz >= highbandDampStartHz
    baselineS22Scale = max(baselineS22Scale, highbandBaselineS22Scale);
    p2OpenShortScale = highbandOpenShortP2Scale;
    p2LoadScale = highbandLoadP2Scale;
end
if freqHz >= highbandStrongS22StartHz
    baselineS22Scale = max(baselineS22Scale, highbandStrongS22Scale);
end

baselinePenalty = ...
    weights.baseline_s11 * baselineS11Scale * abs(s11 - s11Base)^2 + ...
    weights.baseline_s12 * abs(s12 - s12Base)^2 + ...
    weights.baseline_s21 * abs(s21 - s21Base)^2 + ...
    weights.baseline_s22 * baselineS22Scale * abs(s22 - s22Base)^2;

smoothPenalty = 0;
if ~isempty(prevS11)
    smoothPenalty = ...
        weights.smoothness_s11 * abs(s11 - prevS11)^2 + ...
        weights.smoothness_s12 * abs(s12 - prevS12)^2 + ...
        weights.smoothness_s21 * abs(s21 - prevS21)^2 + ...
        weights.smoothness_s22 * abs(s22 - prevS22)^2;
end

trustPenalty = ...
    max(abs(s12 - s12Base) - trustRadius, 0)^2 + ...
    max(abs(s21 - s21Base) - trustRadius, 0)^2;

secondDiffPenalty = 0;
if strcmpi(strtrim(bandName), '0-67') && ~isempty(prevPrevS11) && freqHz >= secondDiffWindowHz(1) && freqHz <= secondDiffWindowHz(2)
    secondDiffPenalty = ...
        secondDiffWeightS11 * abs(s11 - 2 * prevS11 + prevPrevS11)^2 + ...
        secondDiffWeightS22 * abs(s22 - 2 * prevS22 + prevPrevS22)^2;
    if freqHz >= secondDiffS22FocusWindowHz(1) && freqHz <= secondDiffS22FocusWindowHz(2)
        secondDiffPenalty = secondDiffPenalty + ...
            secondDiffS22FocusWeight * abs(s22 - 2 * prevS22 + prevPrevS22)^2;
    end
end

sMatrix = [s11, s12; s21, s22];
sigmaMax = max(svd(sMatrix));
sigmaExcess = max(sigmaMax - 1, 0);
embeddingPassivity = ...
    weights.embedding_passivity_mean * sigmaExcess^2 + ...
    weights.embedding_passivity_local * sigmaExcess^4;
reciprocityPenalty = weights.reciprocity * abs(s12 - s21)^2;

thruRecon = reconstruct_thru_from_embedding(sMatrix, lineS, z0);
thruPenalty = norm(thruRecon - p12Measured, 'fro')^2;

p34Penalty = 0;
if ~isempty(p34Measured)
    p34Actual = deembed_twport_standard_single(p34Measured, sMatrix, z0);
    p34Penalty = ...
        weights.p34_line_match * norm(p34Actual - lineS, 'fro')^2 + ...
        weights.p34_return_loss * (abs(p34Actual(1, 1))^2 + abs(p34Actual(2, 2))^2);
    p34Penalty = p34Scale * p34Penalty;
end

p1DenShort = abs((s12 * s21) + (s22 * (gammaShortMeasuredP1 - s11)));
p1DenOpen = abs((s12 * s21) + (s22 * (gammaOpenMeasuredP1 - s11)));
p1DenLoad = abs((s12 * s21) + (s22 * (gammaLoadMeasuredP1 - s11)));
p2DenShort = abs((s12 * s21) + (s11 * (gammaShortMeasuredP2 - s22)));
p2DenOpen = abs((s12 * s21) + (s11 * (gammaOpenMeasuredP2 - s22)));
p2DenLoad = abs((s12 * s21) + (s11 * (gammaLoadMeasuredP2 - s22)));

denominatorPenalty = ...
    max(denominatorFloor - p1DenShort, 0)^2 + ...
    max(denominatorFloor - p1DenOpen, 0)^2 + ...
    max(denominatorFloor - p1DenLoad, 0)^2 + ...
    max(denominatorFloor - p2DenShort, 0)^2 + ...
    max(denominatorFloor - p2DenOpen, 0)^2 + ...
    max(denominatorFloor - p2DenLoad, 0)^2;

simRefPenalty = 0;
simRawPenalty = 0;
if simBand.enabled
    for jj = 1:numel(simBand.fields)
        key = simBand.fields{jj};
        if strcmp(key, 'short')
            gP1 = gammaShortP1;
            gP2 = gammaShortP2;
        elseif strcmp(key, 'open')
            gP1 = gammaOpenP1;
            gP2 = gammaOpenP2;
        else
            gP1 = gammaLoadP1;
            gP2 = gammaLoadP2;
        end

        refP1 = simBand.ref_p1.(key)(idx);
        refP2 = simBand.ref_p2.(key)(idx);
        simRefPenalty = simRefPenalty + abs(gP1 - refP1)^2 + abs(gP2 - refP2)^2;

        if isfield(simBand, 'raw') && isfield(simBand.raw, key)
            rawVal = simBand.raw.(key)(idx);
            predRawP1 = oneport_forward_formula_port1(refP1, s11, s12, s21, s22);
            predRawP2 = oneport_forward_formula_port2(refP2, s11, s12, s21, s22);
            simRawPenalty = simRawPenalty + abs(predRawP1 - rawVal)^2 + abs(predRawP2 - rawVal)^2;
        end
    end

    normCount = max(2 * numel(simBand.fields), 1);
    simRefPenalty = simRefPenalty / normCount;
    simRawPenalty = simRawPenalty / normCount;
end

simScale = 1;
if strcmpi(strtrim(bandName), '0-67') && freqHz >= disableSimWindowHz(1) && freqHz <= disableSimWindowHz(2)
    simScale = 0;
end

score = ...
    localWeight * weights.open_short_passivity_p1 * openShortPassivityP1 + ...
    localWeight * weights.open_short_target_p1 * openShortTargetP1 + ...
    localWeight * weights.load_magnitude_p1 * loadMagnitudeP1 + ...
    localWeight * weights.open_short_passivity_p2 * p2OpenShortScale * openShortPassivityP2 + ...
    localWeight * weights.open_short_target_p2 * p2OpenShortScale * openShortTargetP2 + ...
    localWeight * weights.load_magnitude_p2 * p2LoadScale * loadMagnitudeP2 + ...
    localWeight * weights.thru_reconstruct * thruPenalty + ...
    localWeight * p34Penalty + ...
    baselinePenalty + ...
    smoothPenalty + ...
    localWeight * embeddingPassivity + ...
    localWeight * reciprocityPenalty + ...
    localWeight * weights.denominator_conditioning * denominatorPenalty + ...
    localWeight * weights.trust_region_s12s21 * trustScale * trustPenalty + ...
    secondDiffPenalty + ...
    simScale * localWeight * weights.simulation_ref * simRefPenalty + ...
    simScale * localWeight * weights.simulation_raw * simRawPenalty;
end

function w = smooth_window_weight(freqHz, windowHz, softHz)
startHz = windowHz(1);
stopHz = windowHz(2);
if stopHz <= startHz
    w = 0;
    return;
end

if softHz <= 0
    w = double(freqHz >= startHz && freqHz <= stopHz);
    return;
end

left0 = startHz - softHz;
left1 = startHz + softHz;
right0 = stopHz - softHz;
right1 = stopHz + softHz;

if freqHz <= left0 || freqHz >= right1
    w = 0;
elseif freqHz >= left1 && freqHz <= right0
    w = 1;
elseif freqHz > left0 && freqHz < left1
    t = (freqHz - left0) / max(left1 - left0, eps);
    w = 0.5 - 0.5 * cos(pi * t);
else
    t = (right1 - freqHz) / max(right1 - right0, eps);
    w = 0.5 - 0.5 * cos(pi * t);
end
end

function diagnostics = update_diagnostics(diagnostics, idx, s11Base, s12Base, s21Base, s22Base, s11Refined, s12Refined, s21Refined, s22Refined, gammaShortMeasuredP1, gammaOpenMeasuredP1, gammaLoadMeasuredP1, gammaShortMeasuredP2, gammaOpenMeasuredP2, gammaLoadMeasuredP2, p12Measured, lineS, p34Measured, simBand, simIdx, z0)
gammaShortBaseP1 = oneport_deembed_formula_port1(gammaShortMeasuredP1, s11Base, s12Base, s21Base, s22Base);
gammaOpenBaseP1 = oneport_deembed_formula_port1(gammaOpenMeasuredP1, s11Base, s12Base, s21Base, s22Base);
gammaLoadBaseP1 = oneport_deembed_formula_port1(gammaLoadMeasuredP1, s11Base, s12Base, s21Base, s22Base);
gammaShortBaseP2 = oneport_deembed_formula_port2(gammaShortMeasuredP2, s11Base, s12Base, s21Base, s22Base);
gammaOpenBaseP2 = oneport_deembed_formula_port2(gammaOpenMeasuredP2, s11Base, s12Base, s21Base, s22Base);
gammaLoadBaseP2 = oneport_deembed_formula_port2(gammaLoadMeasuredP2, s11Base, s12Base, s21Base, s22Base);

gammaShortRefinedP1 = oneport_deembed_formula_port1(gammaShortMeasuredP1, s11Refined, s12Refined, s21Refined, s22Refined);
gammaOpenRefinedP1 = oneport_deembed_formula_port1(gammaOpenMeasuredP1, s11Refined, s12Refined, s21Refined, s22Refined);
gammaLoadRefinedP1 = oneport_deembed_formula_port1(gammaLoadMeasuredP1, s11Refined, s12Refined, s21Refined, s22Refined);
gammaShortRefinedP2 = oneport_deembed_formula_port2(gammaShortMeasuredP2, s11Refined, s12Refined, s21Refined, s22Refined);
gammaOpenRefinedP2 = oneport_deembed_formula_port2(gammaOpenMeasuredP2, s11Refined, s12Refined, s21Refined, s22Refined);
gammaLoadRefinedP2 = oneport_deembed_formula_port2(gammaLoadMeasuredP2, s11Refined, s12Refined, s21Refined, s22Refined);

gammaShortSimpleP1 = oneport_simple_transmission_only(gammaShortMeasuredP1, s12Refined, s21Refined);
gammaOpenSimpleP1 = oneport_simple_transmission_only(gammaOpenMeasuredP1, s12Refined, s21Refined);
gammaLoadSimpleP1 = oneport_simple_transmission_only(gammaLoadMeasuredP1, s12Refined, s21Refined);

sigmaBase = max(svd([s11Base, s12Base; s21Base, s22Base]));
sigmaRefined = max(svd([s11Refined, s12Refined; s21Refined, s22Refined]));

sEmbedRefined = [s11Refined, s12Refined; s21Refined, s22Refined];
sThruRecon = reconstruct_thru_from_embedding(sEmbedRefined, lineS, z0);

diagnostics.abs_s11_shift(idx) = abs(s11Refined - s11Base);
diagnostics.abs_s12_shift(idx) = abs(s12Refined - s12Base);
    diagnostics.abs_s21_shift(idx) = abs(s21Refined - s21Base);
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

if simBand.enabled
    simRefAccum = 0;
    simRawAccum = 0;
    for jj = 1:numel(simBand.fields)
        key = simBand.fields{jj};
        if strcmp(key, 'short')
            gP1 = gammaShortRefinedP1;
            gP2 = gammaShortRefinedP2;
        elseif strcmp(key, 'open')
            gP1 = gammaOpenRefinedP1;
            gP2 = gammaOpenRefinedP2;
        else
            gP1 = gammaLoadRefinedP1;
            gP2 = gammaLoadRefinedP2;
        end

        refP1 = simBand.ref_p1.(key)(simIdx);
        refP2 = simBand.ref_p2.(key)(simIdx);
        simRefAccum = simRefAccum + abs(gP1 - refP1)^2 + abs(gP2 - refP2)^2;

        if isfield(simBand, 'raw') && isfield(simBand.raw, key)
            rawVal = simBand.raw.(key)(simIdx);
            predRawP1 = oneport_forward_formula_port1(refP1, s11Refined, s12Refined, s21Refined, s22Refined);
            predRawP2 = oneport_forward_formula_port2(refP2, s11Refined, s12Refined, s21Refined, s22Refined);
            simRawAccum = simRawAccum + abs(predRawP1 - rawVal)^2 + abs(predRawP2 - rawVal)^2;
        end
    end

    normCount = max(2 * numel(simBand.fields), 1);
    diagnostics.sim_ref_error(idx) = simRefAccum / normCount;
    diagnostics.sim_raw_error(idx) = simRawAccum / normCount;
end
end

function gammaMeasured = oneport_forward_formula_port1(gammaActual, s11, s12, s21, s22)
gammaMeasured = s11 + (s12 * s21 * gammaActual) / (1 - s22 * gammaActual);
end

function gammaMeasured = oneport_forward_formula_port2(gammaActual, s11, s12, s21, s22)
gammaMeasured = s22 + (s12 * s21 * gammaActual) / (1 - s11 * gammaActual);
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







function [a, b] = fit_complex_affine(x, y, robustMode)
x = x(:);
y = y(:);
valid = isfinite(real(x)) & isfinite(imag(x)) & isfinite(real(y)) & isfinite(imag(y));
x = x(valid);
y = y(valid);

if numel(x) < 3
    a = 1;
    b = 0;
    return;
end

A = [x, ones(size(x))];
if robustMode
    r = abs(y - median(y));
    s = median(r) + eps;
    w = 1 ./ max(1, r / (4 * s));
    W = diag(w);
    p = (A' * W * A) \ (A' * W * y);
else
    p = A \ y;
end

a = p(1);
b = p(2);
if ~isfinite(real(a)) || ~isfinite(imag(a)) || ~isfinite(real(b)) || ~isfinite(imag(b))
    a = 1;
    b = 0;
end
end




