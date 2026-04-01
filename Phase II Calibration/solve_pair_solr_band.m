function pairResult = solve_pair_solr_band(config, pairEntry, bandName, bandData, referenceStandards)
%SOLVE_PAIR_SOLR_BAND Solve one active-pair bandwise SOLR model.

pairKey = upper(pairEntry.pair);
geomKey = upper(pairEntry.geometry);
ports = phase2_pair_key_to_ports(pairKey);
portA = sprintf('P%d', ports(1));
portB = sprintf('P%d', ports(2));
bandKey = matlab.lang.makeValidName(strrep(bandName, '-', '_'));

thruNet = bandData.reference_thrus.(pairKey).(geomKey);
freq = thruNet.freq(:);
SthruRaw = phase2_extract_pair_submatrix(thruNet.S, pairKey);
switchTerms = load_phase2_switch_terms( ...
    bandData.switch_terms.(pairKey).(geomKey), freq, config.switch_interp_method, pairKey, config.switch_ratio_mode);
[SthruSwitch, switchDiag] = apply_switch_correction_2port(SthruRaw, switchTerms.gamma_forward, switchTerms.gamma_reverse, config.switch_den_floor);

stdNames = {'Short', 'Open', 'Load'};
measPort1 = struct();
measPort2 = struct();
correctedStandards = struct();
standardResiduals = struct();

for idxStd = 1:numel(stdNames)
    stdName = stdNames{idxStd};

    stdNetA = bandData.standards.(stdName).(portA);
    stdNetB = bandData.standards.(stdName).(portB);
    SstdARaw = phase2_extract_pair_submatrix(stdNetA.S, pairKey);
    SstdBRaw = phase2_extract_pair_submatrix(stdNetB.S, pairKey);

    SstdASwitch = apply_switch_correction_2port(SstdARaw, switchTerms.gamma_forward, switchTerms.gamma_reverse, config.switch_den_floor);
    SstdBSwitch = apply_switch_correction_2port(SstdBRaw, switchTerms.gamma_forward, switchTerms.gamma_reverse, config.switch_den_floor);

    measPort1.(stdName) = squeeze(SstdASwitch(1, 1, :));
    measPort2.(stdName) = squeeze(SstdBSwitch(2, 2, :));

    correctedStandards.(stdName) = nan(2, nFreq(freq));
    standardResiduals.(stdName) = nan(2, nFreq(freq));
end

gammaShort = interpolate_reference_gamma(referenceStandards.Short, freq);
gammaOpen = interpolate_reference_gamma(referenceStandards.Open, freq);
gammaLoad = interpolate_reference_gamma(referenceStandards.Load, freq);

errorTerms = calculate_error_terms_solr_from_gamma(freq, ...
    measPort1.Short, measPort1.Open, measPort1.Load, ...
    measPort2.Short, measPort2.Open, measPort2.Load, ...
    SthruSwitch, gammaShort, gammaOpen, gammaLoad, 0, config.branch_stabilization);

phase1Target = load_phase1_thru_target(config, geomKey, pairKey, bandName);
[errorTerms, anchorDiagnostics] = apply_transmission_anchor_if_enabled( ...
    config, errorTerms, SthruSwitch, phase1Target);
SthruCorrected = correct_network_from_error_terms(SthruSwitch, errorTerms, config.cal_den_floor);

referenceMap = struct('Short', gammaShort, 'Open', gammaOpen, 'Load', gammaLoad);
for idxStd = 1:numel(stdNames)
    stdName = stdNames{idxStd};
    stdNetA = bandData.standards.(stdName).(portA);
    stdNetB = bandData.standards.(stdName).(portB);
    SstdARaw = phase2_extract_pair_submatrix(stdNetA.S, pairKey);
    SstdBRaw = phase2_extract_pair_submatrix(stdNetB.S, pairKey);
    SstdASwitch = apply_switch_correction_2port(SstdARaw, switchTerms.gamma_forward, switchTerms.gamma_reverse, config.switch_den_floor);
    SstdBSwitch = apply_switch_correction_2port(SstdBRaw, switchTerms.gamma_forward, switchTerms.gamma_reverse, config.switch_den_floor);

    correctedA = nan(numel(freq), 1);
    correctedB = nan(numel(freq), 1);
    for idxFreq = 1:numel(freq)
        correctedMatA = apply_error_correction_8term_local(SstdASwitch(:, :, idxFreq), errorTerms, idxFreq, config.cal_den_floor);
        correctedMatB = apply_error_correction_8term_local(SstdBSwitch(:, :, idxFreq), errorTerms, idxFreq, config.cal_den_floor);
        correctedA(idxFreq) = correctedMatA(1, 1);
        correctedB(idxFreq) = correctedMatB(2, 2);
    end
    correctedStandards.(stdName) = [correctedA.'; correctedB.'];
    refGamma = referenceMap.(stdName)(:).';
    standardResiduals.(stdName) = [abs(correctedA.' - refGamma); abs(correctedB.' - refGamma)];
end

pairResult = struct();
pairResult.band = bandName;
pairResult.band_key = bandKey;
pairResult.pair = pairKey;
pairResult.geometry = geomKey;
pairResult.ports = ports;
pairResult.freq = freq;
pairResult.switch_terms = switchTerms;
pairResult.switch_diagnostics = switchDiag;
pairResult.error_terms = errorTerms;
pairResult.transmission_anchor = anchorDiagnostics;
pairResult.reference_raw = SthruRaw;
pairResult.reference_switch_corrected = SthruSwitch;
pairResult.reference_corrected = SthruCorrected;
pairResult.reference_phase1_target = phase1Target;
pairResult.reference_metrics = summarize_pair_metrics(SthruCorrected);
pairResult.standards_measured_port1 = measPort1;
pairResult.standards_measured_port2 = measPort2;
pairResult.standards_corrected = correctedStandards;
pairResult.standard_residuals = standardResiduals;
end

function gamma = interpolate_reference_gamma(reference, freq)
gamma = interp1(reference.freq, reference.gamma, freq, 'pchip', 'extrap');
end

function corrected = correct_network_from_error_terms(Sraw, errorTerms, denominatorFloor)
corrected = nan(size(Sraw));
for idxFreq = 1:size(Sraw, 3)
    corrected(:, :, idxFreq) = apply_error_correction_8term_local(Sraw(:, :, idxFreq), errorTerms, idxFreq, denominatorFloor);
end
end

function target = load_phase1_thru_target(config, geomKey, pairKey, bandName)
switch upper(geomKey)
    case 'STRAIGHT'
        if strcmpi(pairKey, 'P3P4')
            baseName = sprintf('EXTRACTED_THRU_P3P4_STRAIGHT_%s.mat', bandName);
        else
            baseName = sprintf('EXTRACTED_THRU_P1P2_STRAIGHT_%s.mat', bandName);
        end
    case 'ARC'
        baseName = sprintf('EXTRACTED_THRU_P1P4_ARC_%s.mat', bandName);
    case {'DIAG', 'DIAGONAL'}
        baseName = sprintf('EXTRACTED_THRU_P1P4_DIAGONAL_%s.mat', bandName);
    otherwise
        error('Unsupported geometry %s for Phase I target lookup.', geomKey);
end

loaded = load(fullfile(config.phase1_output_mat_dir, baseName));
target = loaded.reciprocalResult;
end

function [errorTerms, diagnostics] = apply_transmission_anchor_if_enabled(config, errorTerms, Sswitch, phase1Target)
diagnostics = struct('enabled', false, 'mode', '', 'anchor_ratio', NaN, 'anchor_ratio_db', NaN, ...
    'anchor_ratio_deg', NaN, 'median_pre_s21_db', NaN, 'median_pre_s12_db', NaN, ...
    'median_target_s21_db', NaN, 'median_target_s12_db', NaN);

if ~isfield(config, 'transmission_anchor') || ~config.transmission_anchor.enabled
    return;
end

Spre = correct_network_from_error_terms(Sswitch, errorTerms, config.cal_den_floor);
preS21 = squeeze(Spre(2, 1, :));
preS12 = squeeze(Spre(1, 2, :));
targetS21 = squeeze(phase1Target.S(2, 1, :));
targetS12 = squeeze(phase1Target.S(1, 2, :));

[ratio, ratioDiag] = estimate_transmission_anchor_ratio(preS21, preS12, targetS21, targetS12, config.transmission_anchor);
if ~isfinite(real(ratio)) || ~isfinite(imag(ratio)) || abs(ratio) < 1e-9
    diagnostics.mode = config.transmission_anchor.mode;
    diagnostics.anchor_ratio = ratio;
    diagnostics.anchor_ratio_db = ratioDiag.ratio_db;
    diagnostics.anchor_ratio_deg = ratioDiag.ratio_deg;
    diagnostics.median_pre_s21_db = ratioDiag.median_pre_s21_db;
    diagnostics.median_pre_s12_db = ratioDiag.median_pre_s12_db;
    diagnostics.median_target_s21_db = ratioDiag.median_target_s21_db;
    diagnostics.median_target_s12_db = ratioDiag.median_target_s12_db;
    return;
end

errorTerms.t21 = errorTerms.t21 .* ratio;
errorTerms.t12 = errorTerms.t12 .* ratio;

diagnostics.enabled = true;
diagnostics.mode = config.transmission_anchor.mode;
diagnostics.anchor_ratio = ratio;
diagnostics.anchor_ratio_db = ratioDiag.ratio_db;
diagnostics.anchor_ratio_deg = ratioDiag.ratio_deg;
diagnostics.median_pre_s21_db = ratioDiag.median_pre_s21_db;
diagnostics.median_pre_s12_db = ratioDiag.median_pre_s12_db;
diagnostics.median_target_s21_db = ratioDiag.median_target_s21_db;
diagnostics.median_target_s12_db = ratioDiag.median_target_s12_db;
end

function [ratio, diagnostics] = estimate_transmission_anchor_ratio(preS21, preS12, targetS21, targetS12, anchorCfg)
ratio = NaN;
diagnostics = struct('ratio_db', NaN, 'ratio_deg', NaN, ...
    'median_pre_s21_db', median(local_mag_db(preS21), 'omitnan'), ...
    'median_pre_s12_db', median(local_mag_db(preS12), 'omitnan'), ...
    'median_target_s21_db', median(local_mag_db(targetS21), 'omitnan'), ...
    'median_target_s12_db', median(local_mag_db(targetS12), 'omitnan'));

minTargetMag = 10^(anchorCfg.min_target_mag_db / 20);
valid21 = isfinite(preS21) & isfinite(targetS21) & abs(targetS21) >= minTargetMag;
valid12 = isfinite(preS12) & isfinite(targetS12) & abs(targetS12) >= minTargetMag;
ratios = [preS21(valid21) ./ targetS21(valid21); preS12(valid12) ./ targetS12(valid12)];
if isempty(ratios)
    return;
end

magDb = 20 * log10(max(abs(ratios), 1e-12));
trimFrac = max(min(anchorCfg.trim_frac, 0.45), 0);
if numel(magDb) > 4 && trimFrac > 0
    qLo = quantile(magDb, trimFrac);
    qHi = quantile(magDb, 1 - trimFrac);
    keep = magDb >= qLo & magDb <= qHi;
    ratios = ratios(keep);
    magDb = magDb(keep);
end
if isempty(ratios)
    return;
end

ratioDb = median(magDb, 'omitnan');
ratioDb = max(min(ratioDb, anchorCfg.max_gain_db), -anchorCfg.max_gain_db);
ratioPhase = angle(mean(exp(1j * angle(ratios)), 'omitnan'));

ratio = 10^(ratioDb / 20) * exp(1j * ratioPhase);
diagnostics.ratio_db = ratioDb;
diagnostics.ratio_deg = ratioPhase * 180 / pi;
end

function values = local_mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function metrics = summarize_pair_metrics(S)
s21 = squeeze(S(2, 1, :));
s12 = squeeze(S(1, 2, :));
s11 = squeeze(S(1, 1, :));
s22 = squeeze(S(2, 2, :));
metrics = struct();
metrics.mean_recip_mismatch = mean(abs(s21 - s12), 'omitnan');
metrics.max_recip_mismatch = max(abs(s21 - s12), [], 'omitnan');
metrics.median_insertion_loss_db = median(-20 * log10(max(abs(s21), 1e-12)), 'omitnan');
metrics.max_return_loss_db = max(-20 * log10(max([abs(s11); abs(s22)], 1e-12)), [], 'omitnan');
end

function n = nFreq(freq)
n = numel(freq);
end
