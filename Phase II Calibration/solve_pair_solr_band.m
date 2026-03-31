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
switchTerms = load_phase2_switch_terms(bandData.switch_terms.(pairKey).(geomKey), freq, config.switch_interp_method);
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
    SthruSwitch, gammaShort, gammaOpen, gammaLoad, 0);

SthruCorrected = nan(size(SthruSwitch));
for idxFreq = 1:numel(freq)
    SthruCorrected(:, :, idxFreq) = apply_error_correction_8term_local(SthruSwitch(:, :, idxFreq), errorTerms, idxFreq, config.cal_den_floor);
end

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
pairResult.reference_raw = SthruRaw;
pairResult.reference_switch_corrected = SthruSwitch;
pairResult.reference_corrected = SthruCorrected;
pairResult.reference_metrics = summarize_pair_metrics(SthruCorrected);
pairResult.standards_measured_port1 = measPort1;
pairResult.standards_measured_port2 = measPort2;
pairResult.standards_corrected = correctedStandards;
pairResult.standard_residuals = standardResiduals;
end

function gamma = interpolate_reference_gamma(reference, freq)
gamma = interp1(reference.freq, reference.gamma, freq, 'pchip', 'extrap');
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
