function switchTerms = load_phase2_switch_terms(paths, targetFreq, interpMethod)
%LOAD_PHASE2_SWITCH_TERMS Load and interpolate pair-specific switch ratios.

if nargin < 3 || isempty(interpMethod)
    interpMethod = 'linear';
end

requiredFields = {'aForward', 'bForward', 'aReverse', 'bReverse'};
for idx = 1:numel(requiredFields)
    if ~isfield(paths, requiredFields{idx}) || isempty(paths.(requiredFields{idx}))
        error('Missing switch-term wave %s', requiredFields{idx});
    end
end

aForward = parse_phase2_wave_csv(paths.aForward);
bForward = parse_phase2_wave_csv(paths.bForward);
aReverse = parse_phase2_wave_csv(paths.aReverse);
bReverse = parse_phase2_wave_csv(paths.bReverse);

gammaForward = bForward.values ./ aForward.values;
gammaReverse = bReverse.values ./ aReverse.values;

switchTerms = struct();
switchTerms.paths = paths;
switchTerms.freq_raw = aForward.freq;
switchTerms.gamma_forward_raw = gammaForward(:);
switchTerms.gamma_reverse_raw = gammaReverse(:);
switchTerms.gamma_forward = interpolate_complex(aForward.freq, gammaForward, targetFreq, interpMethod);
switchTerms.gamma_reverse = interpolate_complex(aReverse.freq, gammaReverse, targetFreq, interpMethod);
end

function values = interpolate_complex(freqIn, valuesIn, freqOut, interpMethod)
realPart = interp1(freqIn, real(valuesIn), freqOut, interpMethod, 'extrap');
imagPart = interp1(freqIn, imag(valuesIn), freqOut, interpMethod, 'extrap');
values = complex(realPart, imagPart);
end
