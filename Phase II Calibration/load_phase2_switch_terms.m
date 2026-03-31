function switchTerms = load_phase2_switch_terms(paths, targetFreq, interpMethod, pairKey, ratioMode)
%LOAD_PHASE2_SWITCH_TERMS Load, validate, and interpolate pair-specific switch ratios.

if nargin < 3 || isempty(interpMethod)
    interpMethod = 'linear';
end
if nargin < 4
    pairKey = '';
end
if nargin < 5 || isempty(ratioMode)
    ratioMode = 'b_over_a';
end

requiredFields = {'aForward', 'bForward', 'aReverse', 'bReverse'};
for idx = 1:numel(requiredFields)
    if ~isfield(paths, requiredFields{idx}) || isempty(paths.(requiredFields{idx}))
        error('Missing switch-term wave %s', requiredFields{idx});
    end
end

meta.aForward = parse_switch_wave_path(paths.aForward);
meta.bForward = parse_switch_wave_path(paths.bForward);
meta.aReverse = parse_switch_wave_path(paths.aReverse);
meta.bReverse = parse_switch_wave_path(paths.bReverse);
validate_switch_semantics(meta, upper(pairKey));

aForward = parse_phase2_wave_csv(paths.aForward);
bForward = parse_phase2_wave_csv(paths.bForward);
aReverse = parse_phase2_wave_csv(paths.aReverse);
bReverse = parse_phase2_wave_csv(paths.bReverse);

assert_frequency_match(aForward.freq, bForward.freq, 'forward');
assert_frequency_match(aReverse.freq, bReverse.freq, 'reverse');

[gammaForward, gammaReverse] = compute_switch_ratios( ...
    aForward.values, bForward.values, aReverse.values, bReverse.values, ratioMode);

switchTerms = struct();
switchTerms.paths = paths;
switchTerms.metadata = meta;
switchTerms.ratio_mode = ratioMode;
switchTerms.freq_forward_raw = aForward.freq(:);
switchTerms.freq_reverse_raw = aReverse.freq(:);
switchTerms.a_forward_raw = aForward.values(:);
switchTerms.b_forward_raw = bForward.values(:);
switchTerms.a_reverse_raw = aReverse.values(:);
switchTerms.b_reverse_raw = bReverse.values(:);
switchTerms.gamma_forward_raw = gammaForward(:);
switchTerms.gamma_reverse_raw = gammaReverse(:);
switchTerms.gamma_forward = interpolate_complex(aForward.freq, gammaForward, targetFreq, interpMethod);
switchTerms.gamma_reverse = interpolate_complex(aReverse.freq, gammaReverse, targetFreq, interpMethod);
end

function meta = parse_switch_wave_path(filePath)
[~, fileName, ext] = fileparts(filePath);
tokens = regexp([fileName, ext], ...
    '^Thru_(.+)_(Straight|Arc|Diag)_(0-67|67-115|110-170)_([ab])([1-4])S([1-4])(forward|reverse)\.csv$', ...
    'tokens', 'once');
if isempty(tokens)
    error('Unrecognized switch-term filename: %s', filePath);
end
pairRaw = upper(regexprep(tokens{1}, '[-_]', ''));
pairTokens = regexp(pairRaw, '^P([1-4])P([1-4])$', 'tokens', 'once');
if isempty(pairTokens)
    error('Could not parse pair key from switch-term filename: %s', filePath);
end

meta = struct();
meta.file = filePath;
meta.pair = sprintf('P%sP%s', pairTokens{1}, pairTokens{2});
meta.geometry = upper(tokens{2});
meta.band = tokens{3};
meta.wave_kind = lower(tokens{4});
meta.receiver = sprintf('P%s', tokens{5});
meta.stimulus = sprintf('P%s', tokens{6});
meta.direction = lower(tokens{7});
end

function validate_switch_semantics(meta, pairKey)
ports = phase2_pair_key_to_ports(pairKey);
portA = sprintf('P%d', ports(1));
portB = sprintf('P%d', ports(2));

validate_one(meta.aForward, pairKey, 'a', portB, portA, 'forward');
validate_one(meta.bForward, pairKey, 'b', portB, portA, 'forward');
validate_one(meta.aReverse, pairKey, 'a', portA, portB, 'reverse');
validate_one(meta.bReverse, pairKey, 'b', portA, portB, 'reverse');
end

function validate_one(meta, pairKey, waveKind, receiver, stimulus, direction)
if ~strcmpi(meta.pair, pairKey)
    error('Switch-term pair mismatch for %s: expected %s, found %s', meta.file, pairKey, meta.pair);
end
if ~strcmpi(meta.wave_kind, waveKind)
    error('Switch-term wave kind mismatch for %s: expected %s, found %s', meta.file, waveKind, meta.wave_kind);
end
if ~strcmpi(meta.receiver, receiver)
    error('Switch-term receiver mismatch for %s: expected %s, found %s', meta.file, receiver, meta.receiver);
end
if ~strcmpi(meta.stimulus, stimulus)
    error('Switch-term stimulus mismatch for %s: expected %s, found %s', meta.file, stimulus, meta.stimulus);
end
if ~strcmpi(meta.direction, direction)
    error('Switch-term direction mismatch for %s: expected %s, found %s', meta.file, direction, meta.direction);
end
end

function assert_frequency_match(freqA, freqB, tag)
if numel(freqA) ~= numel(freqB) || any(abs(freqA(:) - freqB(:)) > 1)
    error('Switch-term %s wave frequencies do not match.', tag);
end
end

function [gammaForward, gammaReverse] = compute_switch_ratios(aForward, bForward, aReverse, bReverse, ratioMode)
switch lower(ratioMode)
    case 'b_over_a'
        gammaForward = bForward ./ aForward;
        gammaReverse = bReverse ./ aReverse;
    case 'a_over_b'
        gammaForward = aForward ./ bForward;
        gammaReverse = aReverse ./ bReverse;
    otherwise
        error('Unsupported switch ratio mode: %s', ratioMode);
end
end

function values = interpolate_complex(freqIn, valuesIn, freqOut, interpMethod)
realPart = interp1(freqIn, real(valuesIn), freqOut, interpMethod, 'extrap');
imagPart = interp1(freqIn, imag(valuesIn), freqOut, interpMethod, 'extrap');
values = complex(realPart, imagPart);
end
