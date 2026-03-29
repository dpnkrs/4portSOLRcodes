function data = read_oneport_touchstone_file(filename)
%READ_ONEPORT_TOUCHSTONE_FILE Read one-port Touchstone from .s1p or .s2p.

[~, ~, ext] = fileparts(filename);
ext = lower(ext);

switch ext
    case '.s1p'
        data = read_s1p_local(filename);
    case '.s2p'
        twoport = read_touchstone_file(filename);
        data = struct();
        data.filename = filename;
        data.freq = twoport.freq;
        data.z0 = twoport.z0;
        data.gamma = squeeze(twoport.S(1, 1, :));
    otherwise
        error('Unsupported one-port Touchstone extension: %s', ext);
end
end

function data = read_s1p_local(filename)
fid = fopen(filename, 'r');
if fid == -1
    error('Unable to open Touchstone file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
lines = lines{1};

fmt = '';
freqUnit = '';
z0 = 50;
numericRows = [];

for idx = 1:numel(lines)
    line = strtrim(lines{idx});
    if isempty(line) || startsWith(line, '!')
        continue;
    end

    if startsWith(line, '#')
        tokens = split(regexprep(line, '\s+', ' '), ' ');
        tokens = tokens(~cellfun(@isempty, tokens));
        if numel(tokens) < 5
            error('Unsupported Touchstone option line in %s', filename);
        end
        freqUnit = upper(tokens{2});
        fmt = upper(tokens{4});
        if numel(tokens) >= 6
            z0 = str2double(tokens{6});
        end
        continue;
    end

    values = sscanf(line, '%f').';
    if ~isempty(values)
        numericRows = [numericRows; values]; %#ok<AGROW>
    end
end

if size(numericRows, 2) ~= 3
    error('Expected 3 numeric columns in %s but found %d.', filename, size(numericRows, 2));
end

freq = convert_frequency_to_hz_local(numericRows(:, 1), freqUnit);
numFreq = numel(freq);
gamma = zeros(numFreq, 1);
for idx = 1:numFreq
    gamma(idx) = decode_touchstone_pair_local(numericRows(idx, 2:3), fmt);
end

data = struct();
data.filename = filename;
data.freq = freq;
data.z0 = z0;
data.gamma = gamma;
end

function freqHz = convert_frequency_to_hz_local(freq, freqUnit)
switch upper(freqUnit)
    case 'HZ'
        scale = 1;
    case 'KHZ'
        scale = 1e3;
    case 'MHZ'
        scale = 1e6;
    case 'GHZ'
        scale = 1e9;
    otherwise
        error('Unsupported Touchstone frequency unit: %s', freqUnit);
end
freqHz = freq * scale;
end

function value = decode_touchstone_pair_local(rawPair, fmt)
switch upper(fmt)
    case 'RI'
        value = rawPair(1) + 1i * rawPair(2);
    case 'MA'
        value = rawPair(1) .* exp(1i * deg2rad(rawPair(2)));
    case 'DB'
        value = 10 .^ (rawPair(1) / 20) .* exp(1i * deg2rad(rawPair(2)));
    otherwise
        error('Unsupported Touchstone data format: %s', fmt);
end
end