function network = read_touchstone_nport(filePath)
%READ_TOUCHSTONE_NPORT Minimal Touchstone v1 reader for S-parameter files.

if ~exist(filePath, 'file')
    error('Touchstone file not found: %s', filePath);
end

[~, ~, ext] = fileparts(filePath);
tokens = regexp(lower(ext), '^\.s(\d+)p$', 'tokens', 'once');
if isempty(tokens)
    error('Unsupported Touchstone extension: %s', ext);
end
nPorts = str2double(tokens{1});

rawLines = string(splitlines(fileread(filePath)));
rawLines = strtrim(rawLines);
rawLines(rawLines == "") = [];

optionLine = "";
dataStart = 0;
for idx = 1:numel(rawLines)
    line = rawLines(idx);
    if startsWith(line, "#")
        optionLine = line;
    elseif ~startsWith(line, "!")
        dataStart = idx;
        break;
    end
end

if optionLine == "" || dataStart == 0
    error('Invalid Touchstone file (missing option/data lines): %s', filePath);
end

opt = regexp(char(optionLine), '#\s+(\S+)\s+(\S+)\s+(\S+)\s+R\s+([0-9eE\+\-\.]+)', 'tokens', 'once');
if isempty(opt)
    error('Unsupported Touchstone option line in %s: %s', filePath, optionLine);
end

freqUnit = upper(opt{1});
paramType = upper(opt{2});
dataFormat = upper(opt{3});
z0 = str2double(opt{4});

if ~strcmp(paramType, 'S')
    error('Only S-parameter Touchstone files are supported: %s', filePath);
end

valuesPerFreq = 2 * nPorts * nPorts;
numericValues = [];
for idx = dataStart:numel(rawLines)
    line = char(rawLines(idx));
    if isempty(line) || startsWith(strtrim(line), '!')
        continue;
    end
    nums = sscanf(line, '%f').';
    if ~isempty(nums)
        numericValues = [numericValues, nums]; %#ok<AGROW>
    end
end

blockLen = valuesPerFreq + 1;
if mod(numel(numericValues), blockLen) ~= 0
    error('Touchstone data length mismatch in %s', filePath);
end

nFreq = numel(numericValues) / blockLen;
numericValues = reshape(numericValues, blockLen, nFreq).';
freq = numericValues(:, 1);
data = numericValues(:, 2:end);

switch upper(freqUnit)
    case 'HZ'
        freqScale = 1;
    case 'KHZ'
        freqScale = 1e3;
    case 'MHZ'
        freqScale = 1e6;
    case 'GHZ'
        freqScale = 1e9;
    otherwise
        error('Unsupported frequency unit %s in %s', freqUnit, filePath);
end
freq = freq * freqScale;

S = zeros(nPorts, nPorts, nFreq);
for idxFreq = 1:nFreq
    pairs = reshape(data(idxFreq, :), 2, []).';
    complexVals = convert_touchstone_pairs(pairs(:, 1), pairs(:, 2), dataFormat);
    S(:, :, idxFreq) = reshape(complexVals, nPorts, nPorts).';
end

network = struct( ...
    'file', filePath, ...
    'nPorts', nPorts, ...
    'freq', freq(:), ...
    'S', S, ...
    'z0', z0, ...
    'format', dataFormat);
end

function values = convert_touchstone_pairs(v1, v2, dataFormat)
switch upper(dataFormat)
    case 'MA'
        values = v1 .* exp(1i * deg2rad(v2));
    case 'DB'
        values = 10.^(v1 / 20) .* exp(1i * deg2rad(v2));
    case 'RI'
        values = complex(v1, v2);
    otherwise
        error('Unsupported Touchstone data format: %s', dataFormat);
end
end
