function wave = parse_phase2_wave_csv(filePath)
%PARSE_PHASE2_WAVE_CSV Parse Keysight wave CSV with BEGIN CH1_DATA section.

if ~exist(filePath, 'file')
    error('Wave CSV not found: %s', filePath);
end

lines = string(splitlines(fileread(filePath)));
idxStart = find(strtrim(lines) == "BEGIN CH1_DATA", 1, 'first');
if isempty(idxStart) || idxStart >= numel(lines)
    error('BEGIN CH1_DATA block not found in %s', filePath);
end

headerLine = strtrim(lines(idxStart + 1));
headerParts = split(headerLine, ',');
if numel(headerParts) < 3
    error('Unexpected wave CSV header in %s: %s', filePath, headerLine);
end

dataLines = lines(idxStart + 2:end);
freq = [];
magDb = [];
phaseDeg = [];
for idx = 1:numel(dataLines)
    line = strtrim(dataLines(idx));
    if line == ""
        continue;
    end
    if isempty(regexp(char(line), '^[\+\-]?\d', 'once'))
        continue;
    end
    parts = textscan(char(line), '%f%f%f', 'Delimiter', ',');
    if isempty(parts{1})
        continue;
    end
    freq(end + 1, 1) = parts{1}; %#ok<AGROW>
    magDb(end + 1, 1) = parts{2}; %#ok<AGROW>
    phaseDeg(end + 1, 1) = parts{3}; %#ok<AGROW>
end

wave = struct();
wave.file = filePath;
wave.freq = freq;
wave.mag_db = magDb;
wave.phase_deg = phaseDeg;
wave.values = 10.^(magDb / 20) .* exp(1i * deg2rad(phaseDeg));
wave.header = headerLine;
end
