function parsed = parse_phase1_filename(filename)
%PARSE_PHASE1_FILENAME Support both agreed and current raw-data naming.

parsed = struct( ...
    'is_supported', false, ...
    'kind', '', ...
    'standard', '', ...
    'geometry', '', ...
    'port1', [], ...
    'port2', [], ...
    'band', '');

[~, stem, ext] = fileparts(filename);
if ~strcmpi(ext, '.s2p')
    return;
end

solExpr = '^(Short|Open|Load)_P(\d)_(0-67|67-115|110-170)$';
solMatch = regexp(stem, solExpr, 'tokens', 'once');
if ~isempty(solMatch)
    parsed.is_supported = true;
    parsed.kind = 'sol';
    parsed.standard = lower(solMatch{1});
    parsed.port1 = str2double(solMatch{2});
    parsed.band = solMatch{3};
    return;
end

thruExprA = '^THRU_P(\d)P(\d)_(STRAIGHT|ARC|DIAG|DIAGONAL)_(0-67|67-115|110-170)$';
thruMatchA = regexp(upper(stem), thruExprA, 'tokens', 'once');
if ~isempty(thruMatchA)
    parsed.is_supported = true;
    parsed.kind = 'thru';
    parsed.port1 = str2double(thruMatchA{1});
    parsed.port2 = str2double(thruMatchA{2});
    parsed.geometry = normalize_geometry(thruMatchA{3});
    parsed.band = thruMatchA{4};
    return;
end

thruExprB = '^Thru_(Straight|Arc|Diagonal|Diag)_P(\d)P(\d)_(0-67|67-115|110-170)$';
thruMatchB = regexp(stem, thruExprB, 'tokens', 'once');
if ~isempty(thruMatchB)
    parsed.is_supported = true;
    parsed.kind = 'thru';
    parsed.port1 = str2double(thruMatchB{2});
    parsed.port2 = str2double(thruMatchB{3});
    parsed.geometry = normalize_geometry(thruMatchB{1});
    parsed.band = thruMatchB{4};
end
end

function geometry = normalize_geometry(rawGeometry)
geometry = lower(rawGeometry);
if strcmp(geometry, 'diag')
    geometry = 'diagonal';
end
end
