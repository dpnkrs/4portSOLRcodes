function inventory = discover_phase2_measurements(config)
%DISCOVER_PHASE2_MEASUREMENTS Discover Phase II standards, switch terms, and DUTs.

inventory = struct();
inventory.raw_standards = struct();
inventory.reference_thrus = struct();
inventory.switch_terms = struct();
inventory.duts = struct();
inventory.holdout_thrus = struct();
inventory.missing = {};

for idxBand = 1:numel(config.bands)
    bandKey = config.band_keys{idxBand};
    inventory.raw_standards.(bandKey) = struct();
    inventory.reference_thrus.(bandKey) = struct();
    inventory.switch_terms.(bandKey) = struct();
    inventory.duts.(bandKey) = struct();
    inventory.holdout_thrus.(bandKey) = struct();
end

rawFiles = dir(fullfile(config.phase2_raw_dir, '*.s4p'));
for idx = 1:numel(rawFiles)
    filePath = fullfile(rawFiles(idx).folder, rawFiles(idx).name);
    stdMatch = regexp(rawFiles(idx).name, '^(Open|Short|Load)_P([1-4])_(0-67|67-115|110-170)\.s4p$', 'tokens', 'once');
    if ~isempty(stdMatch)
        stdName = stdMatch{1};
        portName = sprintf('P%s', stdMatch{2});
        bandKey = local_band_key(stdMatch{3});
        if ~isfield(inventory.raw_standards.(bandKey), stdName)
            inventory.raw_standards.(bandKey).(stdName) = struct();
        end
        inventory.raw_standards.(bandKey).(stdName).(portName) = filePath;
        continue;
    end

    thruMeta = local_parse_thru_name(rawFiles(idx).name);
    if ~isempty(thruMeta)
        bandKey = local_band_key(thruMeta.band);
        pairKey = thruMeta.pair;
        geomKey = thruMeta.geometry;
        if ~isfield(inventory.reference_thrus.(bandKey), pairKey)
            inventory.reference_thrus.(bandKey).(pairKey) = struct();
        end
        inventory.reference_thrus.(bandKey).(pairKey).(geomKey) = filePath;
    end
end

switchFiles = dir(fullfile(config.phase2_switch_dir, '*.csv'));
for idx = 1:numel(switchFiles)
    meta = local_parse_switch_name(switchFiles(idx).name);
    if isempty(meta)
        continue;
    end
    bandKey = local_band_key(meta.band);
    if ~isfield(inventory.switch_terms.(bandKey), meta.pair)
        inventory.switch_terms.(bandKey).(meta.pair) = struct();
    end
    if ~isfield(inventory.switch_terms.(bandKey).(meta.pair), meta.geometry)
        inventory.switch_terms.(bandKey).(meta.pair).(meta.geometry) = struct();
    end
    inventory.switch_terms.(bandKey).(meta.pair).(meta.geometry).(meta.wave_field) = ...
        fullfile(switchFiles(idx).folder, switchFiles(idx).name);
end

dutFiles = dir(fullfile(config.phase2_dut_dir, '*.s4p'));
for idx = 1:numel(dutFiles)
    fileName = dutFiles(idx).name;
    filePath = fullfile(dutFiles(idx).folder, fileName);
    dutMatch = regexp(fileName, '^LangeCoupler_(0-67|67-115|110-170)\.s4p$', 'tokens', 'once');
    if ~isempty(dutMatch)
        bandKey = local_band_key(dutMatch{1});
        inventory.duts.(bandKey).LangeCoupler = filePath;
        continue;
    end

    thruMeta = local_parse_thru_name(fileName);
    if ~isempty(thruMeta)
        bandKey = local_band_key(thruMeta.band);
        if ~isfield(inventory.holdout_thrus.(bandKey), thruMeta.pair)
            inventory.holdout_thrus.(bandKey).(thruMeta.pair) = struct();
        end
        inventory.holdout_thrus.(bandKey).(thruMeta.pair).(thruMeta.geometry) = filePath;
    end
end

inventory.missing = local_validate_inventory(config, inventory);
if ~isempty(inventory.missing)
    error('Phase II data discovery failed:\n%s', strjoin(inventory.missing, newline));
end
end

function missing = local_validate_inventory(config, inventory)
missing = {};

for idxBand = 1:numel(config.bands)
    bandName = config.bands{idxBand};
    bandKey = config.band_keys{idxBand};
    rawStandards = inventory.raw_standards.(bandKey);
    for idxStd = 1:numel(config.standard_names)
        stdName = config.standard_names{idxStd};
        for idxPort = 1:numel(config.port_labels)
            portName = config.port_labels{idxPort};
            if ~isfield(rawStandards, stdName) || ~isfield(rawStandards.(stdName), portName)
                missing{end + 1} = sprintf('Missing %s %s raw standard for band %s', stdName, portName, bandName); %#ok<AGROW>
            end
        end
    end

    for idxPair = 1:numel(config.reference_pairs)
        pairKey = config.reference_pairs(idxPair).pair;
        geomKey = upper(config.reference_pairs(idxPair).geometry);
        if ~isfield(inventory.reference_thrus.(bandKey), pairKey) || ...
                ~isfield(inventory.reference_thrus.(bandKey).(pairKey), geomKey)
            missing{end + 1} = sprintf('Missing reference thru %s %s for band %s', pairKey, geomKey, bandName); %#ok<AGROW>
        end
        if ~isfield(inventory.switch_terms.(bandKey), pairKey) || ...
                ~isfield(inventory.switch_terms.(bandKey).(pairKey), geomKey)
            missing{end + 1} = sprintf('Missing switch terms for %s %s in band %s', pairKey, geomKey, bandName); %#ok<AGROW>
        else
            fields = {'aForward', 'bForward', 'aReverse', 'bReverse'};
            switchStruct = inventory.switch_terms.(bandKey).(pairKey).(geomKey);
            for idxField = 1:numel(fields)
                if ~isfield(switchStruct, fields{idxField})
                    missing{end + 1} = sprintf('Missing switch wave %s for %s %s in band %s', fields{idxField}, pairKey, geomKey, bandName); %#ok<AGROW>
                end
            end
        end
    end

    for idxHold = 1:numel(config.holdout_pairs)
        pairKey = config.holdout_pairs(idxHold).pair;
        geomKey = upper(config.holdout_pairs(idxHold).geometry);
        if ~isfield(inventory.holdout_thrus.(bandKey), pairKey) || ...
                ~isfield(inventory.holdout_thrus.(bandKey).(pairKey), geomKey)
            missing{end + 1} = sprintf('Missing hold-out thru %s %s for band %s', pairKey, geomKey, bandName); %#ok<AGROW>
        end
    end

    if ~isfield(inventory.duts.(bandKey), 'LangeCoupler')
        missing{end + 1} = sprintf('Missing LangeCoupler DUT for band %s', bandName); %#ok<AGROW>
    end
end
end

function meta = local_parse_thru_name(fileName)
meta = [];
tokens = regexp(fileName, '^Thru_(.+)_(Straight|Arc|Diag)_(0-67|67-115|110-170)\.s4p$', 'tokens', 'once');
if isempty(tokens)
    return;
end
pairRaw = upper(regexprep(tokens{1}, '[-_]', ''));
pairTokens = regexp(pairRaw, '^P([1-4])P([1-4])$', 'tokens', 'once');
if isempty(pairTokens)
    return;
end
meta = struct();
meta.pair = sprintf('P%sP%s', pairTokens{1}, pairTokens{2});
meta.geometry = upper(tokens{2});
meta.band = tokens{3};
end

function meta = local_parse_switch_name(fileName)
meta = [];
tokens = regexp(fileName, '^Thru_(.+)_(Straight|Arc|Diag)_(0-67|67-115|110-170)_([ab])([1-4])S([1-4])(forward|reverse)\.csv$', 'tokens', 'once');
if isempty(tokens)
    return;
end
pairRaw = upper(regexprep(tokens{1}, '[-_]', ''));
pairTokens = regexp(pairRaw, '^P([1-4])P([1-4])$', 'tokens', 'once');
if isempty(pairTokens)
    return;
end
waveKind = lower(tokens{4});
recvPort = tokens{5};
stimPort = tokens{6};
direction = lower(tokens{7});
if strcmp(direction, 'forward')
    suffix = 'Forward';
else
    suffix = 'Reverse';
end

meta = struct();
meta.pair = sprintf('P%sP%s', pairTokens{1}, pairTokens{2});
meta.geometry = upper(tokens{2});
meta.band = tokens{3};
meta.wave_field = sprintf('%s%s', waveKind, suffix);
meta.receiver = sprintf('P%s', recvPort);
meta.stimulus = sprintf('P%s', stimPort);
end

function bandKey = local_band_key(bandName)
bandKey = matlab.lang.makeValidName(strrep(bandName, '-', '_'));
end
