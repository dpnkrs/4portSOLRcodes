function inventory = discover_phase1_measurements(rawDataDir)
%DISCOVER_PHASE1_MEASUREMENTS Build a per-band inventory from raw filenames.

files = dir(fullfile(rawDataDir, '*.s2p'));

inventory = struct();
inventory.files = files;
inventory.bands = {};
inventory.band_keys = {};
inventory.by_band = struct();
inventory.unparsed = {};

for idxFile = 1:numel(files)
    file = files(idxFile);
    parsed = parse_phase1_filename(file.name);

    if ~parsed.is_supported
        inventory.unparsed{end + 1} = file.name; %#ok<AGROW>
        continue;
    end

    bandKey = matlab.lang.makeValidName(['BAND_' strrep(parsed.band, '-', '_')]);
    if ~isfield(inventory.by_band, bandKey)
        inventory.by_band.(bandKey) = struct();
        inventory.bands{end + 1} = parsed.band; %#ok<AGROW>
        inventory.band_keys{end + 1} = bandKey; %#ok<AGROW>
    end

    fullPath = fullfile(file.folder, file.name);

    if strcmp(parsed.kind, 'sol')
        solKey = sprintf('%s_P%d', upper(parsed.standard), parsed.port1);
        inventory.by_band.(bandKey).(solKey) = fullPath;
    else
        pairLabel = sprintf('P%dP%d', parsed.port1, parsed.port2);
        geoLabel = upper(parsed.geometry);
        thruKey = sprintf('THRU_%s_%s', pairLabel, geoLabel);
        inventory.by_band.(bandKey).(thruKey) = fullPath;
    end
end
end
