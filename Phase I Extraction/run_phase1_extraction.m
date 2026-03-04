function results = run_phase1_extraction(rawDataDir)
%RUN_PHASE1_EXTRACTION Extract de-embedded Phase I standard definitions.
%   RESULTS = RUN_PHASE1_EXTRACTION() discovers raw measurements under
%   "Measured Data\Raw Standards", extracts one shared embedding network per
%   band from the straight thru, de-embeds available SOL and reciprocal thru
%   standards, and writes Touchstone and MAT outputs under:
%       Phase I Extraction\outputs
%
%   RESULTS = RUN_PHASE1_EXTRACTION(RAWDATADIR) uses a custom raw-data path.

if nargin < 1 || isempty(rawDataDir)
    candidateDirs = {
        fullfile(fileparts(mfilename('fullpath')), '..', 'Measured Data', 'Raw Standards')
        fullfile(fileparts(mfilename('fullpath')), '..', 'Measured Data', 'RawStandards')
    };
    rawDataDir = '';
    for idxDir = 1:numel(candidateDirs)
        if exist(candidateDirs{idxDir}, 'dir')
            rawDataDir = candidateDirs{idxDir};
            break;
        end
    end
end

scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir, 'outputs');
touchstoneDir = fullfile(outputDir, 'touchstone');
matDir = fullfile(outputDir, 'mat');

if isempty(rawDataDir) || ~exist(rawDataDir, 'dir')
    error('Raw data directory not found: %s', rawDataDir);
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
if ~exist(touchstoneDir, 'dir')
    mkdir(touchstoneDir);
end
if ~exist(matDir, 'dir')
    mkdir(matDir);
end

config = phase1_config(rawDataDir, touchstoneDir, matDir);
inventory = discover_phase1_measurements(config.raw_data_dir);

if isempty(inventory.bands)
    error('No supported Phase I measurement files were found under: %s', config.raw_data_dir);
end

results = struct('config', config, 'bands', {cell(1, numel(inventory.bands))});

for idxBand = 1:numel(inventory.bands)
    bandName = inventory.bands{idxBand};
    bandFiles = inventory.by_band.(inventory.band_keys{idxBand});
    fprintf('Processing band %s\n', bandName);
    results.bands{idxBand} = extract_phase1_band(config, bandName, bandFiles);
end

plot_phase1_stitched_results(config, results);

results.inventory = inventory;
save(fullfile(matDir, 'PHASE1_RESULTS_ALL.mat'), 'results');

fprintf('Phase I extraction complete.\n');
fprintf('Touchstone outputs: %s\n', touchstoneDir);
fprintf('MAT outputs: %s\n', matDir);
end
