function sim = load_simulation_priors(scriptDir)
%LOAD_SIMULATION_PRIORS Load optional simulation priors for SOL refinement.

sim = struct();
sim.enabled = false;
sim.raw = struct();
sim.ref = struct();

baseDir = fullfile(scriptDir, '..', 'Simulated Data');
rawDir = fullfile(baseDir, 'Raw');
refDir = fullfile(baseDir, 'virtualmTRL_TUG');

if ~exist(rawDir, 'dir') || ~exist(refDir, 'dir')
    return;
end

rawMap = struct('short', 'short.s1p', 'open', 'open.s1p', 'load', 'load.s1p');
refMap = struct('short', 'Short-tugmtrl.s2p', 'open', 'Open-tugmtrl.s2p', 'load', 'Load-tugmtrl.s2p');

fields = {'short', 'open', 'load'};
for idx = 1:numel(fields)
    key = fields{idx};
    rawFile = fullfile(rawDir, rawMap.(key));
    refFile = fullfile(refDir, refMap.(key));
    if ~exist(rawFile, 'file') || ~exist(refFile, 'file')
        return;
    end

    rawData = read_oneport_touchstone_file(rawFile);
    refData = read_oneport_touchstone_file(refFile);

    sim.raw.(key) = struct('freq', rawData.freq, 'gamma', rawData.gamma);
    sim.ref.(key) = struct('freq', refData.freq, 'gamma', refData.gamma);
end

sim.enabled = true;
end