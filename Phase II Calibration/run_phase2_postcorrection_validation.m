function summary = run_phase2_postcorrection_validation(config, resultsPath, outputBaseDir)
%RUN_PHASE2_POSTCORRECTION_VALIDATION Generate diagnostic-only physical validation outputs.

if nargin < 1 || isempty(config)
    config = phase2_config();
end
if nargin < 2 || isempty(resultsPath)
    resultsPath = fullfile(config.mat_dir, 'PHASE2_RESULTS_ALL.mat');
end
if nargin < 3
    outputBaseDir = '';
end

results = load_results_container(resultsPath);

validationDirs = ensure_validation_dirs(config, outputBaseDir);
referenceStandards = results.references;
summary = struct();
summary.sol = struct();
summary.reciprocal = struct();

for idxPort = 1:numel(config.port_labels)
    portLabel = config.port_labels{idxPort};
    portData = collect_port_validation_data(results, config, portLabel);
    fig = create_port_validation_figure(config, portLabel, portData, referenceStandards);
    save_validation_figure(fig, fullfile(validationDirs.figure_dir, sprintf('PhaseII_PostCorrection_SOL_%s', portLabel)));
    summary.sol.(portLabel) = summarize_port_validation(portData, referenceStandards, config);
end

for idxPair = 1:numel(config.reference_pairs)
    pairKey = upper(config.reference_pairs(idxPair).pair);
    geomKey = upper(config.reference_pairs(idxPair).geometry);
    pairData = collect_pair_validation_data(results, config, pairKey, geomKey);
    fig = create_pair_validation_figure(config, pairData);
    save_validation_figure(fig, fullfile(validationDirs.figure_dir, sprintf('PhaseII_PostCorrection_Thru_%s_%s', pairKey, geomKey)));
    summary.reciprocal.(pairKey).(geomKey) = summarize_pair_validation(pairData, config);
end

summary.results_path = resultsPath;
summary.validation_dirs = validationDirs;
summary.generated_at = char(datetime('now'));
save(fullfile(validationDirs.mat_dir, 'PHASE2_POSTCORRECTION_VALIDATION.mat'), 'summary');

noteLines = build_validation_note(summary, config);
write_text_file(fullfile(validationDirs.note_dir, 'PhaseII_PostCorrection_ValidationSummary.txt'), noteLines);
end

function dirs = ensure_validation_dirs(config, outputBaseDir)
if nargin < 2 || isempty(outputBaseDir)
    baseDir = fullfile(config.interim_dir, 'postCorrectionValidation');
else
    baseDir = fullfile(outputBaseDir, 'postCorrectionValidation');
end
dirs = struct();
dirs.base_dir = baseDir;
dirs.figure_dir = fullfile(baseDir, 'figures');
dirs.note_dir = fullfile(baseDir, 'notes');
dirs.mat_dir = fullfile(baseDir, 'mat');

ensure_dir(baseDir);
ensure_dir(dirs.figure_dir);
ensure_dir(dirs.note_dir);
ensure_dir(dirs.mat_dir);
end

function results = load_results_container(resultsPath)
loaded = load(resultsPath);
if isfield(loaded, 'results')
    results = loaded.results;
elseif isfield(loaded, 'experimentResults')
    results = loaded.experimentResults;
elseif isfield(loaded, 'bundle') && isfield(loaded.bundle, 'experiment')
    results = loaded.bundle.experiment;
elseif isfield(loaded, 'bundle') && isfield(loaded.bundle, 'baseline')
    results = loaded.bundle.baseline;
else
    error('PhaseII:ResultsContainerNotFound', ...
        'Expected results, experimentResults, or bundle.experiment in %s', resultsPath);
end
end

function portData = collect_port_validation_data(results, config, portLabel)
portData = struct();
portData.port = portLabel;
portData.band_entries = struct();
cache = struct();

for idxStd = 1:numel(config.standard_names)
    stdName = config.standard_names{idxStd};
    portData.band_entries.(stdName) = struct();
end

for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    bandName = bandResult.band;
    bandKey = matlab.lang.makeValidName(strrep(bandName, '-', '_'));

    for idxStd = 1:numel(config.standard_names)
        stdName = config.standard_names{idxStd};
        entries = struct([]);
        for idxPair = 1:numel(bandResult.reference_pair_results)
            pairResult = bandResult.reference_pair_results(idxPair);
            ports = pairResult.ports;
            localLabels = {sprintf('P%d', ports(1)), sprintf('P%d', ports(2))};
            sideIndex = find(strcmpi(localLabels, portLabel), 1, 'first');
            if isempty(sideIndex)
                continue;
            end

            rawPath = results.inventory.raw_standards.(bandKey).(stdName).(portLabel);
            [rawTrace, cache] = load_raw_standard_trace(rawPath, pairResult.pair, sideIndex, cache);
            if sideIndex == 1
                switchTrace = pairResult.standards_measured_port1.(stdName);
            else
                switchTrace = pairResult.standards_measured_port2.(stdName);
            end
            correctedTrace = squeeze(pairResult.standards_corrected.(stdName)(sideIndex, :)).';

            entry = struct( ...
                'band', bandName, ...
                'band_key', bandKey, ...
                'freq', pairResult.freq(:), ...
                'pair', pairResult.pair, ...
                'geometry', pairResult.geometry, ...
                'pair_label', sprintf('%s %s', pairResult.pair, pairResult.geometry), ...
                'side_index', sideIndex, ...
                'raw', rawTrace(:), ...
                'switch', switchTrace(:), ...
                'corrected', correctedTrace(:));
            if isempty(entries)
                entries = entry;
            else
                entries(end + 1) = entry; %#ok<AGROW>
            end
        end
        portData.band_entries.(stdName).(bandKey) = entries;
    end
end
end

function [trace, cache] = load_raw_standard_trace(rawPath, pairKey, sideIndex, cache)
cacheKey = matlab.lang.makeValidName(strrep(rawPath, filesep, '_'));
if isfield(cache, cacheKey)
    net = cache.(cacheKey);
else
    net = read_touchstone_nport(rawPath);
    cache.(cacheKey) = net;
end
Sraw = phase2_extract_pair_submatrix(net.S, pairKey);
trace = squeeze(Sraw(sideIndex, sideIndex, :));
end

function fig = create_port_validation_figure(config, portLabel, portData, referenceStandards)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1700, 1150]);
tlo = tiledlayout(fig, 3, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
stdOrder = {'Open', 'Short', 'Load'};
pairColors = lines(6);

for idxStd = 1:numel(stdOrder)
    stdName = stdOrder{idxStd};
    [refFreq, refGamma] = flatten_reference_trace(referenceStandards.(stdName), config.bands);
    [idealGamma, idealMagDb, idealPhaseDeg] = ideal_standard_target(stdName, refFreq);

    axMag = nexttile(tlo);
    hold(axMag, 'on');
    style_standard_axes(axMag, 'Frequency (GHz)', sprintf('%s |\\Gamma| (dB)', stdName));
    title(axMag, sprintf('%s magnitude', stdName));
    plot(axMag, refFreq / 1e9, idealMagDb, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s ideal', stdName));
    plot(axMag, refFreq / 1e9, safe_mag_db(refGamma), '-', 'Color', [0.1 0.5 0.1], 'LineWidth', 2.0, ...
        'DisplayName', sprintf('%s Phase I reference', stdName));

    axPhase = nexttile(tlo);
    hold(axPhase, 'on');
    style_standard_axes(axPhase, 'Frequency (GHz)', sprintf('%s phase (deg)', stdName));
    title(axPhase, sprintf('%s phase', stdName));
    plot(axPhase, refFreq / 1e9, idealPhaseDeg, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%s ideal', stdName));
    plot(axPhase, refFreq / 1e9, align_phase_deg(refGamma, idealGamma), '-', 'Color', [0.1 0.5 0.1], 'LineWidth', 2.0, ...
        'DisplayName', sprintf('%s Phase I reference', stdName));

    axSmith = nexttile(tlo);
    hold(axSmith, 'on');
    setup_smith_axes(axSmith);
    title(axSmith, sprintf('%s Smith chart', stdName));
    plot(axSmith, real(refGamma), imag(refGamma), '-', 'Color', [0.1 0.5 0.1], 'LineWidth', 2.0, ...
        'DisplayName', sprintf('%s Phase I reference', stdName));
    plot(axSmith, real(idealGamma), imag(idealGamma), 'o', 'Color', [0.55 0.55 0.55], ...
        'MarkerFaceColor', [0.55 0.55 0.55], 'DisplayName', sprintf('%s ideal', stdName));

    legendHandles = gobjects(0);
    legendLabels = {};
    for idxBand = 1:numel(config.bands)
        bandKey = matlab.lang.makeValidName(strrep(config.bands{idxBand}, '-', '_'));
        entries = portData.band_entries.(stdName).(bandKey);
        if isempty(entries)
            continue;
        end

        [rawAgg, switchAgg, correctedAgg, freq] = aggregate_standard_entries(entries, idealGamma);

        hRaw = plot(axMag, freq / 1e9, safe_mag_db(rawAgg), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.3);
        plot(axPhase, freq / 1e9, align_phase_deg(rawAgg, idealGamma), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.3);
        plot(axSmith, real(rawAgg), imag(rawAgg), ':', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.3);

        hSwitch = plot(axMag, freq / 1e9, safe_mag_db(switchAgg), '--', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5);
        plot(axPhase, freq / 1e9, align_phase_deg(switchAgg, idealGamma), '--', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5);
        plot(axSmith, real(switchAgg), imag(switchAgg), '--', 'Color', [0.00 0.45 0.74], 'LineWidth', 1.5);

        for idxEntry = 1:numel(entries)
            entryColor = pairColors(mod(idxEntry - 1, size(pairColors, 1)) + 1, :);
            plot(axMag, freq / 1e9, safe_mag_db(entries(idxEntry).corrected), '-', 'Color', entryColor, 'LineWidth', 1.0);
            plot(axPhase, freq / 1e9, align_phase_deg(entries(idxEntry).corrected, idealGamma), '-', 'Color', entryColor, 'LineWidth', 1.0);
            plot(axSmith, real(entries(idxEntry).corrected), imag(entries(idxEntry).corrected), '-', 'Color', entryColor, 'LineWidth', 1.0);
        end

        hCorr = plot(axMag, freq / 1e9, safe_mag_db(correctedAgg), '-', 'Color', [0 0 0], 'LineWidth', 2.1);
        plot(axPhase, freq / 1e9, align_phase_deg(correctedAgg, idealGamma), '-', 'Color', [0 0 0], 'LineWidth', 2.1);
        plot(axSmith, real(correctedAgg), imag(correctedAgg), '-', 'Color', [0 0 0], 'LineWidth', 2.1);

        if isempty(legendHandles)
            legendHandles = [hRaw; hSwitch; hCorr]; %#ok<AGROW>
            legendLabels = {'Raw aggregate', 'Switch-corrected aggregate', 'Corrected aggregate'}; %#ok<AGROW>
            for idxEntry = 1:numel(entries)
                legendHandles(end + 1) = plot(axMag, nan, nan, '-', 'Color', pairColors(mod(idxEntry - 1, size(pairColors, 1)) + 1, :), 'LineWidth', 1.0); %#ok<AGROW>
                legendLabels{end + 1} = entries(idxEntry).pair_label; %#ok<AGROW>
            end
            legendHandles(end + 1) = plot(axMag, nan, nan, '-', 'Color', [0.1 0.5 0.1], 'LineWidth', 2.0); %#ok<AGROW>
            legendLabels{end + 1} = 'Phase I reference'; %#ok<AGROW>
            legendHandles(end + 1) = plot(axMag, nan, nan, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5); %#ok<AGROW>
            legendLabels{end + 1} = 'Ideal target'; %#ok<AGROW>
        end
    end

    if ~isempty(legendHandles)
        lgd = legend(axMag, legendHandles, legendLabels, 'Location', 'eastoutside');
        set(lgd, 'Interpreter', 'none');
    end
end

title(tlo, sprintf('Phase II post-correction SOL physical validation - %s', portLabel), ...
    'FontWeight', 'bold', 'FontSize', 16);
end

function [rawAgg, switchAgg, correctedAgg, freq] = aggregate_standard_entries(entries, idealGamma)
freq = entries(1).freq(:);
rawMatrix = nan(numel(freq), numel(entries));
switchMatrix = nan(numel(freq), numel(entries));
correctedMatrix = nan(numel(freq), numel(entries));
for idxEntry = 1:numel(entries)
    rawMatrix(:, idxEntry) = entries(idxEntry).raw(:);
    switchMatrix(:, idxEntry) = entries(idxEntry).switch(:);
    correctedMatrix(:, idxEntry) = entries(idxEntry).corrected(:);
end
rawAgg = aggregate_complex_traces(rawMatrix, idealGamma);
switchAgg = aggregate_complex_traces(switchMatrix, idealGamma);
correctedAgg = aggregate_complex_traces(correctedMatrix, idealGamma);
end

function aggregateTrace = aggregate_complex_traces(traceMatrix, targetGamma)
magMedian = median(abs(traceMatrix), 2, 'omitnan');
phaseMatrix = nan(size(traceMatrix));
for idxTrace = 1:size(traceMatrix, 2)
    phaseMatrix(:, idxTrace) = align_phase_rad(traceMatrix(:, idxTrace), targetGamma);
end
phaseMedian = median(phaseMatrix, 2, 'omitnan');
aggregateTrace = magMedian .* exp(1j * phaseMedian);
end

function summary = summarize_port_validation(portData, referenceStandards, config)
summary = struct();
stdOrder = {'Open', 'Short', 'Load'};
for idxStd = 1:numel(stdOrder)
    stdName = stdOrder{idxStd};
    summary.(stdName) = struct();
    for idxBand = 1:numel(config.bands)
        bandName = config.bands{idxBand};
        bandKey = matlab.lang.makeValidName(strrep(bandName, '-', '_'));
        entries = portData.band_entries.(stdName).(bandKey);
        if isempty(entries)
            continue;
        end
        refTrace = interpolate_reference_gamma(referenceStandards.(stdName), entries(1).freq, bandName);
        [idealGamma, ~, ~] = ideal_standard_target(stdName, entries(1).freq);
        [~, ~, correctedAgg, ~] = aggregate_standard_entries(entries, idealGamma);
        [failFlag, metricValue, metricLabel] = summarize_standard_metric(stdName, correctedAgg, idealGamma);
        summary.(stdName).(bandKey) = struct( ...
            'band', bandName, ...
            'median_mag_db', median(safe_mag_db(correctedAgg), 'omitnan'), ...
            'median_phase_deg', median(align_phase_deg(correctedAgg, idealGamma), 'omitnan'), ...
            'median_ref_deviation', median(abs(correctedAgg - refTrace), 'omitnan'), ...
            'metric_value', metricValue, ...
            'metric_label', metricLabel, ...
            'fail_flag', failFlag, ...
            'num_pair_traces', numel(entries));
    end
end
end

function [failFlag, metricValue, metricLabel] = summarize_standard_metric(stdName, correctedGamma, idealGamma)
switch lower(stdName)
    case 'open'
        metricValue = median(abs(correctedGamma - idealGamma), 'omitnan');
        metricLabel = 'median |Gamma-(+1)|';
        failFlag = metricValue > 0.25;
    case 'short'
        metricValue = median(abs(correctedGamma - idealGamma), 'omitnan');
        metricLabel = 'median |Gamma-(-1)|';
        failFlag = metricValue > 0.25;
    case 'load'
        metricValue = median(abs(correctedGamma), 'omitnan');
        metricLabel = 'median |Gamma|';
        failFlag = metricValue > 0.25;
    otherwise
        error('Unsupported standard %s', stdName);
end
end

function pairData = collect_pair_validation_data(results, config, pairKey, geomKey)
pairData = struct();
pairData.pair = pairKey;
pairData.geometry = geomKey;
pairData.band_entries = struct();

for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    hit = find(strcmpi({bandResult.reference_pair_results.pair}, pairKey) & ...
        strcmpi({bandResult.reference_pair_results.geometry}, geomKey), 1, 'first');
    if isempty(hit)
        continue;
    end
    pairResult = bandResult.reference_pair_results(hit);
    bandKey = matlab.lang.makeValidName(strrep(bandResult.band, '-', '_'));

    phase1Ref = load_phase1_thru_target_local(config, pairKey, geomKey, bandResult.band);
    pairData.band_entries.(bandKey) = struct( ...
        'band', bandResult.band, ...
        'freq', pairResult.freq(:), ...
        'raw', pairResult.reference_raw, ...
        'switch', pairResult.reference_switch_corrected, ...
        'corrected', pairResult.reference_corrected, ...
        'phase1', phase1Ref);
end
end

function fig = create_pair_validation_figure(config, pairData)
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1750, 1000]);
tlo = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

axMag = nexttile(tlo);
hold(axMag, 'on');
style_standard_axes(axMag, 'Frequency (GHz)', 'Transmission magnitude (dB)');
title(axMag, sprintf('%s %s |S21| and |S12|', pairData.pair, pairData.geometry));
yline(axMag, 0, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.3, 'DisplayName', '0 dB upper guide');

axPhase = nexttile(tlo);
hold(axPhase, 'on');
style_standard_axes(axPhase, 'Frequency (GHz)', 'Transmission phase (deg)');
title(axPhase, sprintf('%s %s transmission phase', pairData.pair, pairData.geometry));

axS11 = nexttile(tlo);
hold(axS11, 'on');
setup_smith_axes(axS11);
title(axS11, sprintf('%s %s Smith chart S11', pairData.pair, pairData.geometry));
plot(axS11, 0, 0, 'o', 'Color', [0.55 0.55 0.55], 'MarkerFaceColor', [0.55 0.55 0.55], 'DisplayName', 'Ideal low reflection');

axS22 = nexttile(tlo);
hold(axS22, 'on');
setup_smith_axes(axS22);
title(axS22, sprintf('%s %s Smith chart S22', pairData.pair, pairData.geometry));
plot(axS22, 0, 0, 'o', 'Color', [0.55 0.55 0.55], 'MarkerFaceColor', [0.55 0.55 0.55], 'DisplayName', 'Ideal low reflection');

legendHandles = gobjects(0);
legendLabels = {};
bandColors = lines(numel(config.bands));

for idxBand = 1:numel(config.bands)
    bandKey = matlab.lang.makeValidName(strrep(config.bands{idxBand}, '-', '_'));
    if ~isfield(pairData.band_entries, bandKey)
        continue;
    end
    entry = pairData.band_entries.(bandKey);
    freqGHz = entry.freq / 1e9;
    bandColor = bandColors(idxBand, :);
    colorRaw = tint_color(bandColor, 0.65);
    colorSwitch = tint_color(bandColor, 0.25);
    colorCorr = bandColor;

    hRaw21 = plot(axMag, freqGHz, safe_mag_db(squeeze(entry.raw(2, 1, :))), ':', 'Color', colorRaw, 'LineWidth', 1.1);
    plot(axMag, freqGHz, safe_mag_db(squeeze(entry.raw(1, 2, :))), ':', 'Color', colorRaw, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    hSwitch21 = plot(axMag, freqGHz, safe_mag_db(squeeze(entry.switch(2, 1, :))), '--', 'Color', colorSwitch, 'LineWidth', 1.3);
    plot(axMag, freqGHz, safe_mag_db(squeeze(entry.switch(1, 2, :))), '--', 'Color', colorSwitch, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    hCorr21 = plot(axMag, freqGHz, safe_mag_db(squeeze(entry.corrected(2, 1, :))), '-', 'Color', colorCorr, 'LineWidth', 1.8);
    plot(axMag, freqGHz, safe_mag_db(squeeze(entry.corrected(1, 2, :))), '-', 'Color', colorCorr, 'LineWidth', 1.8, 'HandleVisibility', 'off');

    phaseTarget = squeeze(entry.corrected(2, 1, :));
    plot(axPhase, freqGHz, align_phase_deg(squeeze(entry.raw(2, 1, :)), phaseTarget), ':', 'Color', colorRaw, 'LineWidth', 1.1);
    plot(axPhase, freqGHz, align_phase_deg(squeeze(entry.raw(1, 2, :)), phaseTarget), ':', 'Color', colorRaw, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    plot(axPhase, freqGHz, align_phase_deg(squeeze(entry.switch(2, 1, :)), phaseTarget), '--', 'Color', colorSwitch, 'LineWidth', 1.3);
    plot(axPhase, freqGHz, align_phase_deg(squeeze(entry.switch(1, 2, :)), phaseTarget), '--', 'Color', colorSwitch, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    plot(axPhase, freqGHz, align_phase_deg(squeeze(entry.corrected(2, 1, :)), phaseTarget), '-', 'Color', colorCorr, 'LineWidth', 1.8);
    plot(axPhase, freqGHz, align_phase_deg(squeeze(entry.corrected(1, 2, :)), phaseTarget), '-', 'Color', colorCorr, 'LineWidth', 1.8, 'HandleVisibility', 'off');

    plot(axS11, real(squeeze(entry.raw(1, 1, :))), imag(squeeze(entry.raw(1, 1, :))), ':', 'Color', colorRaw, 'LineWidth', 1.1);
    plot(axS11, real(squeeze(entry.switch(1, 1, :))), imag(squeeze(entry.switch(1, 1, :))), '--', 'Color', colorSwitch, 'LineWidth', 1.3);
    plot(axS11, real(squeeze(entry.corrected(1, 1, :))), imag(squeeze(entry.corrected(1, 1, :))), '-', 'Color', colorCorr, 'LineWidth', 1.8);

    plot(axS22, real(squeeze(entry.raw(2, 2, :))), imag(squeeze(entry.raw(2, 2, :))), ':', 'Color', colorRaw, 'LineWidth', 1.1);
    plot(axS22, real(squeeze(entry.switch(2, 2, :))), imag(squeeze(entry.switch(2, 2, :))), '--', 'Color', colorSwitch, 'LineWidth', 1.3);
    plot(axS22, real(squeeze(entry.corrected(2, 2, :))), imag(squeeze(entry.corrected(2, 2, :))), '-', 'Color', colorCorr, 'LineWidth', 1.8);

    if ~isempty(entry.phase1)
        plot(axMag, freqGHz, safe_mag_db(squeeze(entry.phase1.S(2, 1, :))), '-.', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.6);
        plot(axMag, freqGHz, safe_mag_db(squeeze(entry.phase1.S(1, 2, :))), '-.', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.6, 'HandleVisibility', 'off');
        phaseRefTarget = squeeze(entry.phase1.S(2, 1, :));
        plot(axPhase, freqGHz, align_phase_deg(squeeze(entry.phase1.S(2, 1, :)), phaseRefTarget), '-.', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.6);
        plot(axPhase, freqGHz, align_phase_deg(squeeze(entry.phase1.S(1, 2, :)), phaseRefTarget), '-.', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.6, 'HandleVisibility', 'off');
        plot(axS11, real(squeeze(entry.phase1.S(1, 1, :))), imag(squeeze(entry.phase1.S(1, 1, :))), '-.', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.6);
        plot(axS22, real(squeeze(entry.phase1.S(2, 2, :))), imag(squeeze(entry.phase1.S(2, 2, :))), '-.', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.6);
    end

    legendHandles(end + 1) = hRaw21; %#ok<AGROW>
    legendLabels{end + 1} = sprintf('Raw (%s)', entry.band); %#ok<AGROW>
    legendHandles(end + 1) = hSwitch21; %#ok<AGROW>
    legendLabels{end + 1} = sprintf('Switch-corrected (%s)', entry.band); %#ok<AGROW>
    legendHandles(end + 1) = hCorr21; %#ok<AGROW>
    legendLabels{end + 1} = sprintf('Corrected (%s)', entry.band); %#ok<AGROW>
end

legendHandles(end + 1) = plot(axMag, nan, nan, '-.', 'Color', [0.1 0.5 0.1], 'LineWidth', 1.6); %#ok<AGROW>
legendLabels{end + 1} = 'Phase I reference'; %#ok<AGROW>
legendHandles(end + 1) = plot(axMag, nan, nan, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.3); %#ok<AGROW>
legendLabels{end + 1} = '0 dB / ideal low-reflection guide'; %#ok<AGROW>
lgd = legend(axMag, legendHandles, legendLabels, 'Location', 'eastoutside');
set(lgd, 'Interpreter', 'none');

title(tlo, sprintf('Phase II post-correction reciprocal-thru physical validation - %s %s', pairData.pair, pairData.geometry), ...
    'FontWeight', 'bold', 'FontSize', 16);
end

function summary = summarize_pair_validation(pairData, config)
summary = struct();
for idxBand = 1:numel(config.bands)
    bandKey = matlab.lang.makeValidName(strrep(config.bands{idxBand}, '-', '_'));
    if ~isfield(pairData.band_entries, bandKey)
        continue;
    end
    entry = pairData.band_entries.(bandKey);
    s21 = squeeze(entry.corrected(2, 1, :));
    s12 = squeeze(entry.corrected(1, 2, :));
    s11 = squeeze(entry.corrected(1, 1, :));
    s22 = squeeze(entry.corrected(2, 2, :));

    phase21 = unwrap(angle(s21)) * 180 / pi;
    phase12 = unwrap(angle(s12)) * 180 / pi;
    [viol21, frac21] = phase_monotonicity_violations(phase21);
    [viol12, frac12] = phase_monotonicity_violations(phase12);

    mag21 = safe_mag_db(s21);
    mag12 = safe_mag_db(s12);
    aboveZero21 = nnz(mag21 > 0);
    aboveZero12 = nnz(mag12 > 0);
    medianReturnRadius = median(max(abs([s11, s22]), [], 2), 'omitnan');

    summary.(bandKey) = struct( ...
        'band', entry.band, ...
        'median_s21_db', median(mag21, 'omitnan'), ...
        'median_s12_db', median(mag12, 'omitnan'), ...
        'median_s11_radius', median(abs(s11), 'omitnan'), ...
        'median_s22_radius', median(abs(s22), 'omitnan'), ...
        'median_return_radius', medianReturnRadius, ...
        's21_above_zero_count', aboveZero21, ...
        's12_above_zero_count', aboveZero12, ...
        's21_max_excess_db', max([mag21; 0]), ...
        's12_max_excess_db', max([mag12; 0]), ...
        'phase21_violations', viol21, ...
        'phase12_violations', viol12, ...
        'phase21_violation_fraction', frac21, ...
        'phase12_violation_fraction', frac12, ...
        'fail_gain', aboveZero21 > 0 || aboveZero12 > 0, ...
        'fail_phase', frac21 > 0.05 || frac12 > 0.05, ...
        'fail_return_loss', medianReturnRadius > 0.316);
end
end

function lines = build_validation_note(summary, config)
lines = {
    'Phase II Post-Correction Physical Validation Summary'
    '==================================================='
    ''
    'Diagnostic basis: raw downstream traces, switch-corrected traces, corrected traces, Phase I adopted references, and ideal textbook guides.'
    ''
    'Decision rule: do not integrate any mitigation unless the corrected physical-response plots improve, not only the target-mismatch summaries.'
    ''
    };

for idxPort = 1:numel(config.port_labels)
    portLabel = config.port_labels{idxPort};
    lines{end + 1} = sprintf('Port %s', portLabel); %#ok<AGROW>
    for idxStd = 1:numel(config.standard_names)
        stdName = config.standard_names{idxStd};
        stdSummary = summary.sol.(portLabel).(stdName);
        bandKeys = fieldnames(stdSummary);
        for idxBand = 1:numel(bandKeys)
            bandKey = bandKeys{idxBand};
            item = stdSummary.(bandKey);
            failText = ternary_text(item.fail_flag, 'FAIL', 'PASS');
            lines{end + 1} = sprintf('  %s %s: %s = %.4f, median ref deviation = %.4f, pair traces = %d [%s]', ... %#ok<AGROW>
                item.band, stdName, item.metric_label, item.metric_value, item.median_ref_deviation, item.num_pair_traces, failText);
        end
    end
    lines{end + 1} = ''; %#ok<AGROW>
end

pairNames = fieldnames(summary.reciprocal);
for idxPair = 1:numel(pairNames)
    pairKey = pairNames{idxPair};
    geomNames = fieldnames(summary.reciprocal.(pairKey));
    for idxGeom = 1:numel(geomNames)
        geomKey = geomNames{idxGeom};
        lines{end + 1} = sprintf('Reciprocal %s %s', pairKey, geomKey); %#ok<AGROW>
        bandKeys = fieldnames(summary.reciprocal.(pairKey).(geomKey));
        for idxBand = 1:numel(bandKeys)
            item = summary.reciprocal.(pairKey).(geomKey).(bandKeys{idxBand});
            failFlags = {};
            if item.fail_gain
                failFlags{end + 1} = '|S21| or |S12| exceeds 0 dB'; %#ok<AGROW>
            end
            if item.fail_phase
                failFlags{end + 1} = 'phase not monotonic'; %#ok<AGROW>
            end
            if item.fail_return_loss
                failFlags{end + 1} = 'return loss not plausibly low'; %#ok<AGROW>
            end
            if isempty(failFlags)
                failText = 'PASS';
            else
                failText = strjoin(failFlags, '; ');
            end
            lines{end + 1} = sprintf(['  %s: median |S21|/|S12| = %.3f / %.3f dB, ' ...
                'S21>0 count = %d, S12>0 count = %d, phase violations = %d / %d, ' ...
                'median |S11|/|S22| = %.3f / %.3f [%s]'], ... %#ok<AGROW>
                item.band, item.median_s21_db, item.median_s12_db, ...
                item.s21_above_zero_count, item.s12_above_zero_count, ...
                item.phase21_violations, item.phase12_violations, ...
                item.median_s11_radius, item.median_s22_radius, failText);
        end
        lines{end + 1} = ''; %#ok<AGROW>
    end
end
end

function [violations, fraction] = phase_monotonicity_violations(phaseDeg)
diffPhase = diff(phaseDeg(:));
trend = sign(median(diffPhase, 'omitnan'));
if trend == 0 || ~isfinite(trend)
    violations = 0;
    fraction = 0;
    return;
end
tolerance = 0.5;
violations = nnz((trend .* diffPhase) < -tolerance);
fraction = violations / max(numel(diffPhase), 1);
end

function [idealGamma, idealMagDb, idealPhaseDeg] = ideal_standard_target(stdName, freq)
freq = freq(:);
switch lower(stdName)
    case 'open'
        idealGamma = ones(size(freq));
        idealPhaseDeg = zeros(size(freq));
    case 'short'
        idealGamma = -ones(size(freq));
        idealPhaseDeg = 180 * ones(size(freq));
    case 'load'
        idealGamma = zeros(size(freq));
        idealPhaseDeg = zeros(size(freq));
    otherwise
        error('Unsupported standard %s', stdName);
end
idealMagDb = safe_mag_db(idealGamma);
idealMagDb(~isfinite(idealMagDb)) = -80;
end

function refGamma = interpolate_reference_gamma(reference, freq, bandName)
if isfield(reference, 'bands')
    bandKey = matlab.lang.makeValidName(strrep(bandName, '-', '_'));
    source = reference.bands.(bandKey);
else
    source = reference;
end
refGamma = interp1(source.freq(:), source.gamma(:), freq(:), 'pchip', 'extrap');
end

function [freq, gamma] = flatten_reference_trace(reference, bands)
if isfield(reference, 'bands')
    freq = [];
    gamma = [];
    for idxBand = 1:numel(bands)
        bandKey = matlab.lang.makeValidName(strrep(bands{idxBand}, '-', '_'));
        source = reference.bands.(bandKey);
        freq = [freq; source.freq(:); nan]; %#ok<AGROW>
        gamma = [gamma; source.gamma(:); nan]; %#ok<AGROW>
    end
else
    freq = reference.freq(:);
    gamma = reference.gamma(:);
end
end

function target = load_phase1_thru_target_local(config, pairKey, geomKey, bandName)
switch upper(geomKey)
    case 'STRAIGHT'
        if strcmpi(pairKey, 'P3P4')
            baseName = sprintf('EXTRACTED_THRU_P3P4_STRAIGHT_%s.mat', bandName);
        else
            baseName = sprintf('EXTRACTED_THRU_P1P2_STRAIGHT_%s.mat', bandName);
        end
    case 'ARC'
        baseName = sprintf('EXTRACTED_THRU_P1P4_ARC_%s.mat', bandName);
    case {'DIAG', 'DIAGONAL'}
        baseName = sprintf('EXTRACTED_THRU_P1P4_DIAGONAL_%s.mat', bandName);
    otherwise
        target = [];
        return;
end

fullPath = fullfile(config.phase1_output_mat_dir, baseName);
if ~exist(fullPath, 'file')
    target = [];
    return;
end
loaded = load(fullPath);
target = loaded.reciprocalResult;
end

function phaseDeg = align_phase_deg(trace, target)
phaseDeg = align_phase_rad(trace, target) * 180 / pi;
end

function phaseRad = align_phase_rad(trace, target)
trace = trace(:);
if isscalar(target)
    target = repmat(target, size(trace));
else
    target = target(:);
    if numel(target) ~= numel(trace)
        target = interp1(linspace(0, 1, numel(target)), target, linspace(0, 1, numel(trace)), 'linear', 'extrap').';
    end
end
phaseRad = unwrap(angle(trace));
targetPhase = unwrap(angle(target));
offsetCycles = round(median((targetPhase - phaseRad) / (2 * pi), 'omitnan'));
if ~isfinite(offsetCycles)
    offsetCycles = 0;
end
phaseRad = phaseRad + 2 * pi * offsetCycles;
end

function values = safe_mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function style_standard_axes(ax, xLabelText, yLabelText)
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, xLabelText);
ylabel(ax, yLabelText);
xline(ax, 67, ':', 'Color', [0.7 0.7 0.7]);
xline(ax, 115, ':', 'Color', [0.7 0.7 0.7]);
set(ax, 'Color', 'w', 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'LineWidth', 1.0, 'XColor', 'k', 'YColor', 'k');
xlim(ax, [0, 170]);
end

function setup_smith_axes(ax)
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, [-1.1, 1.1, -1.1, 1.1]);
set(ax, 'Color', 'w', 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'LineWidth', 1.0, 'XColor', 'none', 'YColor', 'none');
box(ax, 'off');
draw_smith_grid(ax);
end

function save_validation_figure(fig, basePath)
set(fig, 'Color', 'w', 'InvertHardcopy', 'off');
axesList = findall(fig, 'Type', 'Axes');
set(axesList, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
drawnow;
savefig(fig, [basePath, '.fig']);
exportgraphics(fig, [basePath, '.jpg'], 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, [basePath, '.pdf'], 'ContentType', 'image', 'Resolution', 300, 'BackgroundColor', 'white');
close(fig);
end

function ensure_dir(pathStr)
if ~exist(pathStr, 'dir')
    mkdir(pathStr);
end
end

function write_text_file(pathStr, lines)
fid = fopen(pathStr, 'w');
if fid < 0
    error('PhaseII:WriteFailed', 'Unable to write file: %s', pathStr);
end
cleanup = onCleanup(@() fclose(fid));
for idxLine = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idxLine});
end
end

function out = ternary_text(tf, trueText, falseText)
if tf
    out = trueText;
else
    out = falseText;
end
end

function colorOut = tint_color(colorIn, mixFrac)
colorOut = colorIn + mixFrac * (1 - colorIn);
colorOut = min(max(colorOut, 0), 1);
end

function draw_smith_grid(ax)
gridColor = 0.72 * [1 1 1];
gridThin = 0.85 * [1 1 1];
theta = linspace(0, 2*pi, 1200);

plot(ax, cos(theta), sin(theta), 'Color', 0.55 * [1 1 1], 'LineWidth', 1.0);
plot(ax, [-1 1], [0 0], 'Color', 0.55 * [1 1 1], 'LineWidth', 1.0);

rValsMajor = [0.2 0.5 1 2 5];
rValsMinor = [0.1 0.3 0.7 1.5 3];
xValsMajor = [0.2 0.5 1 2 5];
xValsMinor = [0.1 0.3 0.7 1.5 3];

for r = rValsMinor
    draw_const_r(ax, r, gridThin, 0.7);
end
for r = rValsMajor
    draw_const_r(ax, r, gridColor, 0.9);
end
for x = xValsMinor
    draw_const_x(ax, x, gridThin, 0.7);
    draw_const_x(ax, -x, gridThin, 0.7);
end
for x = xValsMajor
    draw_const_x(ax, x, gridColor, 0.9);
    draw_const_x(ax, -x, gridColor, 0.9);
end
end

function draw_const_r(ax, r, colorVal, lineWidth)
cx = r / (1 + r);
rad = 1 / (1 + r);
t = linspace(0, 2*pi, 4000);
x = cx + rad * cos(t);
y = rad * sin(t);
mask = (x.^2 + y.^2) <= 1 + 1e-9;
plot(ax, x(mask), y(mask), 'Color', colorVal, 'LineWidth', lineWidth);
end

function draw_const_x(ax, xConst, colorVal, lineWidth)
cx = 1;
cy = 1 / xConst;
rad = abs(1 / xConst);
t = linspace(0, 2*pi, 4000);
x = cx + rad * cos(t);
y = cy + rad * sin(t);
mask = (x.^2 + y.^2) <= 1 + 1e-9 & x <= 1 + 1e-9;
plot(ax, x(mask), y(mask), 'Color', colorVal, 'LineWidth', lineWidth);
end
