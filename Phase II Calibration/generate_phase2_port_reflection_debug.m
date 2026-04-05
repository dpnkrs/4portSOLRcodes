function generate_phase2_port_reflection_debug(config, bandResults, referenceStandards)
%GENERATE_PHASE2_PORT_REFLECTION_DEBUG Audit which one-port standard drives bad local solves.

debug = struct();
noteLines = {
    'Phase II Port Reflection Debug Summary'
    '===================================='
    ''
    'Interpretation of improvement metric:'
    '  improvement = baseline corrected-thru error - corrected-thru error after'
    '  replacing one measured O/S/L reflection with the stitched Phase I reference.'
    '  Positive improvement means that measurement is likely contributing to the blowup.'
    ''
    };

for idxPort = 1:numel(config.port_labels)
    portLabel = config.port_labels{idxPort};
    entries = collect_port_entries(portLabel, bandResults, referenceStandards, config);
    if isempty(entries)
        continue;
    end

    figRef = create_reference_overlay_figure(config, portLabel, entries, referenceStandards);
    save_debug_figure(figRef, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_PortReflection_%s', portLabel)));

    figSens = create_sensitivity_figure(config, portLabel, entries);
    save_debug_figure(figSens, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_PortSensitivity_%s', portLabel)));

    [summaryLines, aggregateSummary] = summarize_port_entries(portLabel, entries, config.bands);
    noteLines = [noteLines; summaryLines; {''}]; %#ok<AGROW>

    debug.(portLabel) = struct( ...
        'entries', {entries}, ...
        'aggregate', aggregateSummary); %#ok<STRNU>
end

save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_PORT_REFLECTION.mat'), 'debug');
write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_PortReflection_Summary.txt'), noteLines);
end

function entries = collect_port_entries(portLabel, bandResults, referenceStandards, config)
entries = struct([]);

for idxBand = 1:numel(bandResults)
    bandResult = bandResults{idxBand};
    for idxPair = 1:numel(bandResult.reference_pair_results)
        pairResult = bandResult.reference_pair_results(idxPair);
        ports = pairResult.ports;
        localLabels = {sprintf('P%d', ports(1)), sprintf('P%d', ports(2))};
        sideIndex = find(strcmpi(localLabels, portLabel), 1, 'first');
        if isempty(sideIndex)
            continue;
        end

        freq = pairResult.freq(:);
        refs = struct( ...
            'Open', interpolate_reference_gamma(referenceStandards.Open, freq, bandResult.band), ...
            'Short', interpolate_reference_gamma(referenceStandards.Short, freq, bandResult.band), ...
            'Load', interpolate_reference_gamma(referenceStandards.Load, freq, bandResult.band));

        measPort1 = pairResult.standards_measured_port1;
        measPort2 = pairResult.standards_measured_port2;
        if sideIndex == 1
            meas = measPort1;
            baselineLocal = solve_local_reflection_terms(meas.Short, meas.Open, meas.Load, refs.Short, refs.Open, refs.Load);
            tBase = pairResult.error_terms.t11(:);
        else
            meas = measPort2;
            baselineLocal = solve_local_reflection_terms(meas.Short, meas.Open, meas.Load, refs.Short, refs.Open, refs.Load);
            tBase = pairResult.error_terms.t22(:);
        end

        phase1Target = load_phase1_thru_target(config, pairResult.geometry, pairResult.pair, bandResult.band);
        targetS21 = interpolate_network_element(phase1Target.freq(:), squeeze(phase1Target.S(2, 1, :)), freq);
        correctedS21 = squeeze(pairResult.reference_corrected(2, 1, :));
        baselineErr = abs(correctedS21 - targetS21);

        replacements = struct();
        stdNames = {'Open', 'Short', 'Load'};
        for idxStd = 1:numel(stdNames)
            stdName = stdNames{idxStd};
            measPort1Mod = measPort1;
            measPort2Mod = measPort2;
            if sideIndex == 1
                measPort1Mod.(stdName) = refs.(stdName);
            else
                measPort2Mod.(stdName) = refs.(stdName);
            end

            modifiedTerms = calculate_error_terms_solr_from_gamma(freq, ...
                measPort1Mod.Short, measPort1Mod.Open, measPort1Mod.Load, ...
                measPort2Mod.Short, measPort2Mod.Open, measPort2Mod.Load, ...
                pairResult.reference_switch_corrected, refs.Short, refs.Open, refs.Load, 0);

            correctedRef = apply_error_terms_to_network(pairResult.reference_switch_corrected, modifiedTerms, config.cal_den_floor);
            modifiedS21 = squeeze(correctedRef(2, 1, :));
            modifiedErr = abs(modifiedS21 - targetS21);

            if sideIndex == 1
                localSolve = solve_local_reflection_terms(measPort1Mod.Short, measPort1Mod.Open, measPort1Mod.Load, refs.Short, refs.Open, refs.Load);
                tModified = modifiedTerms.t11(:);
            else
                localSolve = solve_local_reflection_terms(measPort2Mod.Short, measPort2Mod.Open, measPort2Mod.Load, refs.Short, refs.Open, refs.Load);
                tModified = modifiedTerms.t22(:);
            end

            replacements.(stdName) = struct( ...
                'local_solve', localSolve, ...
                'error_terms', modifiedTerms, ...
                'corrected_reference', correctedRef, ...
                'corrected_s21', modifiedS21, ...
                'corrected_s21_error', modifiedErr, ...
                'error_improvement', baselineErr - modifiedErr, ...
                't_modified', tModified, ...
                't_change_db', local_mag_db(tModified) - local_mag_db(tBase));
        end

        entry = struct();
        entry.port = portLabel;
        entry.side_index = sideIndex;
        entry.band = bandResult.band;
        entry.freq = freq;
        entry.pair = pairResult.pair;
        entry.geometry = pairResult.geometry;
        entry.pair_label = sprintf('%s %s (%s side)', pairResult.pair, pairResult.geometry, portLabel);
        entry.measured = meas;
        entry.reference = refs;
        entry.phase1_target = phase1Target;
        entry.baseline = struct( ...
            'local_solve', baselineLocal, ...
            't_local', tBase, ...
            'target_s21', targetS21, ...
            'corrected_s21', correctedS21, ...
            'thru_error', baselineErr);
        entry.replacements = replacements;

        entries = [entries, entry]; %#ok<AGROW>
    end
end
end

function fig = create_reference_overlay_figure(config, portLabel, entries, referenceStandards)
fig = figure('Visible', config.figure_visible, 'Color', 'w');
tlo = tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
stdNames = {'Open', 'Short', 'Load'};
entryColors = lines(numel(entries));

legendHandles = gobjects(numel(entries) + 1, 1);
legendLabels = cell(numel(entries) + 1, 1);
legendReady = false;
legendAxes = gobjects(1, 1);

for idxStd = 1:numel(stdNames)
    stdName = stdNames{idxStd};
    [refFreq, refGamma] = flatten_reference_trace(referenceStandards.(stdName), config.bands);
    refFreqGHz = refFreq(:) / 1e9;

    axMag = nexttile(tlo);
    hold(axMag, 'on');
    grid(axMag, 'on');
    title(axMag, sprintf('%s magnitude', stdName));
    xlabel(axMag, 'Frequency (GHz)');
    ylabel(axMag, sprintf('%s |\\Gamma| (dB)', stdName));
    xline(axMag, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(axMag, 115, ':', 'Color', [0.7 0.7 0.7]);
    hRef = plot(axMag, refFreqGHz, local_mag_db(refGamma), 'k-', 'LineWidth', 2.0);

    axPh = nexttile(tlo);
    hold(axPh, 'on');
    grid(axPh, 'on');
    title(axPh, sprintf('%s phase', stdName));
    xlabel(axPh, 'Frequency (GHz)');
    ylabel(axPh, sprintf('%s phase (deg)', stdName));
    xline(axPh, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(axPh, 115, ':', 'Color', [0.7 0.7 0.7]);
    plot(axPh, refFreqGHz, unwrap(angle(refGamma)) * 180 / pi, 'k-', 'LineWidth', 2.0);

    if ~legendReady
        legendLabels{1} = sprintf('%s Phase I reference', stdName);
        legendAxes = axMag;
    end

    for idxEntry = 1:numel(entries)
        trace = entries(idxEntry).measured.(stdName);
        plot(axMag, entries(idxEntry).freq / 1e9, local_mag_db(trace), '-', ...
            'Color', entryColors(idxEntry, :), 'LineWidth', 1.1);

        alignedPhase = align_phase_to_reference(trace, entries(idxEntry).reference.(stdName));
        plot(axPh, entries(idxEntry).freq / 1e9, alignedPhase, '-', ...
            'Color', entryColors(idxEntry, :), 'LineWidth', 1.1);

        if ~legendReady && idxStd == 1
            legendLabels{idxEntry + 1} = entries(idxEntry).pair_label;
        end
    end
    legendReady = true;
end

legendHandles(1) = plot(legendAxes, nan, nan, 'k-', 'LineWidth', 2.0);
for idxEntry = 1:numel(entries)
    legendHandles(idxEntry + 1) = plot(legendAxes, nan, nan, '-', 'Color', entryColors(idxEntry, :), 'LineWidth', 1.2);
end
lgd = legend(legendAxes, legendHandles, legendLabels, 'Location', 'eastoutside');
set(lgd, 'Interpreter', 'none');
title(tlo, sprintf('Phase II reflection inputs vs Phase I references - %s', portLabel));
end

function fig = create_sensitivity_figure(config, portLabel, entries)
fig = figure('Visible', config.figure_visible, 'Color', 'w');
tlo = tiledlayout(fig, 4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
entryColors = lines(numel(entries));
stdNames = {'Open', 'Short', 'Load'};

legendHandles = gobjects(numel(entries), 1);
legendLabels = cell(numel(entries), 1);
for idxEntry = 1:numel(entries)
    legendLabels{idxEntry} = entries(idxEntry).pair_label;
end
legendAxes = gobjects(1, 1);

ax = nexttile(tlo);
hold(ax, 'on');
for idxEntry = 1:numel(entries)
    plot(ax, entries(idxEntry).freq / 1e9, entries(idxEntry).baseline.thru_error, '-', ...
        'Color', entryColors(idxEntry, :), 'LineWidth', 1.2);
end
grid(ax, 'on');
xline(ax, 67, ':', 'Color', [0.7 0.7 0.7]);
xline(ax, 115, ':', 'Color', [0.7 0.7 0.7]);
xlabel(ax, 'Frequency (GHz)');
ylabel(ax, '|S_{21,corr}-S_{21,target}|');
title(ax, sprintf('%s baseline corrected-thru error', portLabel));
legendAxes = ax;

for idxStd = 1:numel(stdNames)
    stdName = stdNames{idxStd};
    ax = nexttile(tlo);
    hold(ax, 'on');
    for idxEntry = 1:numel(entries)
        improvement = entries(idxEntry).replacements.(stdName).error_improvement;
        plot(ax, entries(idxEntry).freq / 1e9, improvement, '-', ...
            'Color', entryColors(idxEntry, :), 'LineWidth', 1.2);
    end
    yline(ax, 0, '--', 'Color', [0.35 0.35 0.35]);
    grid(ax, 'on');
    xline(ax, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax, 115, ':', 'Color', [0.7 0.7 0.7]);
    xlabel(ax, 'Frequency (GHz)');
    ylabel(ax, '\Delta error');
    title(ax, sprintf('%s improvement if %s is replaced by stitched reference', portLabel, stdName));
end

for idxEntry = 1:numel(entries)
    legendHandles(idxEntry) = plot(legendAxes, nan, nan, '-', 'Color', entryColors(idxEntry, :), 'LineWidth', 1.2);
end
lgd = legend(legendAxes, legendHandles, legendLabels, 'Location', 'eastoutside');
set(lgd, 'Interpreter', 'none');
title(tlo, sprintf('Phase II one-port measurement sensitivity through corrected-thru error - %s', portLabel));
end

function [linesOut, aggregateSummary] = summarize_port_entries(portLabel, entries, bandList)
linesOut = {
    sprintf('%s', portLabel)
    repmat('-', 1, numel(portLabel))
    };
aggregateSummary = struct();
stdNames = {'Open', 'Short', 'Load'};

for idxBand = 1:numel(bandList)
    bandName = bandList{idxBand};
    bandEntries = entries(strcmpi({entries.band}, bandName));
    if isempty(bandEntries)
        continue;
    end

    linesOut{end + 1} = sprintf('  %s', bandName); %#ok<AGROW>
    aggregateImprovements = nan(numel(bandEntries), numel(stdNames));
    baselineErrors = nan(numel(bandEntries), 1);

    for idxEntry = 1:numel(bandEntries)
        entry = bandEntries(idxEntry);
        baselineErrors(idxEntry) = median(entry.baseline.thru_error, 'omitnan');
        linesOut{end + 1} = sprintf('    %s', entry.pair_label); %#ok<AGROW>
        linesOut{end + 1} = sprintf('      Median baseline |t_local| (dB): %.3f', median(local_mag_db(entry.baseline.t_local), 'omitnan')); %#ok<AGROW>
        linesOut{end + 1} = sprintf('      Median baseline corrected-thru error: %.4f', baselineErrors(idxEntry)); %#ok<AGROW>
        linesOut{end + 1} = sprintf('      Median rcond(A): %.4e', median(entry.baseline.local_solve.rcond, 'omitnan')); %#ok<AGROW>

        bestStd = 'N/A';
        bestVal = -Inf;
        for idxStd = 1:numel(stdNames)
            stdName = stdNames{idxStd};
            improvement = median(entry.replacements.(stdName).error_improvement, 'omitnan');
            aggregateImprovements(idxEntry, idxStd) = improvement;
            tShift = median(entry.replacements.(stdName).t_change_db, 'omitnan');
            linesOut{end + 1} = sprintf('      Replace %-5s -> median thru-error improvement %.4f, median \\Delta|t| %.3f dB', ...
                stdName, improvement, tShift); %#ok<AGROW>
            if improvement > bestVal
                bestVal = improvement;
                bestStd = stdName;
            end
        end
        linesOut{end + 1} = sprintf('      Dominant suspect for this pair-side: %s', bestStd); %#ok<AGROW>
    end

    [~, worstEntryIdx] = max(baselineErrors);
    meanImprovements = mean(aggregateImprovements, 1, 'omitnan');
    [bestMeanVal, bestMeanIdx] = max(meanImprovements);
    if isempty(bestMeanIdx) || ~isfinite(bestMeanVal)
        bestMeanStd = 'N/A';
    else
        bestMeanStd = stdNames{bestMeanIdx};
    end

    dominantWorstStd = 'N/A';
    if ~isempty(worstEntryIdx) && all(isfinite(aggregateImprovements(worstEntryIdx, :)))
        [~, bestWorstIdx] = max(aggregateImprovements(worstEntryIdx, :));
        dominantWorstStd = stdNames{bestWorstIdx};
    elseif ~isempty(worstEntryIdx)
        row = aggregateImprovements(worstEntryIdx, :);
        [~, bestWorstIdx] = max(replace_nonfinite(row, -Inf));
        if isfinite(row(bestWorstIdx))
            dominantWorstStd = stdNames{bestWorstIdx};
        end
    end

    linesOut{end + 1} = sprintf('    Worst baseline pair-side in %s: %s', bandName, bandEntries(worstEntryIdx).pair_label); %#ok<AGROW>
    linesOut{end + 1} = sprintf('    Dominant suspect on worst pair-side: %s', dominantWorstStd); %#ok<AGROW>
    linesOut{end + 1} = sprintf('    Aggregate dominant suspect for %s in %s: %s', portLabel, bandName, bestMeanStd); %#ok<AGROW>

    ranking = rank_band_suspects(bandEntries, stdNames);
    if ~isempty(ranking)
        linesOut{end + 1} = '    Ranked suspect measurements:'; %#ok<AGROW>
        topN = min(3, numel(ranking));
        for idxRank = 1:topN
            item = ranking(idxRank);
            linesOut{end + 1} = sprintf('      %d. %s | replace %-5s | median improvement %.4f', ...
                idxRank, item.pair_label, item.standard, item.improvement); %#ok<AGROW>
        end
    end

    aggregateSummary.(matlab.lang.makeValidName(strrep(bandName, '-', '_'))) = struct( ...
        'worst_pair_label', bandEntries(worstEntryIdx).pair_label, ...
        'dominant_worst_pair_standard', dominantWorstStd, ...
        'mean_improvements', meanImprovements, ...
        'aggregate_dominant_standard', bestMeanStd, ...
        'ranking', ranking);
end
end

function ranking = rank_band_suspects(bandEntries, stdNames)
ranking = struct('pair_label', {}, 'standard', {}, 'improvement', {});
for idxEntry = 1:numel(bandEntries)
    entry = bandEntries(idxEntry);
    for idxStd = 1:numel(stdNames)
        stdName = stdNames{idxStd};
        improvement = median(entry.replacements.(stdName).error_improvement, 'omitnan');
        if ~isfinite(improvement)
            continue;
        end
        ranking(end + 1) = struct( ... %#ok<AGROW>
            'pair_label', entry.pair_label, ...
            'standard', stdName, ...
            'improvement', improvement);
    end
end
if isempty(ranking)
    return;
end
[~, order] = sort([ranking.improvement], 'descend');
ranking = ranking(order);
end

function [freq, gamma] = flatten_reference_trace(reference, bandNames)
if isfield(reference, 'bands')
    freq = [];
    gamma = [];
    for idxBand = 1:numel(bandNames)
        bandKey = matlab.lang.makeValidName(strrep(bandNames{idxBand}, '-', '_'));
        refBand = reference.bands.(bandKey);
        freq = [freq; refBand.freq(:)]; %#ok<AGROW>
        gamma = [gamma; refBand.gamma(:)]; %#ok<AGROW>
    end
else
    freq = reference.freq(:);
    gamma = reference.gamma(:);
end
end

function gamma = interpolate_reference_gamma(reference, freq, bandName)
if isfield(reference, 'bands')
    bandKey = matlab.lang.makeValidName(strrep(bandName, '-', '_'));
    source = reference.bands.(bandKey);
else
    source = reference;
end
gamma = interp1(source.freq(:), source.gamma(:), freq, 'pchip', 'extrap');
end

function phaseDeg = align_phase_to_reference(trace, ref)
phaseTrace = unwrap(angle(trace(:)));
phaseRef = unwrap(angle(ref(:)));
offsetCycles = round(median((phaseRef - phaseTrace) / (2 * pi), 'omitnan'));
phaseDeg = (phaseTrace + 2 * pi * offsetCycles) * 180 / pi;
end

function metrics = solve_local_reflection_terms(measShort, measOpen, measLoad, gammaShort, gammaOpen, gammaLoad)
nFreq = numel(measShort);
metrics = struct( ...
    'rcond', nan(nFreq, 1), ...
    'e00', nan(nFreq, 1), ...
    'e11', nan(nFreq, 1), ...
    't', nan(nFreq, 1));

for idx = 1:nFreq
    A = [ ...
        1, measShort(idx) * gammaShort(idx), gammaShort(idx); ...
        1, measOpen(idx) * gammaOpen(idx), gammaOpen(idx); ...
        1, measLoad(idx) * gammaLoad(idx), gammaLoad(idx)];
    b = [measShort(idx); measOpen(idx); measLoad(idx)];

    metrics.rcond(idx) = rcond(A);
    if metrics.rcond(idx) < 1e-12
        continue;
    end

    sol = A \ b;
    metrics.e00(idx) = sol(1);
    metrics.e11(idx) = sol(2);
    metrics.t(idx) = sol(3) + sol(1) * sol(2);
end
end

function Scorrected = apply_error_terms_to_network(Sraw, errorTerms, denominatorFloor)
Scorrected = nan(size(Sraw));
for idxFreq = 1:size(Sraw, 3)
    Scorrected(:, :, idxFreq) = apply_error_correction_8term_local(Sraw(:, :, idxFreq), errorTerms, idxFreq, denominatorFloor);
end
end

function target = load_phase1_thru_target(config, geomKey, pairKey, bandName)
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
        error('Unsupported geometry %s for Phase I target lookup.', geomKey);
end

loaded = load(fullfile(config.phase1_output_mat_dir, baseName));
target = loaded.reciprocalResult;
end

function values = local_mag_db(trace)
values = 20 * log10(max(abs(trace), 1e-12));
end

function traceOut = interpolate_network_element(freqIn, traceIn, freqOut)
realPart = interp1(freqIn, real(traceIn), freqOut, 'linear', 'extrap');
imagPart = interp1(freqIn, imag(traceIn), freqOut, 'linear', 'extrap');
traceOut = complex(realPart, imagPart);
end

function out = replace_nonfinite(values, replacement)
out = values;
out(~isfinite(out)) = replacement;
end

function save_debug_figure(fig, basePath)
set(fig, 'Color', 'w');
ax = findall(fig, 'type', 'axes');
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
drawnow;
exportgraphics(fig, [basePath, '.jpg'], 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, [basePath, '.pdf'], 'BackgroundColor', 'white', 'ContentType', 'vector');
savefig(fig, [basePath, '.fig']);
close(fig);
end

function write_text_file(filePath, lines)
fid = fopen(filePath, 'w');
cleanup = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end
