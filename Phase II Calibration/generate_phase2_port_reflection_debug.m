function generate_phase2_port_reflection_debug(config, bandResults, referenceStandards)
%GENERATE_PHASE2_PORT_REFLECTION_DEBUG Audit one-port measurements driving local SOLR solves.

debug = struct();
noteLines = {};
noteLines{end + 1} = 'Phase II Port Reflection Debug Summary';
noteLines{end + 1} = '====================================';
noteLines{end + 1} = '';

for idxPort = 1:numel(config.port_labels)
    portLabel = config.port_labels{idxPort};
    entries = collect_port_entries(portLabel, bandResults, referenceStandards, config);
    if isempty(entries)
        continue;
    end

    [figRef, refSummary] = create_reference_overlay_figure(config, portLabel, entries, referenceStandards);
    save_debug_figure(figRef, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_PortReflection_%s', portLabel)));

    [figSens, sensSummary] = create_sensitivity_figure(config, portLabel, entries);
    save_debug_figure(figSens, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_PortSensitivity_%s', portLabel)));

    [summaryLines, aggregateSummary] = summarize_port_entries(portLabel, entries, config.bands);

    debug.(portLabel) = struct( ...
        'entries', {entries}, ...
        'reference_summary', refSummary, ...
        'sensitivity_summary', sensSummary, ...
        'aggregate', aggregateSummary); %#ok<STRNU>

    noteLines = [noteLines, summaryLines, {''}]; %#ok<AGROW>
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
            'Open', interpolate_reference_gamma(referenceStandards.Open, freq), ...
            'Short', interpolate_reference_gamma(referenceStandards.Short, freq), ...
            'Load', interpolate_reference_gamma(referenceStandards.Load, freq));

        if sideIndex == 1
            meas = pairResult.standards_measured_port1;
            baselineLocalSolve = solve_local_reflection_terms(meas.Short, meas.Open, meas.Load, refs.Short, refs.Open, refs.Load);
            tBase = pairResult.error_terms.t11(:);
        else
            meas = pairResult.standards_measured_port2;
            baselineLocalSolve = solve_local_reflection_terms(meas.Short, meas.Open, meas.Load, refs.Short, refs.Open, refs.Load);
            tBase = pairResult.error_terms.t22(:);
        end

        phase1Target = load_phase1_thru_target(config, pairResult.geometry, pairResult.pair, bandResult.band);
        targetS21 = squeeze(phase1Target.S(2, 1, :));
        corrS21 = squeeze(pairResult.reference_corrected(2, 1, :));
        baselineErr = abs(corrS21 - targetS21);

        replacements = struct();
        stdNames = {'Open', 'Short', 'Load'};
        for idxStd = 1:numel(stdNames)
            stdName = stdNames{idxStd};
            [modifiedTerms, modifiedMetrics] = replace_standard_and_resolve(pairResult, refs, sideIndex, stdName);
            correctedRef = apply_error_terms_to_network(pairResult.reference_switch_corrected, modifiedTerms, config.cal_den_floor);
            modifiedErr = abs(squeeze(correctedRef(2, 1, :)) - targetS21);
            if sideIndex == 1
                tModified = modifiedTerms.t11(:);
            else
                tModified = modifiedTerms.t22(:);
            end

            replacements.(stdName) = struct( ...
                'error_terms', modifiedTerms, ...
                'metrics', modifiedMetrics, ...
                'corrected_reference', correctedRef, ...
                'corrected_s21_error', modifiedErr, ...
                'error_improvement', baselineErr - modifiedErr, ...
                't_modified', tModified, ...
                't_change_db', local_mag_db(tModified) - local_mag_db(tBase));
        end

        newEntry = struct();
        newEntry.pair = pairResult.pair;
        newEntry.geometry = pairResult.geometry;
        newEntry.band = bandResult.band;
        newEntry.freq = freq;
        newEntry.port = portLabel;
        newEntry.side_index = sideIndex;
        newEntry.pair_label = sprintf('%s %s (%s side)', pairResult.pair, pairResult.geometry, portLabel);
        newEntry.measured = meas;
        newEntry.reference = refs;
        newEntry.baseline = struct( ...
            't_local', tBase, ...
            'thru_error', baselineErr, ...
            'target_s21', targetS21, ...
            'corrected_s21', corrS21, ...
            'local_solve', baselineLocalSolve);
        newEntry.replacements = replacements;
        newEntry.phase1_target = phase1Target;

        entries = [entries, newEntry]; %#ok<AGROW>
    end
end
end

function [fig, summary] = create_reference_overlay_figure(config, portLabel, entries, referenceStandards)
fig = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
stdNames = {'Open', 'Short', 'Load'};
entryColors = lines(numel(entries));
summary = struct();

for idxStd = 1:numel(stdNames)
    stdName = stdNames{idxStd};

    axMag = nexttile;
    hold(axMag, 'on');
    title(axMag, sprintf('%s magnitude', stdName));
    xlabel(axMag, 'Frequency (GHz)');
    ylabel(axMag, sprintf('%s |\\Gamma| (dB)', stdName));
    grid(axMag, 'on');
    xline(axMag, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(axMag, 115, ':', 'Color', [0.7 0.7 0.7]);

    refFreqGHz = referenceStandards.(stdName).freq / 1e9;
    plot(axMag, refFreqGHz, local_mag_db(referenceStandards.(stdName).gamma), ...
        'k-', 'LineWidth', 2.0, 'DisplayName', sprintf('%s stitched reference', stdName));

    magOffsets = nan(1, numel(entries));
    for idxEntry = 1:numel(entries)
        trace = entries(idxEntry).measured.(stdName);
        ref = entries(idxEntry).reference.(stdName);
        plot(axMag, entries(idxEntry).freq / 1e9, local_mag_db(trace), ...
            '-', 'Color', entryColors(idxEntry, :), 'LineWidth', 1.1, 'DisplayName', entries(idxEntry).pair_label);
        magOffsets(idxEntry) = median(local_mag_db(trace) - local_mag_db(ref), 'omitnan');
    end
    summary.(stdName).median_mag_offset_db = median(magOffsets, 'omitnan');

    axPh = nexttile;
    hold(axPh, 'on');
    title(axPh, sprintf('%s phase', stdName));
    xlabel(axPh, 'Frequency (GHz)');
    ylabel(axPh, sprintf('%s phase (deg)', stdName));
    grid(axPh, 'on');
    xline(axPh, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(axPh, 115, ':', 'Color', [0.7 0.7 0.7]);

    refPhase = unwrap(angle(referenceStandards.(stdName).gamma)) * 180 / pi;
    plot(axPh, refFreqGHz, refPhase, 'k-', 'LineWidth', 2.0, 'DisplayName', sprintf('%s stitched reference', stdName));

    phaseOffsets = nan(1, numel(entries));
    for idxEntry = 1:numel(entries)
        alignedPhase = align_phase_to_reference(entries(idxEntry).measured.(stdName), entries(idxEntry).reference.(stdName));
        plot(axPh, entries(idxEntry).freq / 1e9, alignedPhase, ...
            '-', 'Color', entryColors(idxEntry, :), 'LineWidth', 1.1, 'DisplayName', entries(idxEntry).pair_label);
        refPhaseLocal = unwrap(angle(entries(idxEntry).reference.(stdName))) * 180 / pi;
        phaseOffsets(idxEntry) = median(alignedPhase - refPhaseLocal, 'omitnan');
    end
    summary.(stdName).median_phase_offset_deg = median(phaseOffsets, 'omitnan');
end

lgd = legend(findall(fig, 'Type', 'Line', '-not', 'LineStyle', ':'), 'Location', 'bestoutside');
set(lgd, 'Interpreter', 'none');
sgtitle(fig, sprintf('Phase II reflection inputs vs Phase I stitched references - %s', portLabel));
end

function [fig, summary] = create_sensitivity_figure(config, portLabel, entries)
fig = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig, 5, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
entryColors = lines(numel(entries));
summary = struct();

ax1 = nexttile;
hold(ax1, 'on');
for idxEntry = 1:numel(entries)
    plot(ax1, entries(idxEntry).freq / 1e9, local_mag_db(entries(idxEntry).baseline.t_local), ...
        '-', 'Color', entryColors(idxEntry, :), 'LineWidth', 1.3, 'DisplayName', entries(idxEntry).pair_label);
end
grid(ax1, 'on');
xline(ax1, 67, ':', 'Color', [0.7 0.7 0.7]);
xline(ax1, 115, ':', 'Color', [0.7 0.7 0.7]);
xlabel(ax1, 'Frequency (GHz)');
ylabel(ax1, '|t_{local}| (dB)');
title(ax1, sprintf('%s local tracking term used in one-port solve', portLabel));
legend(ax1, 'Location', 'bestoutside');

ax2 = nexttile;
hold(ax2, 'on');
baselineMedians = nan(1, numel(entries));
for idxEntry = 1:numel(entries)
    plot(ax2, entries(idxEntry).freq / 1e9, entries(idxEntry).baseline.thru_error, ...
        '-', 'Color', entryColors(idxEntry, :), 'LineWidth', 1.3);
    baselineMedians(idxEntry) = median(entries(idxEntry).baseline.thru_error, 'omitnan');
end
grid(ax2, 'on');
xline(ax2, 67, ':', 'Color', [0.7 0.7 0.7]);
xline(ax2, 115, ':', 'Color', [0.7 0.7 0.7]);
xlabel(ax2, 'Frequency (GHz)');
ylabel(ax2, '|S_{21,corr}-S_{21,target}|');
title(ax2, sprintf('%s corrected-thru error before replacement', portLabel));
summary.baseline_thru_error_median = median(baselineMedians, 'omitnan');

stdNames = {'Open', 'Short', 'Load'};
for idxStd = 1:numel(stdNames)
    stdName = stdNames{idxStd};
    ax = nexttile;
    hold(ax, 'on');
    improvementMedians = nan(1, numel(entries));
    for idxEntry = 1:numel(entries)
        improvement = entries(idxEntry).replacements.(stdName).error_improvement;
        plot(ax, entries(idxEntry).freq / 1e9, improvement, ...
            '-', 'Color', entryColors(idxEntry, :), 'LineWidth', 1.3);
        improvementMedians(idxEntry) = median(improvement, 'omitnan');
    end
    yline(ax, 0, '--', 'Color', [0.35 0.35 0.35]);
    grid(ax, 'on');
    xline(ax, 67, ':', 'Color', [0.7 0.7 0.7]);
    xline(ax, 115, ':', 'Color', [0.7 0.7 0.7]);
    xlabel(ax, 'Frequency (GHz)');
    ylabel(ax, '\Delta thru error');
    title(ax, sprintf('%s improvement if %s is replaced by Phase I reference', portLabel, stdName));
    summary.(stdName).median_improvement = median(improvementMedians, 'omitnan');
end

sgtitle(fig, sprintf('Phase II replacement sensitivity through corrected reference-thru error - %s', portLabel));
end

function [linesOut, aggregateSummary] = summarize_port_entries(portLabel, entries, bandList)
linesOut = {};
linesOut{end + 1} = sprintf('%s', portLabel);
linesOut{end + 1} = repmat('-', 1, numel(portLabel));

aggregateSummary = struct();
stdNames = {'Open', 'Short', 'Load'};

for idxBand = 1:numel(bandList)
    bandName = bandList{idxBand};
    bandEntries = entries(strcmpi({entries.band}, bandName));
    if isempty(bandEntries)
        continue;
    end

    linesOut{end + 1} = sprintf('  %s', bandName);
    aggregateImprovement = nan(1, numel(stdNames));

    for idxEntry = 1:numel(bandEntries)
        entry = bandEntries(idxEntry);
        linesOut{end + 1} = sprintf('    %s', entry.pair_label);
        linesOut{end + 1} = sprintf('      Median baseline |t_local| (dB): %.3f', median(local_mag_db(entry.baseline.t_local), 'omitnan'));
        linesOut{end + 1} = sprintf('      Median baseline corrected-thru error: %.4f', median(entry.baseline.thru_error, 'omitnan'));
        linesOut{end + 1} = sprintf('      Median rcond(A): %.4e', median(entry.baseline.local_solve.rcond, 'omitnan'));

        dominantStd = '';
        dominantVal = -Inf;
        for idxStd = 1:numel(stdNames)
            stdName = stdNames{idxStd};
            improvement = median(entry.replacements.(stdName).error_improvement, 'omitnan');
            tShift = median(entry.replacements.(stdName).t_change_db, 'omitnan');
            linesOut{end + 1} = sprintf('      Replace %-5s -> median thru-error improvement %.4f, median \\Delta|t| %.3f dB', ...
                stdName, improvement, tShift);
            aggregateImprovement(idxStd) = nansum_local([aggregateImprovement(idxStd), improvement]);
            if improvement > dominantVal
                dominantVal = improvement;
                dominantStd = stdName;
            end
        end
        linesOut{end + 1} = sprintf('      Dominant suspect for this pair-side: %s', dominantStd);
    end

    [bestVal, bestIdx] = max(aggregateImprovement);
    if isempty(bestIdx) || ~isfinite(bestVal)
        bestStd = 'N/A';
    else
        bestStd = stdNames{bestIdx};
    end

    linesOut{end + 1} = sprintf('    Aggregate dominant suspect for %s in %s: %s', portLabel, bandName, bestStd);
    aggregateSummary.(matlab.lang.makeValidName(strrep(bandName, '-', '_'))) = struct( ...
        'dominant_standard', bestStd, ...
        'aggregate_improvement', aggregateImprovement);
end
end

function refs = interpolate_reference_gamma(reference, freq)
refs = interp1(reference.freq, reference.gamma, freq, 'pchip', 'extrap');
end

function phaseDeg = align_phase_to_reference(trace, ref)
phaseTrace = unwrap(angle(trace));
phaseRef = unwrap(angle(ref));
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

function [modifiedTerms, sideMetrics] = replace_standard_and_resolve(pairResult, refs, sideIndex, stdName)
measPort1 = pairResult.standards_measured_port1;
measPort2 = pairResult.standards_measured_port2;

if sideIndex == 1
    measPort1.(stdName) = refs.(stdName);
else
    measPort2.(stdName) = refs.(stdName);
end

sideMetrics = solve_local_reflection_terms( ...
    measPort1.Short, measPort1.Open, measPort1.Load, refs.Short, refs.Open, refs.Load);
if sideIndex == 2
    sideMetrics = solve_local_reflection_terms( ...
        measPort2.Short, measPort2.Open, measPort2.Load, refs.Short, refs.Open, refs.Load);
end

modifiedTerms = calculate_error_terms_solr_from_gamma(pairResult.freq, ...
    measPort1.Short, measPort1.Open, measPort1.Load, ...
    measPort2.Short, measPort2.Open, measPort2.Load, ...
    pairResult.reference_switch_corrected, refs.Short, refs.Open, refs.Load, 0);
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

function total = nansum_local(values)
finiteMask = isfinite(values);
if ~any(finiteMask)
    total = NaN;
else
    total = sum(values(finiteMask));
end
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
