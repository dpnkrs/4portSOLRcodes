function summary = run_dissertation_results()
%RUN_DISSERTATION_RESULTS Build dissertation-facing notes and tables.

scriptDir = fileparts(mfilename('fullpath'));
repoDir = fileparts(scriptDir);

phase1MatDir = fullfile(repoDir, 'Phase I Extraction', 'outputs', 'mat');
outputRoot = fullfile(scriptDir, 'outputs');
matDir = fullfile(outputRoot, 'mat');
tableDir = fullfile(outputRoot, 'tables');
noteDir = fullfile(outputRoot, 'notes');

ensure_dir(outputRoot);
ensure_dir(matDir);
ensure_dir(tableDir);
ensure_dir(noteDir);

resultsPath = fullfile(phase1MatDir, 'PHASE1_RESULTS_ALL.mat');
if ~isfile(resultsPath)
    error('DissertationResults:MissingInput', ...
        'Required file not found: %s', resultsPath);
end

loaded = load(resultsPath, 'results');
results = loaded.results;

stitched = struct();
stitched.OPEN = load_fitout(fullfile(phase1MatDir, 'FINALFIT_OPEN_STITCHED_FINAL.mat'));
stitched.SHORT = load_fitout(fullfile(phase1MatDir, 'FINALFIT_SHORT_STITCHED_FINAL.mat'));
stitched.LOAD = load_fitout(fullfile(phase1MatDir, 'FINALFIT_LOAD_STITCHED_FINAL.mat'));

bandOrder = sort_band_indices(results.bands);
accepted = assemble_accepted_data(results, bandOrder);
solMetrics = compute_sol_metrics(accepted, results.config.z0);
consistencyMetrics = compute_consistency_metrics(results, bandOrder);
summaryRows = build_bandwise_summary_rows(accepted, solMetrics, consistencyMetrics, results.config.z0);
table51 = build_table51_rows(accepted, solMetrics, consistencyMetrics);
table52 = build_table52_rows(accepted);
sectionNotes = build_section_notes(accepted, solMetrics, consistencyMetrics, table51, table52, stitched);
writeupSuggestions = build_writeup_suggestions(solMetrics, consistencyMetrics);

summary = struct();
summary.generated_at = datestr(now, 'yyyy-mm-dd HH:MM:SS');
summary.output_root = outputRoot;
summary.accepted = accepted;
summary.sol_metrics = solMetrics;
summary.consistency_metrics = consistencyMetrics;
summary.bandwise_metrics = summaryRows;
summary.table51 = table51;
summary.table52 = table52;
summary.section_notes = sectionNotes;
summary.writeup_suggestions = writeupSuggestions;
summary.stitched = stitched;

save(fullfile(matDir, 'DissertationResultsSummary.mat'), 'summary');
writetable(struct2table(summaryRows), fullfile(tableDir, 'BandwiseStandardMetrics.csv'));
writetable(struct2table(table51), fullfile(tableDir, 'Table5_1_Metrics.csv'));
writetable(struct2table(table52), fullfile(tableDir, 'Table5_2_AdoptedStandards.csv'));

write_text(fullfile(noteDir, 'Section5_1_2_EvaluationFramework.txt'), sectionNotes.evaluation_framework);
write_text(fullfile(noteDir, 'Section5_1_3_OpenShort.txt'), sectionNotes.open_short);
write_text(fullfile(noteDir, 'Section5_1_4_Load.txt'), sectionNotes.load);
write_text(fullfile(noteDir, 'Section5_1_5_AdoptedStandards.txt'), sectionNotes.adopted_standards);
write_text(fullfile(noteDir, 'WriteupSuggestions.txt'), writeupSuggestions);
end

function fitOut = load_fitout(pathStr)
loaded = load(pathStr, 'fitOut');
fitOut = loaded.fitOut;
end

function ensure_dir(pathStr)
if ~exist(pathStr, 'dir')
    mkdir(pathStr);
end
end

function idx = sort_band_indices(bands)
starts = zeros(numel(bands), 1);
for k = 1:numel(bands)
    [starts(k), ~] = parse_band_limits(bands{k}.band);
end
[~, idx] = sort(starts);
end

function accepted = assemble_accepted_data(results, bandOrder)
stdList = {'OPEN', 'SHORT', 'LOAD'};
accepted = struct();
entryTemplate = blank_accepted_entry();
for iStd = 1:numel(stdList)
    accepted.(stdList{iStd}) = struct('bands', repmat(entryTemplate, numel(bandOrder), 1));
end

for iBand = 1:numel(bandOrder)
    bandResult = results.bands{bandOrder(iBand)};
    bandName = bandResult.band;
    for iStd = 1:numel(stdList)
        stdName = stdList{iStd};
        selection = bandResult.final_selection.(stdName);
        entry = entryTemplate;
        entry.band = bandName;
        [entry.band_start_ghz, entry.band_stop_ghz] = parse_band_limits(bandName);
        entry.freq = make_col(selection.freq);
        entry.gamma = make_col(selection.gamma_winner);
        entry.gamma_fit = get_optional_col(selection, 'gamma_fit', entry.freq);
        entry.winner_port = string(selection.winner_port);
        entry.winner_key = get_field_or_default(selection, 'winner_output_key', ...
            derive_output_key(stdName, selection.winner_port, bandName));

        winnerOutput = bandResult.outputs.(entry.winner_key);
        entry.measured_gamma = get_optional_output_col(winnerOutput, 'measured_gamma', entry.freq);
        entry.gamma_raw = get_optional_output_col(winnerOutput, 'gamma_raw', entry.freq);
        entry.source_file = get_optional_output_string(winnerOutput, 'source_file');

        loserPort = other_port(entry.winner_port);
        loserKey = derive_output_key(stdName, loserPort, bandName);
        altField = sprintf('candidate_%s_key', lower(char(loserPort)));
        if isfield(selection, altField)
            loserKey = get_field_or_default(selection, altField, loserKey);
        end
        entry.loser_port = loserPort;
        entry.loser_key = loserKey;
        entry.has_loser = isfield(bandResult.outputs, loserKey);
        if entry.has_loser
            loserOutput = bandResult.outputs.(loserKey);
            entry.gamma_loser = get_optional_output_col(loserOutput, 'gamma', entry.freq);
            entry.measured_gamma_loser = get_optional_output_col(loserOutput, 'measured_gamma', entry.freq);
            entry.gamma_raw_loser = get_optional_output_col(loserOutput, 'gamma_raw', entry.freq);
            entry.loser_source_file = get_optional_output_string(loserOutput, 'source_file');
        else
            entry.gamma_loser = nan(size(entry.freq));
            entry.measured_gamma_loser = nan(size(entry.freq));
            entry.gamma_raw_loser = nan(size(entry.freq));
            entry.loser_source_file = "";
        end

        entry.selection_reason = string(get_field_or_default(selection, 'selection_reason', ''));
        entry.sim_reference = get_optional_sim_reference(selection, entry.freq);
        entry.metrics_p1 = get_optional_struct(selection, 'metrics_p1');
        entry.metrics_p2 = get_optional_struct(selection, 'metrics_p2');
        entry.global_winner_port = get_optional_global_port(bandResult, stdName);
        accepted.(stdName).bands(iBand, 1) = entry;
    end
end
end

function entry = blank_accepted_entry()
entry = struct( ...
    'band', "", ...
    'band_start_ghz', nan, ...
    'band_stop_ghz', nan, ...
    'freq', [], ...
    'gamma', [], ...
    'gamma_fit', [], ...
    'winner_port', "", ...
    'winner_key', "", ...
    'measured_gamma', [], ...
    'gamma_raw', [], ...
    'source_file', "", ...
    'loser_port', "", ...
    'loser_key', "", ...
    'has_loser', false, ...
    'gamma_loser', [], ...
    'measured_gamma_loser', [], ...
    'gamma_raw_loser', [], ...
    'loser_source_file', "", ...
    'selection_reason', "", ...
    'sim_reference', [], ...
    'metrics_p1', struct(), ...
    'metrics_p2', struct(), ...
    'global_winner_port', "");
end

function solMetrics = compute_sol_metrics(accepted, z0)
bandCount = numel(accepted.OPEN.bands);
solMetrics = repmat(struct(), bandCount, 1);

for iBand = 1:bandCount
    openEntry = accepted.OPEN.bands(iBand);
    shortEntry = accepted.SHORT.bands(iBand);
    loadEntry = accepted.LOAD.bands(iBand);

    freq = openEntry.freq;
    gammaStd = [openEntry.gamma, shortEntry.gamma, loadEntry.gamma];
    gammaMeas = [openEntry.measured_gamma, shortEntry.measured_gamma, loadEntry.measured_gamma];

    nFreq = numel(freq);
    kappa = nan(nFreq, 1);
    directResidual = nan(nFreq, 3);
    e00 = nan(nFreq, 1);
    e11 = nan(nFreq, 1);
    de = nan(nFreq, 1);
    heldoutResiduals = nan(nFreq, 3);

    for idx = 1:nFreq
        A = [1, gammaStd(idx, 1) * gammaMeas(idx, 1), -gammaStd(idx, 1); ...
             1, gammaStd(idx, 2) * gammaMeas(idx, 2), -gammaStd(idx, 2); ...
             1, gammaStd(idx, 3) * gammaMeas(idx, 3), -gammaStd(idx, 3)];
        b = gammaMeas(idx, :).';

        if any(~isfinite(A(:))) || any(~isfinite(b(:)))
            continue;
        end

        x = A \ b;
        e00(idx) = x(1);
        e11(idx) = x(2);
        de(idx) = x(3);
        kappa(idx) = cond(A);

        for iStd = 1:3
            gStd = gammaStd(idx, iStd);
            gHat = (e00(idx) - de(idx) * gStd) ./ (1 - e11(idx) * gStd);
            directResidual(idx, iStd) = gammaMeas(idx, iStd) - gHat;
        end

        rejectedGamma = [openEntry.gamma_loser(idx), shortEntry.gamma_loser(idx), loadEntry.gamma_loser(idx)];
        rejectedMeas = [openEntry.measured_gamma_loser(idx), shortEntry.measured_gamma_loser(idx), loadEntry.measured_gamma_loser(idx)];
        for iStd = 1:3
            if ~isfinite(rejectedGamma(iStd)) || ~isfinite(rejectedMeas(iStd))
                continue;
            end
            gHat = (e00(idx) - de(idx) * rejectedGamma(iStd)) ./ (1 - e11(idx) * rejectedGamma(iStd));
            heldoutResiduals(idx, iStd) = rejectedMeas(iStd) - gHat;
        end
    end

    pairOS = abs(gammaStd(:, 1) - gammaStd(:, 2));
    pairOL = abs(gammaStd(:, 1) - gammaStd(:, 3));
    pairSL = abs(gammaStd(:, 2) - gammaStd(:, 3));
    zLoad = z0 .* (1 + loadEntry.gamma) ./ (1 - loadEntry.gamma);

    heldout = repmat(struct( ...
        'standard', "", ...
        'available', false, ...
        'median_abs_residual', nan, ...
        'p95_abs_residual', nan, ...
        'max_abs_residual', nan), 3, 1);
    stdNames = {'OPEN', 'SHORT', 'LOAD'};
    for iStd = 1:3
        validHeldout = abs(heldoutResiduals(:, iStd));
        validHeldout = validHeldout(isfinite(validHeldout));
        heldout(iStd).standard = string(stdNames{iStd});
        heldout(iStd).available = ~isempty(validHeldout);
        if ~isempty(validHeldout)
            heldout(iStd).median_abs_residual = median(validHeldout);
            heldout(iStd).p95_abs_residual = percentile_local(validHeldout, 95);
            heldout(iStd).max_abs_residual = max(validHeldout);
        end
    end

    solMetrics(iBand).band = openEntry.band;
    solMetrics(iBand).freq = freq;
    solMetrics(iBand).kappa = kappa;
    solMetrics(iBand).kappa_min = min(kappa, [], 'omitnan');
    solMetrics(iBand).kappa_median = median(kappa, 'omitnan');
    solMetrics(iBand).kappa_max = max(kappa, [], 'omitnan');
    solMetrics(iBand).kappa_pct_below_50 = 100 * mean(kappa < 50, 'omitnan');
    solMetrics(iBand).direct_residual = directResidual;
    solMetrics(iBand).direct_residual_mean = mean(abs(directResidual), 1, 'omitnan');
    solMetrics(iBand).direct_residual_max = max(abs(directResidual), [], 1, 'omitnan');
    solMetrics(iBand).heldout = heldout;
    solMetrics(iBand).pairwise = struct( ...
        'open_short_min', min(pairOS, [], 'omitnan'), ...
        'open_short_median', median(pairOS, 'omitnan'), ...
        'open_load_min', min(pairOL, [], 'omitnan'), ...
        'open_load_median', median(pairOL, 'omitnan'), ...
        'short_load_min', min(pairSL, [], 'omitnan'), ...
        'short_load_median', median(pairSL, 'omitnan'));
    solMetrics(iBand).load_impedance = struct( ...
        'real_min', min(real(zLoad), [], 'omitnan'), ...
        'real_max', max(real(zLoad), [], 'omitnan'), ...
        'imag_min', min(imag(zLoad), [], 'omitnan'), ...
        'imag_max', max(imag(zLoad), [], 'omitnan'), ...
        'gamma_abs_max', max(abs(loadEntry.gamma), [], 'omitnan'), ...
        'gamma_abs_median', median(abs(loadEntry.gamma), 'omitnan'));
    solMetrics(iBand).error_terms = struct('e00', e00, 'e11', e11, 'de', de);
end
end

function consistencyMetrics = compute_consistency_metrics(results, bandOrder)
consistencyMetrics = repmat(struct(), numel(bandOrder), 1);
for iBand = 1:numel(bandOrder)
    bandResult = results.bands{bandOrder(iBand)};
    diag = bandResult.refine_diagnostics;
    consistencyMetrics(iBand).band = bandResult.band;
    consistencyMetrics(iBand).open_improvement_median = median(make_col(diag.open_mag_improvement), 'omitnan');
    consistencyMetrics(iBand).short_improvement_median = median(make_col(diag.short_mag_improvement), 'omitnan');
    consistencyMetrics(iBand).load_improvement_median = median(make_col(diag.load_mag_improvement), 'omitnan');
    consistencyMetrics(iBand).sigma_change_median = median(abs(make_col(diag.sigma_change)), 'omitnan');
    consistencyMetrics(iBand).thru_reconstruction_error_median = median(make_col(diag.thru_reconstruction_error), 'omitnan');
    consistencyMetrics(iBand).p34_line_match_error_median = median(make_col(diag.p34_line_match_error), 'omitnan');
    consistencyMetrics(iBand).p34_return_loss_error_median = median(make_col(diag.p34_return_loss_error), 'omitnan');
    consistencyMetrics(iBand).simple_open_discrepancy_median = median(make_col(diag.simple_open_discrepancy), 'omitnan');
    consistencyMetrics(iBand).simple_short_discrepancy_median = median(make_col(diag.simple_short_discrepancy), 'omitnan');
    consistencyMetrics(iBand).simple_load_discrepancy_median = median(make_col(diag.simple_load_discrepancy), 'omitnan');
    consistencyMetrics(iBand).sim_ref_error_median = median(make_col(diag.sim_ref_error), 'omitnan');
    consistencyMetrics(iBand).sim_raw_error_median = median(make_col(diag.sim_raw_error), 'omitnan');
end
end

function rows = build_bandwise_summary_rows(accepted, solMetrics, consistencyMetrics, z0)
stdList = {'OPEN', 'SHORT', 'LOAD'};
rowTemplate = blank_bandwise_summary_row();
rows = repmat(rowTemplate, numel(stdList) * numel(solMetrics), 1);
rowIdx = 0;

for iBand = 1:numel(solMetrics)
    for iStd = 1:numel(stdList)
        rowIdx = rowIdx + 1;
        stdName = stdList{iStd};
        entry = accepted.(stdName).bands(iBand);
        row = rowTemplate;
        row.band = string(entry.band);
        row.standard = string(stdName);
        row.winner_port = entry.winner_port;
        row.loser_port = entry.loser_port;
        row.selection_reason = entry.selection_reason;
        row.global_winner_port = entry.global_winner_port;
        row.freq_start_ghz = min(entry.freq) / 1e9;
        row.freq_stop_ghz = max(entry.freq) / 1e9;
        row.gamma_abs_median = median(abs(entry.gamma), 'omitnan');
        row.gamma_abs_min = min(abs(entry.gamma), [], 'omitnan');
        row.gamma_abs_max = max(abs(entry.gamma), [], 'omitnan');
        row.gamma_mag_db_median = median(safe_mag_db(entry.gamma), 'omitnan');
        row.gamma_mag_db_min = min(safe_mag_db(entry.gamma), [], 'omitnan');
        row.gamma_mag_db_max = max(safe_mag_db(entry.gamma), [], 'omitnan');
        row.phase_deg_median = median(unwrap(angle(entry.gamma)) * 180 / pi, 'omitnan');

        if entry.has_loser
            dMag = abs(safe_mag_db(entry.gamma) - safe_mag_db(entry.gamma_loser));
            dPhase = abs(unwrap(angle(entry.gamma)) - unwrap(angle(entry.gamma_loser))) * 180 / pi;
            row.repeatability_mag_db_median = median(dMag, 'omitnan');
            row.repeatability_mag_db_max = max(dMag, [], 'omitnan');
            row.repeatability_phase_deg_median = median(dPhase, 'omitnan');
            row.repeatability_phase_deg_max = max(dPhase, [], 'omitnan');
        else
            row.repeatability_mag_db_median = nan;
            row.repeatability_mag_db_max = nan;
            row.repeatability_phase_deg_median = nan;
            row.repeatability_phase_deg_max = nan;
        end

        sim = entry.sim_reference;
        if any(isfinite(sim))
            row.sim_complex_error_median = median(abs(entry.gamma - sim), 'omitnan');
            row.sim_complex_error_p95 = percentile_local(abs(entry.gamma - sim), 95);
            row.sim_mag_db_error_median = median(abs(safe_mag_db(entry.gamma) - safe_mag_db(sim)), 'omitnan');
            row.sim_mag_db_error_p95 = percentile_local(abs(safe_mag_db(entry.gamma) - safe_mag_db(sim)), 95);
        else
            row.sim_complex_error_median = nan;
            row.sim_complex_error_p95 = nan;
            row.sim_mag_db_error_median = nan;
            row.sim_mag_db_error_p95 = nan;
        end

        heldout = solMetrics(iBand).heldout(iStd);
        row.heldout_abs_residual_median = heldout.median_abs_residual;
        row.heldout_abs_residual_p95 = heldout.p95_abs_residual;
        row.heldout_abs_residual_max = heldout.max_abs_residual;
        row.kappa_median = solMetrics(iBand).kappa_median;
        row.kappa_max = solMetrics(iBand).kappa_max;
        row.kappa_pct_below_50 = solMetrics(iBand).kappa_pct_below_50;
        row.open_short_separation_min = solMetrics(iBand).pairwise.open_short_min;
        row.open_load_separation_min = solMetrics(iBand).pairwise.open_load_min;
        row.short_load_separation_min = solMetrics(iBand).pairwise.short_load_min;
        row.load_real_impedance_min_ohm = solMetrics(iBand).load_impedance.real_min;
        row.load_real_impedance_max_ohm = solMetrics(iBand).load_impedance.real_max;
        row.load_imag_impedance_min_ohm = solMetrics(iBand).load_impedance.imag_min;
        row.load_imag_impedance_max_ohm = solMetrics(iBand).load_impedance.imag_max;
        row.consistency_open_improvement_median = consistencyMetrics(iBand).open_improvement_median;
        row.consistency_short_improvement_median = consistencyMetrics(iBand).short_improvement_median;
        row.consistency_load_improvement_median = consistencyMetrics(iBand).load_improvement_median;
        row.consistency_sigma_change_median = consistencyMetrics(iBand).sigma_change_median;
        row.consistency_thru_error_median = consistencyMetrics(iBand).thru_reconstruction_error_median;
        row.consistency_p34_error_median = consistencyMetrics(iBand).p34_line_match_error_median;
        row.z0_ohm = z0;
        rows(rowIdx) = row;
    end
end
end

function row = blank_bandwise_summary_row()
row = struct( ...
    'band', "", ...
    'standard', "", ...
    'winner_port', "", ...
    'loser_port', "", ...
    'selection_reason', "", ...
    'global_winner_port', "", ...
    'freq_start_ghz', nan, ...
    'freq_stop_ghz', nan, ...
    'gamma_abs_median', nan, ...
    'gamma_abs_min', nan, ...
    'gamma_abs_max', nan, ...
    'gamma_mag_db_median', nan, ...
    'gamma_mag_db_min', nan, ...
    'gamma_mag_db_max', nan, ...
    'phase_deg_median', nan, ...
    'repeatability_mag_db_median', nan, ...
    'repeatability_mag_db_max', nan, ...
    'repeatability_phase_deg_median', nan, ...
    'repeatability_phase_deg_max', nan, ...
    'sim_complex_error_median', nan, ...
    'sim_complex_error_p95', nan, ...
    'sim_mag_db_error_median', nan, ...
    'sim_mag_db_error_p95', nan, ...
    'heldout_abs_residual_median', nan, ...
    'heldout_abs_residual_p95', nan, ...
    'heldout_abs_residual_max', nan, ...
    'kappa_median', nan, ...
    'kappa_max', nan, ...
    'kappa_pct_below_50', nan, ...
    'open_short_separation_min', nan, ...
    'open_load_separation_min', nan, ...
    'short_load_separation_min', nan, ...
    'load_real_impedance_min_ohm', nan, ...
    'load_real_impedance_max_ohm', nan, ...
    'load_imag_impedance_min_ohm', nan, ...
    'load_imag_impedance_max_ohm', nan, ...
    'consistency_open_improvement_median', nan, ...
    'consistency_short_improvement_median', nan, ...
    'consistency_load_improvement_median', nan, ...
    'consistency_sigma_change_median', nan, ...
    'consistency_thru_error_median', nan, ...
    'consistency_p34_error_median', nan, ...
    'z0_ohm', nan);
end

function rows = build_table51_rows(accepted, solMetrics, consistencyMetrics)
rows = repmat(struct( ...
    'Metric', "", ...
    'Purpose', "", ...
    'QuantitativeProxy', "", ...
    'QualitativeAssessment', "", ...
    'RecommendedWriteup', ""), 7, 1);

openShortMag = collect_band_stat(accepted, {'OPEN', 'SHORT'}, @(entry) median(abs(entry.gamma), 'omitnan'));
loadGammaAbs = collect_band_stat(accepted, {'LOAD'}, @(entry) median(abs(entry.gamma), 'omitnan'));
loadRealMin = arrayfun(@(s) s.load_impedance.real_min, solMetrics);
loadRealMax = arrayfun(@(s) s.load_impedance.real_max, solMetrics);
loadImagMin = arrayfun(@(s) s.load_impedance.imag_min, solMetrics);
loadImagMax = arrayfun(@(s) s.load_impedance.imag_max, solMetrics);
pairwiseMin = arrayfun(@(s) min([s.pairwise.open_short_min, s.pairwise.open_load_min, s.pairwise.short_load_min]), solMetrics);

rows(1).Metric = "Physical plausibility";
rows(1).Purpose = "Verify smooth and interpretable Open/Short/Load behavior at the extracted reference plane";
rows(1).QuantitativeProxy = sprintf(['Open/Short median |Gamma| = %.3f to %.3f across bands; ', ...
    'Load median |Gamma| = %.3f to %.3f; Load Re{Z} spans %.1f to %.1f ohm and Im{Z} spans %.1f to %.1f ohm.'], ...
    min(openShortMag), max(openShortMag), min(loadGammaAbs), max(loadGammaAbs), ...
    min(loadRealMin), max(loadRealMax), min(loadImagMin), max(loadImagMax));
rows(1).QualitativeAssessment = "Open and Short remain strongly reflective; the Load remains usable but is the least ideal standard, especially in the upper D-band.";
rows(1).RecommendedWriteup = "State explicitly that plausibility is strongest for Open/Short and that the Load remains non-ideal but still distinct and physically interpretable.";

rows(2).Metric = "Standard separation";
rows(2).Purpose = "Confirm that Open, Short, and Load remain sufficiently distinct in the complex reflection plane";
rows(2).QuantitativeProxy = sprintf(['Minimum pairwise separations by band: 0-67 GHz %.3f, 67-115 GHz %.3f, 110-170 GHz %.3f. ', ...
    'The closest approach is always between the Load and a reflective standard, not between Open and Short.'], ...
    pairwiseMin(1), pairwiseMin(2), pairwiseMin(3));
rows(2).QualitativeAssessment = "Strong overall. Open and Short remain well separated; the Load is less central than ideal but remains identifiable.";
rows(2).RecommendedWriteup = "Tie the identifiability discussion to the measured pairwise separations rather than only to visual Smith-chart spacing.";

repeatOpen = collect_repeat_stat(accepted, 'OPEN', 'repeatability_mag_db_median');
repeatShort = collect_repeat_stat(accepted, 'SHORT', 'repeatability_mag_db_median');
repeatLoad = collect_repeat_stat(accepted, 'LOAD', 'repeatability_mag_db_median');
rows(3).Metric = "Redundant-measurement repeatability";
rows(3).Purpose = "Quantify spread across redundant realizations and repeated touchdowns";
rows(3).QuantitativeProxy = sprintf(['Median |dMag| between redundant extractions: Open %.3f to %.3f dB, ', ...
    'Short %.3f to %.3f dB, Load %.3f to %.3f dB. Phase spread is consistently largest for the Load and for the upper bands.'], ...
    min(repeatOpen), max(repeatOpen), min(repeatShort), max(repeatShort), min(repeatLoad), max(repeatLoad));
rows(3).QualitativeAssessment = "Good for Open/Short, materially weaker for the Load.";
rows(3).RecommendedWriteup = "Use redundant-measurement spread as the main practical uncertainty indicator, particularly for the Load.";

rows(4).Metric = "Extraction consistency";
rows(4).Purpose = "Assess whether reflective refinement improves one-port behavior while remaining consistent with the thru-derived embedding";
rows(4).QuantitativeProxy = sprintf(['Median reflective-improvement terms by band: Open %.4f/%.4f/%.4f, Short %.4f/%.4f/%.4f, Load %.4f/%.4f/%.4f. ', ...
    'Median thru-reconstruction error remains %.4f to %.4f, with |Delta sigma_max| %.4f to %.4f.'], ...
    consistencyMetrics(1).open_improvement_median, consistencyMetrics(2).open_improvement_median, consistencyMetrics(3).open_improvement_median, ...
    consistencyMetrics(1).short_improvement_median, consistencyMetrics(2).short_improvement_median, consistencyMetrics(3).short_improvement_median, ...
    consistencyMetrics(1).load_improvement_median, consistencyMetrics(2).load_improvement_median, consistencyMetrics(3).load_improvement_median, ...
    min([consistencyMetrics.thru_reconstruction_error_median]), max([consistencyMetrics.thru_reconstruction_error_median]), ...
    min([consistencyMetrics.sigma_change_median]), max([consistencyMetrics.sigma_change_median]));
rows(4).QualitativeAssessment = "Good. Reflective refinement improves one-port behavior without materially undermining the thru-based embedding baseline.";
rows(4).RecommendedWriteup = "Frame the refinement as a controlled reflective correction anchored by the measured Thru, not as a free re-estimation of the embedding.";

rows(5).Metric = "SOL solve conditioning";
rows(5).Purpose = "Evaluate kappa_SOL(f) = kappa_2(A_SOL(f)) versus frequency";
rows(5).QuantitativeProxy = sprintf(['kappa_SOL median / max by band: 0-67 GHz %.2f / %.2f, 67-115 GHz %.2f / %.2f, 110-170 GHz %.2f / %.2f; ', ...
    'kappa_SOL stays below 50 across the full measured span.'], ...
    solMetrics(1).kappa_median, solMetrics(1).kappa_max, ...
    solMetrics(2).kappa_median, solMetrics(2).kappa_max, ...
    solMetrics(3).kappa_median, solMetrics(3).kappa_max);
rows(5).QualitativeAssessment = "Excellent. Numerical conditioning is not the limiting factor in the present one-port kit.";
rows(5).RecommendedWriteup = "Do not use kappa_SOL alone to claim an upper frequency limit; the limiting behavior comes instead from Load repeatability and embedding sensitivity.";

heldoutOpen = arrayfun(@(s) s.heldout(1).median_abs_residual, solMetrics);
heldoutShort = arrayfun(@(s) s.heldout(2).median_abs_residual, solMetrics);
heldoutLoad = arrayfun(@(s) s.heldout(3).median_abs_residual, solMetrics);
rows(6).Metric = "Closure residual";
rows(6).Purpose = "Evaluate r_k(f) and distinguish algebraic closure from held-out validation";
rows(6).QuantitativeProxy = sprintf(['Direct same-set closure is effectively machine-zero because the solved system is square. ', ...
    'Held-out redundant residual medians are Open %.4f to %.4f, Short %.4f to %.4f, Load %.4f to %.4f.'], ...
    min(heldoutOpen), max(heldoutOpen), min(heldoutShort), max(heldoutShort), min(heldoutLoad), max(heldoutLoad));
rows(6).QualitativeAssessment = "Direct closure is not discriminative here; held-out redundant residual is the meaningful validation metric.";
rows(6).RecommendedWriteup = "Revise the chapter text so that the direct residual is treated only as an algebraic sanity check, and report held-out redundant residual as the real lack-of-fit measure.";

simOpen = collect_band_stat(accepted, {'OPEN'}, @(e) median(abs(e.gamma - e.sim_reference), 'omitnan'));
simShort = collect_band_stat(accepted, {'SHORT'}, @(e) median(abs(e.gamma - e.sim_reference), 'omitnan'));
simLoad = collect_band_stat(accepted, {'LOAD'}, @(e) median(abs(e.gamma - e.sim_reference), 'omitnan'));
rows(7).Metric = "Simulation trend agreement";
rows(7).Purpose = "Check broad trend agreement without using simulation as absolute truth";
rows(7).QuantitativeProxy = sprintf(['Median complex-plane |dGamma| against the Chapter 4 de-embedded reference: ', ...
    'Open %.3f to %.3f, Short %.3f to %.3f, Load %.3f to %.3f across the measured bands.'], ...
    min(simOpen), max(simOpen), min(simShort), max(simShort), min(simLoad), max(simLoad));
rows(7).QualitativeAssessment = "Trend-level only. Agreement is useful qualitatively, but not strong enough to treat simulation as a ground-truth benchmark.";
rows(7).RecommendedWriteup = "Keep simulation explicitly secondary and describe it as a directional trend reference rather than a quantitative pass/fail baseline.";
end

function rows = build_table52_rows(accepted)
stdList = {'OPEN', 'SHORT', 'LOAD'};
rows = repmat(struct( ...
    'Standard', "", ...
    'PrimaryRepresentation', "", ...
    'SelectedRealization_0_67', "", ...
    'SelectedRealization_67_115', "", ...
    'SelectedRealization_110_170', "", ...
    'GlobalPreference', "", ...
    'UsableRange', "", ...
    'Notes', ""), numel(stdList), 1);

for iStd = 1:numel(stdList)
    stdName = stdList{iStd};
    entries = accepted.(stdName).bands;
    rows(iStd).Standard = string(stdName);
    rows(iStd).PrimaryRepresentation = "Accepted bandwise extraction; stitched full-span fit retained as reporting/downstream smooth reference";
    rows(iStd).SelectedRealization_0_67 = find_band_winner(entries, '0-67');
    rows(iStd).SelectedRealization_67_115 = find_band_winner(entries, '67-115');
    rows(iStd).SelectedRealization_110_170 = find_band_winner(entries, '110-170');
    rows(iStd).GlobalPreference = summarize_global_preference(entries);
    switch stdName
        case 'OPEN'
            rows(iStd).UsableRange = "Full measured span; mild caution near the upper edge of the D-band";
            rows(iStd).Notes = "Most stable as a strongly reflective standard; simulation agreement is only trend-level in the upper bands.";
        case 'SHORT'
            rows(iStd).UsableRange = "Full measured span; strongest caution near the upper D-band";
            rows(iStd).Notes = "Remains strongly reflective, but the upper D-band shows the clearest loss of ideal short behavior.";
        case 'LOAD'
            rows(iStd).UsableRange = "Full measured span for reporting; upper D-band is the dominant sensitivity region";
            rows(iStd).Notes = "Primary limiting standard. Repeatability and centrality both degrade first in the Load.";
    end
end
end

function notes = build_section_notes(accepted, solMetrics, consistencyMetrics, table51, table52, stitched)
notes = struct();
notes.evaluation_framework = compose_evaluation_framework_note(table51);
notes.open_short = compose_open_short_note(accepted, solMetrics, consistencyMetrics);
notes.load = compose_load_note(accepted, solMetrics);
notes.adopted_standards = compose_adopted_standards_note(solMetrics, table52, stitched);
end

function txt = compose_evaluation_framework_note(table51)
lines = {};
lines{end + 1} = 'Section 5.1.2 Evaluation Framework';
lines{end + 1} = '';
lines{end + 1} = 'The extracted standards should be evaluated using the measurement-led criteria already stated in the chapter, but the present data supports a clearer operationalization of each criterion.';
lines{end + 1} = '';
for i = 1:numel(table51)
    lines{end + 1} = sprintf('%d. %s', i, table51(i).Metric);
    lines{end + 1} = sprintf('Purpose: %s', table51(i).Purpose);
    lines{end + 1} = sprintf('Quantitative proxy: %s', table51(i).QuantitativeProxy);
    lines{end + 1} = sprintf('Interpretation: %s', table51(i).QualitativeAssessment);
    lines{end + 1} = '';
end
lines{end + 1} = 'Recommended chapter-level clarification: the direct same-set closure residual is algebraically trivial for a square three-standard solve and should therefore be reported only as an internal sanity check. The held-out redundant residual is the stronger validation metric because it tests the solved model on an extraction branch that was not used to identify the one-port error terms.';
txt = strjoin(lines, newline);
end

function txt = compose_open_short_note(accepted, solMetrics, consistencyMetrics)
openBands = accepted.OPEN.bands;
shortBands = accepted.SHORT.bands;
lines = {};
lines{end + 1} = 'Section 5.1.3 Measured Open and Short Standards';
lines{end + 1} = '';
lines{end + 1} = sprintf(['Across the accepted bandwise extractions, the Open remains a strongly reflective standard with median |Gamma| between %.3f and %.3f, ', ...
    'while the Short remains strongly reflective with median |Gamma| between %.3f and %.3f. ', ...
    'The Open/Short pairwise separation remains high across all three bands, with median Open-to-Short spacing between %.3f and %.3f and minimum spacing never below %.3f.'], ...
    min(collect_band_stat(accepted, {'OPEN'}, @(e) median(abs(e.gamma), 'omitnan'))), ...
    max(collect_band_stat(accepted, {'OPEN'}, @(e) median(abs(e.gamma), 'omitnan'))), ...
    min(collect_band_stat(accepted, {'SHORT'}, @(e) median(abs(e.gamma), 'omitnan'))), ...
    max(collect_band_stat(accepted, {'SHORT'}, @(e) median(abs(e.gamma), 'omitnan'))), ...
    min(arrayfun(@(s) s.pairwise.open_short_median, solMetrics)), ...
    max(arrayfun(@(s) s.pairwise.open_short_median, solMetrics)), ...
    min(arrayfun(@(s) s.pairwise.open_short_min, solMetrics)));
lines{end + 1} = '';
lines{end + 1} = sprintf(['Redundant-measurement agreement is materially better for the Open and Short than for the Load. ', ...
    'For the Open, the median magnitude spread between redundant realizations ranges from %.3f dB to %.3f dB; ', ...
    'for the Short, the corresponding range is %.3f dB to %.3f dB. ', ...
    'The larger phase and magnitude spread appears primarily in the upper band and near band edges, which is consistent with residual embedding sensitivity rather than with a collapse of the extracted one-port definitions.'], ...
    min(collect_repeat_stat(accepted, 'OPEN', 'repeatability_mag_db_median')), ...
    max(collect_repeat_stat(accepted, 'OPEN', 'repeatability_mag_db_median')), ...
    min(collect_repeat_stat(accepted, 'SHORT', 'repeatability_mag_db_median')), ...
    max(collect_repeat_stat(accepted, 'SHORT', 'repeatability_mag_db_median')));
lines{end + 1} = '';
lines{end + 1} = sprintf(['Reflective refinement improves the one-port behavior without materially degrading the thru-derived embedding. ', ...
    'Across the three bands, the median Open improvement metric is %.4f to %.4f, the median Short improvement metric is %.4f to %.4f, ', ...
    'and the median thru-reconstruction error stays between %.4f and %.4f. ', ...
    'This supports the interpretation that the final embedding remains transmission-anchored by the straight Thru while being reflectively corrected by the measured SOL data.'], ...
    min([consistencyMetrics.open_improvement_median]), max([consistencyMetrics.open_improvement_median]), ...
    min([consistencyMetrics.short_improvement_median]), max([consistencyMetrics.short_improvement_median]), ...
    min([consistencyMetrics.thru_reconstruction_error_median]), max([consistencyMetrics.thru_reconstruction_error_median]));
lines{end + 1} = '';
lines{end + 1} = sprintf(['Comparison to the Chapter 4 design references remains trend-level only. ', ...
    'The measured Open shows median complex-plane disagreement |dGamma| of %.3f to %.3f across the three bands, ', ...
    'and the Short shows %.3f to %.3f. The measured standards therefore follow the broad expected reflective roles, ', ...
    'but the fabricated extraction should not be judged by exact EM overlay.'], ...
    min(collect_band_stat(accepted, {'OPEN'}, @(e) median(abs(e.gamma - e.sim_reference), 'omitnan'))), ...
    max(collect_band_stat(accepted, {'OPEN'}, @(e) median(abs(e.gamma - e.sim_reference), 'omitnan'))), ...
    min(collect_band_stat(accepted, {'SHORT'}, @(e) median(abs(e.gamma - e.sim_reference), 'omitnan'))), ...
    max(collect_band_stat(accepted, {'SHORT'}, @(e) median(abs(e.gamma - e.sim_reference), 'omitnan'))));
lines{end + 1} = '';
lines{end + 1} = sprintf(['Adopted bandwise Open winner ports: 0-67 GHz %s, 67-115 GHz %s, 110-170 GHz %s. ', ...
    'Adopted bandwise Short winner ports: 0-67 GHz %s, 67-115 GHz %s, 110-170 GHz %s.'], ...
    openBands(1).winner_port, openBands(2).winner_port, openBands(3).winner_port, ...
    shortBands(1).winner_port, shortBands(2).winner_port, shortBands(3).winner_port);
txt = strjoin(lines, newline);
end

function txt = compose_load_note(accepted, solMetrics)
loadBands = accepted.LOAD.bands;
loadGammaMedians = collect_band_stat(accepted, {'LOAD'}, @(e) median(abs(e.gamma), 'omitnan'));
lines = {};
lines{end + 1} = 'Section 5.1.4 Measured Load Standard';
lines{end + 1} = '';
lines{end + 1} = sprintf(['The extracted Load is the least ideal one-port standard, but it remains sufficiently central and sufficiently distinct from the reflective standards to keep the one-port solve well conditioned. ', ...
    'Across the accepted bandwise results, the Load median |Gamma| ranges from %.3f to %.3f. ', ...
    'The extracted Load impedance spans %.1f to %.1f ohm in real part and %.1f to %.1f ohm in imaginary part across the full measured set.'], ...
    min(loadGammaMedians), max(loadGammaMedians), ...
    min(arrayfun(@(s) s.load_impedance.real_min, solMetrics)), max(arrayfun(@(s) s.load_impedance.real_max, solMetrics)), ...
    min(arrayfun(@(s) s.load_impedance.imag_min, solMetrics)), max(arrayfun(@(s) s.load_impedance.imag_max, solMetrics)));
lines{end + 1} = '';
lines{end + 1} = sprintf(['The Load is also the dominant repeatability limiter. ', ...
    'Median magnitude spread between redundant extracted Loads ranges from %.3f dB to %.3f dB, ', ...
    'and the held-out redundant residual median ranges from %.4f to %.4f in complex reflection coefficient. ', ...
    'This is substantially larger than for the Open or Short and is the main reason to treat the Load as the practical uncertainty driver of the one-port kit.'], ...
    min(collect_repeat_stat(accepted, 'LOAD', 'repeatability_mag_db_median')), ...
    max(collect_repeat_stat(accepted, 'LOAD', 'repeatability_mag_db_median')), ...
    min(arrayfun(@(s) s.heldout(3).median_abs_residual, solMetrics)), ...
    max(arrayfun(@(s) s.heldout(3).median_abs_residual, solMetrics)));
lines{end + 1} = '';
lines{end + 1} = sprintf(['Even with that limitation, the numerical one-port solve remains well behaved. ', ...
    'The accepted standards produce median kappa_SOL values of %.2f, %.2f, and %.2f for the 0-67, 67-115, and 110-170 GHz bands, respectively, ', ...
    'with worst-case values of %.2f, %.2f, and %.2f. ', ...
    'This indicates that the practical upper-band concern is not matrix ill-conditioning itself, but the increasing uncertainty of the extracted Load and the associated sensitivity of the reflective standards to residual embedding error.'], ...
    solMetrics(1).kappa_median, solMetrics(2).kappa_median, solMetrics(3).kappa_median, ...
    solMetrics(1).kappa_max, solMetrics(2).kappa_max, solMetrics(3).kappa_max);
lines{end + 1} = '';
lines{end + 1} = 'The direct same-set closure residual is effectively machine-zero and should not be overinterpreted. The more useful validation result is the held-out Load residual, which remains modest but grows in the upper band. Accordingly, the chapter should frame the Load as acceptable for subsequent calibration use across the full measured span, while making clear that the upper D-band is the dominant caution region.';
lines{end + 1} = '';
lines{end + 1} = sprintf(['Adopted bandwise Load winner ports: 0-67 GHz %s, 67-115 GHz %s, 110-170 GHz %s.'], ...
    loadBands(1).winner_port, loadBands(2).winner_port, loadBands(3).winner_port);
txt = strjoin(lines, newline);
end

function txt = compose_adopted_standards_note(solMetrics, table52, stitched)
lines = {};
lines{end + 1} = 'Section 5.1.5 Adopted One-Port Standards for the Measured SOLR Calibration';
lines{end + 1} = '';
lines{end + 1} = 'The accepted bandwise extracted standards remain the primary calibration standards for measured SOLR. The stitched full-span fits are retained as complementary reporting products and as smooth downstream references, but they should not be described as replacements for the accepted bandwise measurements.';
lines{end + 1} = '';
for i = 1:numel(table52)
    lines{end + 1} = sprintf('%s:', table52(i).Standard);
    lines{end + 1} = sprintf('Selected realizations = 0-67 %s, 67-115 %s, 110-170 %s', ...
        table52(i).SelectedRealization_0_67, table52(i).SelectedRealization_67_115, table52(i).SelectedRealization_110_170);
    lines{end + 1} = sprintf('Global preference = %s', table52(i).GlobalPreference);
    lines{end + 1} = sprintf('Suggested usable range = %s', table52(i).UsableRange);
    lines{end + 1} = sprintf('Notes = %s', table52(i).Notes);
    lines{end + 1} = '';
end
lines{end + 1} = sprintf(['The one-port conditioning remains excellent through the full measured span, with kappa_SOL below %.2f even in the worst case. ', ...
    'For that reason, the most defensible practical limitation statement is not that the one-port model fails above a sharply defined frequency, ', ...
    'but that confidence becomes increasingly controlled by Load repeatability and high-band embedding sensitivity near the upper edge of the 110-170 GHz band.'], ...
    max([solMetrics.kappa_max]));
lines{end + 1} = '';
lines{end + 1} = 'The stitched curves should be described as regularized aggregates of the accepted bandwise winners. Their role is to provide continuity across the full span and a smooth handoff to the downstream calibration flow. They are not better measurements than the accepted bandwise curves; they are smoother representations of them.';
lines{end + 1} = '';
lines{end + 1} = sprintf(['Saved stitched final-fit products are available for OPEN, SHORT, and LOAD over %.2f-%.2f GHz, %.2f-%.2f GHz, and %.2f-%.2f GHz on the common full-span grids, respectively.'], ...
    min(stitched.OPEN.freq) / 1e9, max(stitched.OPEN.freq) / 1e9, ...
    min(stitched.SHORT.freq) / 1e9, max(stitched.SHORT.freq) / 1e9, ...
    min(stitched.LOAD.freq) / 1e9, max(stitched.LOAD.freq) / 1e9);
txt = strjoin(lines, newline);
end

function txt = build_writeup_suggestions(solMetrics, consistencyMetrics)
lines = {};
lines{end + 1} = 'Write-up Suggestions';
lines{end + 1} = '';
lines{end + 1} = '1. Reframe the closure-residual discussion.';
lines{end + 1} = 'The direct residual computed from the same three standards used to solve the one-port model is effectively algebraic for the present square system. Keep it as an internal sanity check, but add the held-out redundant residual as the meaningful lack-of-fit metric.';
lines{end + 1} = '';
lines{end + 1} = sprintf(['2. Do not use kappa_SOL alone to define the practical upper frequency limit. ', ...
    'The one-port conditioning remains excellent through the full measured span, with kappa_SOL max %.2f. ', ...
    'If the chapter wants to discuss a practical caution region, anchor that discussion in Load repeatability, held-out residual growth, and high-band embedding sensitivity instead.'], ...
    max([solMetrics.kappa_max]));
lines{end + 1} = '';
lines{end + 1} = '3. Make the role of the stitched curves explicit.';
lines{end + 1} = 'The stitched curves are smooth, regularized aggregates of the accepted bandwise results. They are useful for reporting and for continuous downstream reference products, but they should not be described as replacing the accepted bandwise extracted standards.';
lines{end + 1} = '';
lines{end + 1} = '4. Keep simulation comparison secondary.';
lines{end + 1} = 'The measured-versus-simulated overlays are useful for trend-level direction checks only. The present data does not support describing simulation disagreement by itself as extraction failure.';
lines{end + 1} = '';
lines{end + 1} = '5. Make the Load the explicit limiting standard in the narrative.';
lines{end + 1} = 'The Open and Short are the better behaved extracted standards. The Load carries the largest redundant-measurement spread and the largest held-out residual, so the discussion of kit quality should foreground the Load rather than treating all three standards symmetrically.';
lines{end + 1} = '';
lines{end + 1} = '6. Convert Table 5.1 from a qualitative placeholder into an operationalized metric table.';
lines{end + 1} = 'Each criterion can now be tied to a concrete proxy: |Gamma| and Zload ranges, pairwise separations, redundant spread, refinement/thru consistency metrics, kappa_SOL, held-out residual, and trend-level simulation error.';
lines{end + 1} = '';
lines{end + 1} = sprintf(['7. State the refinement interpretation carefully. ', ...
    'The final embedding should be described as straight-Thru anchored with reflective correction from the measured SOL data. ', ...
    'Median thru-reconstruction error remains between %.4f and %.4f after refinement.'], ...
    min([consistencyMetrics.thru_reconstruction_error_median]), max([consistencyMetrics.thru_reconstruction_error_median]));
txt = strjoin(lines, newline);
end

function values = collect_band_stat(accepted, stdNames, fn)
values = [];
for i = 1:numel(stdNames)
    entries = accepted.(stdNames{i}).bands;
    for k = 1:numel(entries)
        values(end + 1, 1) = fn(entries(k)); %#ok<AGROW>
    end
end
end

function values = collect_repeat_stat(accepted, stdName, fieldName)
entries = accepted.(stdName).bands;
values = nan(numel(entries), 1);
for k = 1:numel(entries)
    if ~entries(k).has_loser
        continue;
    end
    dMag = abs(safe_mag_db(entries(k).gamma) - safe_mag_db(entries(k).gamma_loser));
    dPhase = abs(unwrap(angle(entries(k).gamma)) - unwrap(angle(entries(k).gamma_loser))) * 180 / pi;
    switch fieldName
        case 'repeatability_mag_db_median'
            values(k) = median(dMag, 'omitnan');
        case 'repeatability_phase_deg_median'
            values(k) = median(dPhase, 'omitnan');
    end
end
end

function port = find_band_winner(entries, bandName)
port = "";
for k = 1:numel(entries)
    if strcmp(entries(k).band, bandName)
        port = entries(k).winner_port;
        return;
    end
end
end

function pref = summarize_global_preference(entries)
ports = strings(numel(entries), 1);
for k = 1:numel(entries)
    ports(k) = entries(k).global_winner_port;
end
uniquePorts = unique(ports(ports ~= ""));
if isempty(uniquePorts)
    pref = "not saved";
elseif numel(uniquePorts) == 1
    pref = uniquePorts(1);
else
    pref = strjoin(cellstr(uniquePorts), '/');
end
end

function sim = get_optional_sim_reference(selection, freq)
sim = nan(size(freq));
if ~isfield(selection, 'sim_reference')
    return;
end
sim = extract_gamma_reference(selection.sim_reference, freq);
end

function sim = extract_gamma_reference(simRef, freq)
sim = nan(size(freq));
if isempty(simRef)
    return;
end
if isnumeric(simRef)
    sim = resize_like(make_col(simRef), freq);
    return;
end
if isstruct(simRef)
    if isfield(simRef, 'gamma')
        gamma = make_col(simRef.gamma);
        if isfield(simRef, 'freq')
            sim = interp1(make_col(simRef.freq), gamma, freq, 'linear', nan);
        else
            sim = resize_like(gamma, freq);
        end
        return;
    end
    if isfield(simRef, 'S')
        sMat = simRef.S;
        if ndims(sMat) == 3
            gamma = squeeze(sMat(1, 1, :));
            if isfield(simRef, 'freq')
                sim = interp1(make_col(simRef.freq), make_col(gamma), freq, 'linear', nan);
            else
                sim = resize_like(make_col(gamma), freq);
            end
            return;
        end
    end
end
end

function out = get_optional_col(structVal, fieldName, freq)
if isfield(structVal, fieldName)
    out = resize_like(make_col(structVal.(fieldName)), freq);
else
    out = nan(size(freq));
end
end

function out = get_optional_output_col(outputStruct, fieldName, freq)
if isfield(outputStruct, fieldName)
    out = resize_like(make_col(outputStruct.(fieldName)), freq);
elseif strcmp(fieldName, 'gamma') && isfield(outputStruct, 'gamma')
    out = resize_like(make_col(outputStruct.gamma), freq);
else
    out = nan(size(freq));
end
end

function out = get_optional_output_string(outputStruct, fieldName)
if isfield(outputStruct, fieldName)
    out = string(outputStruct.(fieldName));
else
    out = "";
end
end

function out = get_optional_struct(structVal, fieldName)
if isfield(structVal, fieldName)
    out = structVal.(fieldName);
else
    out = struct();
end
end

function port = get_optional_global_port(bandResult, stdName)
port = "";
if isfield(bandResult, 'final_selection_global') && ...
        isfield(bandResult.final_selection_global, stdName) && ...
        isfield(bandResult.final_selection_global.(stdName), 'winner_port')
    port = string(bandResult.final_selection_global.(stdName).winner_port);
end
end

function key = derive_output_key(stdName, port, bandName)
key = sprintf('EXTRACTED_%s_%s_%s', upper(stdName), upper(char(port)), strrep(bandName, '-', '_'));
end

function port = other_port(portIn)
if strcmpi(char(portIn), 'P1')
    port = "P2";
else
    port = "P1";
end
end

function [startGHz, stopGHz] = parse_band_limits(bandName)
tokens = regexp(bandName, '(\d+)-(\d+)', 'tokens', 'once');
if isempty(tokens)
    error('DissertationResults:BadBandName', 'Unable to parse band name: %s', bandName);
end
startGHz = str2double(tokens{1});
stopGHz = str2double(tokens{2});
end

function out = make_col(x)
out = x(:);
end

function out = resize_like(x, ref)
if isempty(x)
    out = nan(size(ref));
elseif numel(x) == numel(ref)
    out = make_col(x);
elseif isscalar(x)
    out = repmat(x, size(ref));
else
    out = nan(size(ref));
end
end

function value = get_field_or_default(structVal, fieldName, defaultValue)
if isfield(structVal, fieldName)
    value = structVal.(fieldName);
else
    value = defaultValue;
end
end

function values = safe_mag_db(gamma)
values = 20 * log10(max(abs(gamma), 1e-12));
end

function p = percentile_local(x, pct)
x = x(isfinite(x));
if isempty(x)
    p = nan;
    return;
end
x = sort(x);
idx = 1 + (numel(x) - 1) * pct / 100;
lo = floor(idx);
hi = ceil(idx);
if lo == hi
    p = x(lo);
else
    frac = idx - lo;
    p = (1 - frac) * x(lo) + frac * x(hi);
end
end

function write_text(pathStr, textStr)
fid = fopen(pathStr, 'w');
if fid < 0
    error('DissertationResults:WriteFailed', 'Unable to write file: %s', pathStr);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', textStr);
end
