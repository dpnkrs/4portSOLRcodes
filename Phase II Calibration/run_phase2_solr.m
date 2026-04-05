function results = run_phase2_solr(config)
%RUN_PHASE2_SOLR Orchestrate downstream Phase II bandwise SOLR calibration.

if nargin < 1 || isempty(config)
    config = phase2_config();
end
ensure_output_dirs(config);
apply_plot_defaults();

inventory = discover_phase2_measurements(config);
referenceStandards = load_reference_standards(config);
generate_phase2_reference_import_debug(config, inventory, referenceStandards);

results = struct();
results.config = config;
results.inventory = inventory;
results.references = referenceStandards;
results.bands = cell(1, numel(config.bands));

for idxBand = 1:numel(config.bands)
    bandName = config.bands{idxBand};
    bandKey = config.band_keys{idxBand};
    fprintf('Phase II processing band %s\n', bandName);
    bandData = preload_band_data(inventory, bandKey);
    bandResult = solve_phase2_band(config, bandName, bandData, referenceStandards);
    results.bands{idxBand} = bandResult;

    save(fullfile(config.mat_dir, sprintf('PHASE2_BAND_%s.mat', bandName)), 'bandResult');
end

results.stitched = build_stitched_outputs(config, results.bands);
generate_phase2_switch_ratio_debug(config, results.bands);
generate_phase2_pair_waterfall_debug(config, results.bands);
generate_phase2_error_term_debug(config, results.bands);
generate_phase2_port_reflection_debug(config, results.bands, referenceStandards);
generate_phase2_transmission_scale_audit(config, results.bands, referenceStandards);
generate_phase2_thru_target_consistency_debug(config, results.bands);
generate_phase2_midband_arc_debug(config, results.bands);
save_phase2_outputs(config, results);
write_phase2_walkthrough(config, results);

save(fullfile(config.mat_dir, 'PHASE2_RESULTS_ALL.mat'), 'results');
fprintf('Phase II calibration complete.\n');
fprintf('Outputs: %s\n', config.output_dir);
end

function bandResult = solve_phase2_band(config, bandName, bandData, referenceStandards)
pairResults = solve_pair_solr_band(config, config.reference_pairs(1), bandName, bandData, referenceStandards);
for idxPair = 2:numel(config.reference_pairs)
    pairResults(idxPair) = solve_pair_solr_band(config, config.reference_pairs(idxPair), bandName, bandData, referenceStandards);
end

holdCfg = config.holdout_pairs(1);
referenceIdx = find(strcmpi({config.reference_pairs.pair}, holdCfg.pair), 1, 'first');
holdNet = bandData.holdout_thrus.(holdCfg.pair).(upper(holdCfg.geometry));
holdoutResults = correct_pair_network(config, holdNet, pairResults(referenceIdx), holdCfg.pair, upper(holdCfg.geometry), true);
for idxHold = 2:numel(config.holdout_pairs)
    holdCfg = config.holdout_pairs(idxHold);
    referenceIdx = find(strcmpi({config.reference_pairs.pair}, holdCfg.pair), 1, 'first');
    holdNet = bandData.holdout_thrus.(holdCfg.pair).(upper(holdCfg.geometry));
    holdoutResults(idxHold) = correct_pair_network(config, holdNet, pairResults(referenceIdx), holdCfg.pair, upper(holdCfg.geometry), true);
end

pairBlocks = correct_pair_network(config, bandData.dut.LangeCoupler, pairResults(1), pairResults(1).pair, pairResults(1).geometry, false);
for idxPair = 2:numel(pairResults)
    pairBlocks(idxPair) = correct_pair_network(config, bandData.dut.LangeCoupler, pairResults(idxPair), pairResults(idxPair).pair, pairResults(idxPair).geometry, false);
end
assembledDUT = assemble_corrected_4port_from_pairs(config, pairBlocks);
dutSummary = summarize_dut_metrics(assembledDUT.S);

solResidualSummary = summarize_sol_residuals(pairResults);

bandResult = struct();
bandResult.band = bandName;
bandResult.freq = pairResults(1).freq;
bandResult.reference_pair_results = pairResults;
bandResult.holdout_results = holdoutResults;
bandResult.sol_residual_summary = solResidualSummary;
bandResult.dut = struct( ...
    'file', bandData.dut.LangeCoupler.file, ...
    'raw', bandData.dut.LangeCoupler, ...
    'pair_corrected_blocks', pairBlocks, ...
    'corrected', assembledDUT, ...
    'summary', dutSummary);

plot_phase2_band_results(config, bandResult);
write_band_summary_note(config, bandResult);
end

function corrected = correct_pair_network(config, network, pairModel, pairKey, geometry, isHoldout)
ports = phase2_pair_key_to_ports(pairKey);
Sraw = phase2_extract_pair_submatrix(network.S, pairKey);
switchTerms = pairModel.switch_terms;

if numel(network.freq) ~= numel(pairModel.freq) || any(abs(network.freq(:) - pairModel.freq(:)) > 1)
    switchTerms = struct();
    switchTerms.gamma_forward = interp1(pairModel.freq, pairModel.switch_terms.gamma_forward, network.freq, 'linear', 'extrap');
    switchTerms.gamma_reverse = interp1(pairModel.freq, pairModel.switch_terms.gamma_reverse, network.freq, 'linear', 'extrap');
    errorTerms = interpolate_error_terms(pairModel.error_terms, pairModel.freq, network.freq);
else
    errorTerms = pairModel.error_terms;
end

Sswitch = apply_switch_correction_2port(Sraw, switchTerms.gamma_forward, switchTerms.gamma_reverse, config.switch_den_floor);
Scorrected = nan(size(Sswitch));
for idxFreq = 1:numel(network.freq)
    Scorrected(:, :, idxFreq) = apply_error_correction_8term_local(Sswitch(:, :, idxFreq), errorTerms, idxFreq, config.cal_den_floor);
end

metrics = struct();
metrics.reciprocity_mismatch = abs(squeeze(Scorrected(2, 1, :)) - squeeze(Scorrected(1, 2, :)));
metrics.return_loss = min(-20 * log10(max(abs(squeeze(Scorrected(1, 1, :))), 1e-12)), ...
    -20 * log10(max(abs(squeeze(Scorrected(2, 2, :))), 1e-12)));
metrics.insertion_loss_db = -20 * log10(max(abs(squeeze(Scorrected(2, 1, :))), 1e-12));
metrics.phase_s21_deg = unwrap(angle(squeeze(Scorrected(2, 1, :)))) * 180 / pi;
metrics.phase_s12_deg = unwrap(angle(squeeze(Scorrected(1, 2, :)))) * 180 / pi;

corrected = struct();
corrected.file = network.file;
corrected.pair = pairKey;
corrected.geometry = geometry;
corrected.ports = ports;
corrected.freq = network.freq(:);
corrected.S_raw = Sraw;
corrected.S_switch = Sswitch;
corrected.S_corrected = Scorrected;
corrected.corrected = Scorrected;
corrected.metrics = metrics;
corrected.is_holdout = isHoldout;
end

function summary = summarize_sol_residuals(pairResults)
stdNames = {'Open', 'Short', 'Load'};
summary = struct();
for idxStd = 1:numel(stdNames)
    rows = [];
    for idxPair = 1:numel(pairResults)
        rows = [rows; pairResults(idxPair).standard_residuals.(stdNames{idxStd})]; %#ok<AGROW>
    end
    summary.(stdNames{idxStd}) = rows;
end
end

function metrics = summarize_dut_metrics(S)
nFreq = size(S, 3);
offDiag = [];
recip = [];
for row = 1:4
    for col = 1:4
        if row == col
            continue;
        end
        offDiag = [offDiag; squeeze(S(row, col, :)).']; %#ok<AGROW>
        if row < col
            recip = [recip; abs(squeeze(S(row, col, :)).' - squeeze(S(col, row, :)).')]; %#ok<AGROW>
        end
    end
end

diagVals = zeros(4, nFreq);
for idxPort = 1:4
    diagVals(idxPort, :) = squeeze(S(idxPort, idxPort, :)).';
end

port1 = squeeze(S(2:4, 1, :));
medianPort1 = median(-20 * log10(max(abs(port1), 1e-12)), 2, 'omitnan');
sortedMedian = sort(medianPort1);

metrics = struct();
metrics.median_return_loss_db = median(-20 * log10(max(abs(diagVals(:)), 1e-12)), 'omitnan');
metrics.median_reciprocity_mismatch = median(recip(:), 'omitnan');
metrics.max_reciprocity_mismatch = max(recip(:), [], 'omitnan');
metrics.isolation_db = max(median(-20 * log10(max(abs(offDiag), 1e-12)), 2, 'omitnan'));
if numel(sortedMedian) >= 2
    metrics.coupling_balance_db = abs(sortedMedian(1) - sortedMedian(2));
else
    metrics.coupling_balance_db = NaN;
end
magDb = 20 * log10(max(abs(offDiag), 1e-12));
metrics.smoothness_rms_db = rms(diff(magDb, 2, 2), 'all');
end

function stitched = build_stitched_outputs(config, bandResults)
stitched = struct();

dutNetworks = cell(1, numel(bandResults));
for idxBand = 1:numel(bandResults)
    dutNetworks{idxBand} = bandResults{idxBand}.dut.corrected;
end
stitched.dut = stitch_phase2_networks(dutNetworks, config.stitch_interp_method);

for idxHold = 1:numel(config.holdout_pairs)
    pairKey = config.holdout_pairs(idxHold).pair;
    geomKey = upper(config.holdout_pairs(idxHold).geometry);
    holdNetworks = cell(1, numel(bandResults));
    for idxBand = 1:numel(bandResults)
        hit = find(strcmpi({bandResults{idxBand}.holdout_results.pair}, pairKey), 1, 'first');
        holdNetworks{idxBand} = struct('freq', bandResults{idxBand}.holdout_results(hit).freq, ...
            'S', bandResults{idxBand}.holdout_results(hit).S_corrected, 'z0', config.z0);
    end
    stitched.holdout.(pairKey).(geomKey) = stitch_phase2_networks(holdNetworks, config.stitch_interp_method);
end
end

function save_phase2_outputs(config, results)
write_touchstone_nport(fullfile(config.touchstone_dir, 'CORRECTED_LangeCoupler_STITCHED.s4p'), ...
    results.stitched.dut.freq, results.stitched.dut.S, config.z0);
save(fullfile(config.mat_dir, 'CORRECTED_LangeCoupler_STITCHED.mat'), 'results');

for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    write_touchstone_nport(fullfile(config.touchstone_dir, sprintf('CORRECTED_LangeCoupler_%s.s4p', bandResult.band)), ...
        bandResult.dut.corrected.freq, bandResult.dut.corrected.S, config.z0);
    save(fullfile(config.mat_dir, sprintf('CORRECTED_LangeCoupler_%s.mat', bandResult.band)), 'bandResult');

    for idxHold = 1:numel(bandResult.holdout_results)
        holdout = bandResult.holdout_results(idxHold);
        write_touchstone_nport(fullfile(config.touchstone_dir, sprintf('CORRECTED_%s_%s_%s.s2p', holdout.pair, holdout.geometry, bandResult.band)), ...
            holdout.freq, holdout.S_corrected, config.z0);
    end
end

for idxHold = 1:numel(config.holdout_pairs)
    pairKey = config.holdout_pairs(idxHold).pair;
    geomKey = upper(config.holdout_pairs(idxHold).geometry);
    net = results.stitched.holdout.(pairKey).(geomKey);
    write_touchstone_nport(fullfile(config.touchstone_dir, sprintf('CORRECTED_%s_%s_STITCHED.s2p', pairKey, geomKey)), ...
        net.freq, net.S, config.z0);
end

plot_stitched_results(config, results);
write_summary_note(config, results);
end

function plot_stitched_results(config, results)
fig = figure('Visible', config.figure_visible, 'Color', 'w');
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
freqGHz = results.stitched.dut.freq / 1e9;
colorMap = lines(max(6, numel(config.holdout_pairs) + 3));

nexttile;
hold on;
for idxPort = 2:4
    sVal = squeeze(results.stitched.dut.S(idxPort, 1, :));
    plot(freqGHz, 20 * log10(max(abs(sVal), 1e-12)), 'LineWidth', 1.3, ...
        'Color', colorMap(idxPort - 1, :), 'DisplayName', sprintf('S%d1', idxPort));
end
grid on;
xlabel('Frequency (GHz)');
ylabel('|S_{x1}| (dB)');
title('Stitched corrected Lange coupler transmission');
legend('Location', 'best');

nexttile;
hold on;
for idxHold = 1:numel(config.holdout_pairs)
    pairKey = config.holdout_pairs(idxHold).pair;
    geomKey = upper(config.holdout_pairs(idxHold).geometry);
    net = results.stitched.holdout.(pairKey).(geomKey);
    s21 = squeeze(net.S(2, 1, :));
    plot(freqGHz, 20 * log10(max(abs(interp1(net.freq, s21, results.stitched.dut.freq, 'linear', 'extrap')), 1e-12)), ...
        'LineWidth', 1.2, 'Color', colorMap(idxHold, :), 'DisplayName', sprintf('%s %s', pairKey, geomKey));
end
grid on;
xlabel('Frequency (GHz)');
ylabel('|S_{21}| (dB)');
title('Stitched corrected DIAG hold-out transmission');
legend('Location', 'bestoutside');

save_figure(fig, fullfile(config.figure_dir, 'PhaseII_Stitched_DUT_And_Holdouts'));
end

function write_summary_note(config, results)
lines = {};
lines{end + 1} = 'Phase II Summary';
lines{end + 1} = '===============';
lines{end + 1} = '';
for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    lines{end + 1} = sprintf('Band %s', bandResult.band);
    lines{end + 1} = sprintf('  Median DUT return loss: %.3f dB', bandResult.dut.summary.median_return_loss_db);
    lines{end + 1} = sprintf('  Median DUT reciprocity mismatch: %.5f', bandResult.dut.summary.median_reciprocity_mismatch);
    lines{end + 1} = sprintf('  Coupling balance surrogate: %.3f dB', bandResult.dut.summary.coupling_balance_db);
    lines{end + 1} = sprintf('  Isolation surrogate: %.3f dB', bandResult.dut.summary.isolation_db);
    for idxHold = 1:numel(bandResult.holdout_results)
        holdout = bandResult.holdout_results(idxHold);
        lines{end + 1} = sprintf('  Holdout %s %s mean |S21-S12|: %.5f', ...
            holdout.pair, holdout.geometry, mean(holdout.metrics.reciprocity_mismatch, 'omitnan'));
    end
    lines{end + 1} = '';
end
write_text_file(fullfile(config.note_dir, 'PhaseII_Summary.txt'), lines);
end

function write_band_summary_note(config, bandResult)
lines = {};
lines{end + 1} = sprintf('Band %s summary', bandResult.band);
lines{end + 1} = sprintf('Median DUT return loss: %.3f dB', bandResult.dut.summary.median_return_loss_db);
lines{end + 1} = sprintf('Median DUT reciprocity mismatch: %.5f', bandResult.dut.summary.median_reciprocity_mismatch);
lines{end + 1} = sprintf('Isolation surrogate: %.3f dB', bandResult.dut.summary.isolation_db);
for idxPair = 1:numel(bandResult.reference_pair_results)
    pairResult = bandResult.reference_pair_results(idxPair);
    lines{end + 1} = sprintf('%s %s mean reciprocity mismatch: %.5f', ...
        pairResult.pair, pairResult.geometry, pairResult.reference_metrics.mean_recip_mismatch);
end
write_text_file(fullfile(config.note_dir, sprintf('PhaseII_BandSummary_%s.txt', bandResult.band)), lines);
end

function referenceStandards = load_reference_standards(config)
stdNames = fieldnames(config.reference_standard_files);
referenceStandards = struct();
for idx = 1:numel(stdNames)
    stdName = stdNames{idx};
    switch lower(config.reference_standard_mode)
        case 'stitched_final'
            loaded = load(config.reference_standard_files.(stdName));
            fitOut = loaded.fitOut;
            referenceStandards.(stdName) = struct( ...
                'mode', 'stitched_final', ...
                'freq', fitOut.freq(:), ...
                'gamma', fitOut.gamma_fit(:));
        case 'bandwise_final'
            bandwise = struct();
            for idxBand = 1:numel(config.bands)
                bandName = config.bands{idxBand};
                loaded = load(fullfile(config.phase1_output_mat_dir, sprintf('FINALFIT_%s_%s.mat', upper(stdName), bandName)));
                finalResult = loaded.finalResult;
                bandwise.(matlab.lang.makeValidName(strrep(bandName, '-', '_'))) = struct( ...
                    'freq', finalResult.freq(:), ...
                    'gamma', finalResult.(config.reference_bandwise_field)(:), ...
                    'winner_port', finalResult.winner_port, ...
                    'winner_output_key', finalResult.winner_output_key);
            end
            referenceStandards.(stdName) = struct('mode', 'bandwise_final', 'bands', bandwise);
        otherwise
            error('Unsupported reference_standard_mode: %s', config.reference_standard_mode);
    end
end
end

function bandData = preload_band_data(inventory, bandKey)
bandData = struct();
bandData.standards = struct();
stdNames = fieldnames(inventory.raw_standards.(bandKey));
for idxStd = 1:numel(stdNames)
    stdName = stdNames{idxStd};
    ports = fieldnames(inventory.raw_standards.(bandKey).(stdName));
    for idxPort = 1:numel(ports)
        bandData.standards.(stdName).(ports{idxPort}) = read_touchstone_nport(inventory.raw_standards.(bandKey).(stdName).(ports{idxPort}));
    end
end

bandData.reference_thrus = struct();
refPairs = fieldnames(inventory.reference_thrus.(bandKey));
for idxPair = 1:numel(refPairs)
    pairKey = refPairs{idxPair};
    geoms = fieldnames(inventory.reference_thrus.(bandKey).(pairKey));
    for idxGeom = 1:numel(geoms)
        geomKey = geoms{idxGeom};
        bandData.reference_thrus.(pairKey).(geomKey) = read_touchstone_nport(inventory.reference_thrus.(bandKey).(pairKey).(geomKey));
    end
end

bandData.switch_terms = inventory.switch_terms.(bandKey);
bandData.dut = struct('LangeCoupler', read_touchstone_nport(inventory.duts.(bandKey).LangeCoupler));
bandData.holdout_thrus = struct();
holdPairs = fieldnames(inventory.holdout_thrus.(bandKey));
for idxPair = 1:numel(holdPairs)
    pairKey = holdPairs{idxPair};
    geoms = fieldnames(inventory.holdout_thrus.(bandKey).(pairKey));
    for idxGeom = 1:numel(geoms)
        geomKey = geoms{idxGeom};
        bandData.holdout_thrus.(pairKey).(geomKey) = read_touchstone_nport(inventory.holdout_thrus.(bandKey).(pairKey).(geomKey));
    end
end
end

function termsInterp = interpolate_error_terms(errorTerms, freqIn, freqOut)
fields = fieldnames(errorTerms);
termsInterp = struct();
for idx = 1:numel(fields)
    values = errorTerms.(fields{idx});
    realPart = interp1(freqIn, real(values), freqOut, 'linear', 'extrap');
    imagPart = interp1(freqIn, imag(values), freqOut, 'linear', 'extrap');
    termsInterp.(fields{idx}) = complex(realPart, imagPart);
end
end

function ensure_output_dirs(config)
dirs = {config.output_dir, config.touchstone_dir, config.mat_dir, config.figure_dir, config.note_dir, ...
    config.interim_dir, config.interim_mat_dir, config.interim_figure_dir, config.interim_note_dir};
for idx = 1:numel(dirs)
    if ~exist(dirs{idx}, 'dir')
        mkdir(dirs{idx});
    end
end
end

function apply_plot_defaults()
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultTextColor', 'k');
set(groot, 'defaultAxesColorOrder', lines(12));
set(groot, 'defaultAxesLineStyleOrder', '-');
set(groot, 'defaultLineLineWidth', 1.2);
end

function save_figure(fig, basePath)
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
