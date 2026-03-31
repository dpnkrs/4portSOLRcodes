function generate_phase2_reference_import_debug(config, inventory, referenceStandards)
%GENERATE_PHASE2_REFERENCE_IMPORT_DEBUG Plot imported Phase I references against raw downstream standards.

portLabels = config.port_labels;
stdNames = config.standard_names;
bandEdgesGHz = [67, 115];
debugData = struct();
noteLines = {
    'Phase II reference-import sanity check'
    '===================================='
    ''
    };

for idxPort = 1:numel(portLabels)
    portName = portLabels{idxPort};
    portIdx = idxPort;
    fig = figure('Visible', config.figure_visible, 'Color', 'w', 'Position', [100, 100, 1280, 1080]);
    tiledlayout(fig, numel(stdNames), 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    sgtitle(fig, sprintf('Phase II imported references vs raw measured standards - %s', portName), ...
        'FontWeight', 'bold');

    debugData.(portName) = struct();
    for idxStd = 1:numel(stdNames)
        stdName = stdNames{idxStd};
        [stitchedFreqGHz, stitchedMagDb, stitchedPhaseDeg] = local_prepare_stitched_trace(referenceStandards.(stdName));

        axMag = nexttile((idxStd - 1) * 2 + 1);
        hold(axMag, 'on');
        plot(axMag, stitchedFreqGHz, stitchedMagDb, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Phase I stitched reference');
        local_add_band_markers(axMag, bandEdgesGHz);
        ylabel(axMag, sprintf('%s |\\Gamma| (dB)', stdName));
        title(axMag, sprintf('%s magnitude', stdName));
        grid(axMag, 'on');
        xlim(axMag, [0, 170]);

        axPh = nexttile((idxStd - 1) * 2 + 2);
        hold(axPh, 'on');
        plot(axPh, stitchedFreqGHz, stitchedPhaseDeg, 'k-', 'LineWidth', 1.8, 'DisplayName', 'Phase I stitched reference');
        local_add_band_markers(axPh, bandEdgesGHz);
        ylabel(axPh, sprintf('%s phase (deg)', stdName));
        title(axPh, sprintf('%s phase', stdName));
        grid(axPh, 'on');
        xlim(axPh, [0, 170]);

        standardData = struct();
        for idxBand = 1:numel(config.bands)
            bandName = config.bands{idxBand};
            bandKey = config.band_keys{idxBand};
            rawFile = inventory.raw_standards.(bandKey).(stdName).(portName);
            rawNet = read_touchstone_nport(rawFile);
            gammaMeasured = squeeze(rawNet.S(portIdx, portIdx, :));
            gammaMeasured = gammaMeasured(:);
            gammaReference = interp1(referenceStandards.(stdName).freq(:), referenceStandards.(stdName).gamma(:), ...
                rawNet.freq(:), 'pchip', 'extrap');
            gammaReference = gammaReference(:);

            freqGHz = rawNet.freq(:) / 1e9;
            magMeasuredDb = 20 * log10(max(abs(gammaMeasured), 1e-12));
            magReferenceDb = 20 * log10(max(abs(gammaReference), 1e-12));
            phaseReferenceDeg = unwrap(angle(gammaReference)) * 180 / pi;
            phaseMeasuredDeg = unwrap(angle(gammaMeasured)) * 180 / pi;
            phaseMeasuredDeg = local_align_phase_branch(phaseMeasuredDeg, phaseReferenceDeg);

            magDiffDb = magMeasuredDb - magReferenceDb;
            phaseDiffDeg = phaseMeasuredDeg - phaseReferenceDeg;

            standardData.(bandKey) = struct( ...
                'freq_hz', rawNet.freq(:), ...
                'gamma_measured', gammaMeasured, ...
                'gamma_reference', gammaReference, ...
                'mag_measured_db', magMeasuredDb, ...
                'mag_reference_db', magReferenceDb, ...
                'phase_measured_deg', phaseMeasuredDeg, ...
                'phase_reference_deg', phaseReferenceDeg, ...
                'mag_diff_db', magDiffDb, ...
                'phase_diff_deg', phaseDiffDeg, ...
                'raw_file', rawFile);

            plot(axMag, freqGHz, magMeasuredDb, 'LineWidth', 1.2, 'DisplayName', sprintf('%s raw measured', bandName));
            plot(axPh, freqGHz, phaseMeasuredDeg, 'LineWidth', 1.2, 'DisplayName', sprintf('%s raw measured', bandName));

            noteLines{end + 1, 1} = sprintf('%s %s %s: median |dMag| = %.3f dB, median |dPhase| = %.3f deg', ...
                portName, stdName, bandName, median(abs(magDiffDb), 'omitnan'), median(abs(phaseDiffDeg), 'omitnan'));
        end

        if idxStd == numel(stdNames)
            xlabel(axMag, 'Frequency (GHz)');
            xlabel(axPh, 'Frequency (GHz)');
        end
        legend(axMag, 'Location', 'best');
        legend(axPh, 'Location', 'best');

        debugData.(portName).(stdName) = standardData;
    end

    local_save_figure(fig, fullfile(config.interim_figure_dir, sprintf('PhaseII_Debug_RefVsMeasured_%s', portName)));
end

save(fullfile(config.interim_mat_dir, 'PHASE2_DEBUG_REF_IMPORT.mat'), 'debugData');
local_write_text_file(fullfile(config.interim_note_dir, 'PhaseII_Debug_RefImport_Summary.txt'), noteLines);
end

function [freqGHz, magDb, phaseDeg] = local_prepare_stitched_trace(reference)
freqGHz = reference.freq(:) / 1e9;
magDb = 20 * log10(max(abs(reference.gamma(:)), 1e-12));
phaseDeg = unwrap(angle(reference.gamma(:))) * 180 / pi;
end

function phaseAligned = local_align_phase_branch(phaseMeasured, phaseReference)
phaseAligned = phaseMeasured;
offsetCycles = round(median((phaseReference - phaseMeasured) / 360, 'omitnan'));
phaseAligned = phaseAligned + 360 * offsetCycles;
end

function local_add_band_markers(ax, bandEdgesGHz)
for idxEdge = 1:numel(bandEdgesGHz)
    xline(ax, bandEdgesGHz(idxEdge), '--', 'Color', [0.6, 0.6, 0.6], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end
end

function local_save_figure(fig, basePath)
set(fig, 'Color', 'w');
ax = findall(fig, 'type', 'axes');
set(ax, 'Color', 'w', 'XColor', 'k', 'YColor', 'k');
drawnow;
exportgraphics(fig, [basePath, '.jpg'], 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(fig, [basePath, '.pdf'], 'BackgroundColor', 'white', 'ContentType', 'vector');
savefig(fig, [basePath, '.fig']);
close(fig);
end

function local_write_text_file(filePath, lines)
fid = fopen(filePath, 'w');
cleanup = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end
