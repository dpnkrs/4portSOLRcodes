function write_phase2_walkthrough(config, results)
%WRITE_PHASE2_WALKTHROUGH Write implementation-facing notes for Phase II.

lines = {};
lines{end + 1} = 'Phase II Calibration Walkthrough';
lines{end + 1} = '==============================';
lines{end + 1} = '';
lines{end + 1} = sprintf('Root folder: %s', config.script_dir);
lines{end + 1} = sprintf('Raw standards folder: %s', config.phase2_raw_dir);
lines{end + 1} = sprintf('Switch terms folder: %s', config.phase2_switch_dir);
lines{end + 1} = sprintf('DUT folder: %s', config.phase2_dut_dir);
lines{end + 1} = '';
lines{end + 1} = 'Processing order';
lines{end + 1} = '----------------';
lines{end + 1} = '1. Discover Phase II raw standards, switch-term CSVs, and DUTs by band.';
lines{end + 1} = '2. Load stitched Phase I Open/Short/Load references from the Phase I outputs folder.';
lines{end + 1} = '3. For each reference pair, extract the active-pair 2x2 submatrix from the measured .s4p.';
lines{end + 1} = '4. Parse the four wave CSV files for that pair and compute forward/reverse switch ratios as b/a at the receiver sampler.';
lines{end + 1} = '5. Apply switch correction to the active-pair through and the relevant one-port standard submatrices.';
lines{end + 1} = '6. Solve an 8-term SOLR model for that pair using the stitched Phase I reflection standards as the downstream references.';
lines{end + 1} = '7. Correct the active-pair submatrices of the DUT using the pair model, then assemble the corrected 4-port DUT by averaging diagonal estimates and inserting pair-specific off-diagonal terms.';
lines{end + 1} = '8. Correct the hold-out DIAG structures with the same pair model for validation.';
lines{end + 1} = '9. Stitch corrected DUT and hold-out outputs across the three bands for reporting.';
lines{end + 1} = '';
lines{end + 1} = 'Important implementation details';
lines{end + 1} = '-------------------------------';
lines{end + 1} = 'Reciprocal-thru files are treated as active-pair only. The named pair submatrix is used; the remaining 4-port entries are ignored.';
lines{end + 1} = 'The orthogonal ARC measurements are the calibration references. The DIAG measurements from PhaseII DUTs are held out for validation.';
lines{end + 1} = 'No GLOBALWIN references are used in Phase II. Only FINALFIT_*_STITCHED_FINAL.mat from Phase I are used as the standard definitions.';
lines{end + 1} = 'DIAG hold-outs do not have dedicated switch-wave files in the provided dataset; they are corrected with the switch ratios from the corresponding ARC pair solve.';
lines{end + 1} = '';
lines{end + 1} = 'Files to inspect first';
lines{end + 1} = '---------------------';
lines{end + 1} = sprintf('1. %s', fullfile(config.mat_dir, 'PHASE2_RESULTS_ALL.mat'));
lines{end + 1} = sprintf('2. %s', fullfile(config.note_dir, 'PhaseII_Summary.txt'));
lines{end + 1} = sprintf('3. %s', fullfile(config.figure_dir, 'PhaseII_DUT_Diagnostics_67-115.jpg'));
lines{end + 1} = sprintf('4. %s', fullfile(config.touchstone_dir, 'CORRECTED_LangeCoupler_STITCHED.s4p'));
lines{end + 1} = '';

for idxBand = 1:numel(results.bands)
    bandResult = results.bands{idxBand};
    refMetricValues = arrayfun(@(x) x.reference_metrics.mean_recip_mismatch, bandResult.reference_pair_results);
    lines{end + 1} = sprintf('Band %s', bandResult.band);
    lines{end + 1} = sprintf('  DUT file: %s', bandResult.dut.file);
    lines{end + 1} = sprintf('  Median reciprocity mismatch of used references: %.4g', median(refMetricValues, 'omitnan'));
    lines{end + 1} = sprintf('  Hold-out count: %d', numel(bandResult.holdout_results));
    lines{end + 1} = '';
end

write_text_file(config.walkthrough_filename, lines);
write_text_file(fullfile(config.note_dir, 'PhaseII_Walkthrough.txt'), lines);
end

function write_text_file(filePath, lines)
[outDir, ~, ~] = fileparts(filePath);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end
fid = fopen(filePath, 'w');
cleanup = onCleanup(@() fclose(fid));
for idx = 1:numel(lines)
    fprintf(fid, '%s\n', lines{idx});
end
end
