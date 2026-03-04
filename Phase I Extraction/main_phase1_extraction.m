% main_phase1_extraction.m
% Top-level script for Phase I calibration standard extraction.

clear variables;
clc;
close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

results = run_phase1_extraction();

disp('Phase I extraction finished.');
disp('Processed bands:');
disp(cellfun(@(band) band.band, results.bands, 'UniformOutput', false));
