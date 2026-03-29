% main_dissertation_results.m
% Top-level driver for dissertation-facing Phase I result summaries.

clear variables;
clc;
close all;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

summary = run_dissertation_results();
%generate_figure5_2(summary);
%generate_figure5_3(summary);
generate_figure5_4(summary);

disp('Dissertation result module finished.');
disp('Generated outputs in:');
disp(summary.output_root);
