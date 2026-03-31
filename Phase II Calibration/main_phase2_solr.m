function results = main_phase2_solr()
%MAIN_PHASE2_SOLR Entry point for downstream bandwise four-port SOLR.

clearvars;
close all force;
clc;

scriptDir = fileparts(mfilename('fullpath'));
addpath(scriptDir);

results = run_phase2_solr();
end
