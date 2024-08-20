function init(APGAA_path)
% Initialization process
clearvars -except APGAA_path;
close all;
clc;

% Adding sub-folder to MATLAB path
addpath(fullfile(APGAA_path,'examples'))
addpath(fullfile(APGAA_path,'panel_tools'))
addpath(fullfile(APGAA_path,'solver'))