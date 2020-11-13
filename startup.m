%% Add_local_path

% Include all folders under the current path and the path is no longer valid after restarting MATLAB.
% Adding savepath will permanently add the current path to the MATLAB search path;
% However, it is not recommended to do so in order to avoid naming conflicts.
% You can open pathdef.m to see the added paths.


clc;clear;close all
warning off
addpath(genpath(pwd)); 
% savepath; % not recommended
% open pathdef.m % check the path

