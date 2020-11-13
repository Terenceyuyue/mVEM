clear sysStruct probStruct
Double_Integrator
close all

echo on
% Demonstration of soft constraints
%
% Turn off stabilizing target set
probStruct.Tconstraint = 0;

% Set prediction horizon
probStruct.N = 3;

% Set output constraints
sysStruct.ymax = [5; 5];
sysStruct.ymin = [-5; -5];

% Soft constraints can also be used together with mpt_ownmpc():
probStruct.Sy = 100;         % soften both output constraints
probStruct.symax = [5; 5];   % maximum allowed violation of each constraint

% Prepare constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct);

% Notice that the structure "V" now contains the variables "sy" which correspond
% to slack variables used to soften respective constraints. Therefore we can add
% arbitrary constraints based on these variables. Say we want to limit the
% maximum violation of the output constaints at the first prediction step by +1:
C = C + set(V.sy{1} <= 1);

% Compute an explicit controller
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V);

% Plot the controller partition
plot(ctrl);

echo off
