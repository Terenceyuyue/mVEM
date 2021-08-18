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

% Here we want to consider the constraint on second output as a hard constraint,
% while we allow the first output constraint to be violated by at most 2.
% Moreover, we want to penalize every such violation by a penalty of 1000.
probStruct.Sy = diag([1000 0]);   % no penalty on violation of 2nd output constraint
probStruct.symax = [2; 0];        % a zero indicates that a respective constraint is hard

% Now calculate an explicit controller
ctrl = mpt_control(sysStruct, probStruct);

% Plot the controller partition
plot(ctrl);

% Notice that now only the constraint on second system output is softened.

echo off
