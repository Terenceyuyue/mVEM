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

% Our system has 2 outputs with constraints set on both of them. Now we would
% like to to consider both constraints as soft, where the maximum violation of
% these constraints can be at most 1.
probStruct.Sy = 100;
probStruct.symax = [1; 1];

% If you use probStruct.symax, it is no longer necessary to set sysStruct.Pbnd,
% since the feasible set will be automatically bounded.
%
% Now calculate an explicit controller
ctrl = mpt_control(sysStruct, probStruct);

% Plot the controller partition
plot(ctrl);

% Notice that the output constraints originally set to +/- 5 are now violated by
% at most 1, where each violation is penalized by the value of probStruct.Sy.

echo off
