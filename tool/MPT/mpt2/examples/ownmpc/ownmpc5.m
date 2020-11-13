Double_Integrator
close all

echo on

% Switch off stabilizing target set
probStruct.Tconstraint = 0;

% Prediction horizon 4
probStruct.N = 3;

% Construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct, 'online');

% Add a constraint that the sum of all predicted control actions along the
% prediction horizon should be equal to zero 

C = C + set(sum([V.u{:}]) == 0);

% Compute an on-line controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V, 'online');

% Simulate the open-loop system
x0 = [4; 0];
[X, U] = sim(ctrl, x0, struct('openloop', 1)); U

% Verify that the open-loop control actions sum up to zero (with some tolerance)
sum(U)

echo off
