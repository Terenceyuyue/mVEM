Double_Integrator
close all

echo on

% Switch off stabilizing target set
probStruct.Tconstraint = 0;

% Prediction horizon 4
probStruct.N = 3;

% Construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct, 'online');

% Add a constraint that sum of absolute values of all control inputs should be
% less than some bounds

C = C + set(sum(abs([V.u{:}])) <= 1.5);

% Compute an on-line controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V, 'online');

% Simulate the open-loop system
x0 = [4; 0];
[X, U] = sim(ctrl, x0, struct('openloop', 1)); U

% Verify that the sum of absolute values is bounded by 1.5
sum(abs(U))

echo off
