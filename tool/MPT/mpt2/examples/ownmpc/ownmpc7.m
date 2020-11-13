Double_Integrator
close all

echo on

% Switch off stabilizing target set
probStruct.Tconstraint = 0;

% Prediction horizon 4
probStruct.N = 3;

% Use linear norm
probStruct.norm = 1;

% Construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct, 'online');

% Add contraction constraints, i.e. force state x(k+1) to be closer to the
% origin (in 1/Inf-norm sense) than the state x(k) has been 
for k = 1:length(V.x)-1,
    C = C + set(norm(V.x{k+1}, 1) <= norm(V.x{k}, 1));
end

% Compute an on-line controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V, 'online');

% Simulate the closed-loop system for 5 steps
x0 = [4; 0];
[X, U] = sim(ctrl, x0, 5); X

% Verify that the contraction constraint is satisfied
for k = 1:size(X, 1)-1,
    norm(X(k+1, :), 1) <= norm(X(k, :), 1)
end

echo off
