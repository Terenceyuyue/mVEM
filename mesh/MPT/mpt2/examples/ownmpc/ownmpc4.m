Double_Integrator
close all

echo on

% Switch off stabilizing target set
probStruct.Tconstraint = 0;

% Prediction horizon 4
probStruct.N = 3;

% Construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct, 'online');

% We want to add a "mixed" constraint of the form:
%    y_k(1) + u_k <= 3
%    y_k(1) + u_k >= -3
%
% where y_k(1) denotes the first element of the output vector 
% predicted at time "k" (our system has 2 outputs)

for k = 1:length(V.u),
    C = C + set(-3 <= V.y{k}(1) + V.u{k} <= 3);
end

% Compute an on-line controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V, 'online');

% Simulate the closed-loop system
x0 = [4; 0];
[X, U, Y] = sim(ctrl, x0, 5); U, Y

% Verify that constraints are satisfied (each element must be between -3 and 3)
UplusY = U + Y(:, 1)

echo off
