Double_Integrator
close all

echo on

% Switch off stabilizing target set
probStruct.Tconstraint = 0;

% Use shorter prediction horizon => less complex solution
probStruct.N = 2;

% Construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct);

% Matrices of polytopic constraints H*x<=K
H = [eye(2); -eye(2)];
K = [2; 2; 2; 2];

% Add the polytopic constraint on each predicted state
C = C + set(H * V.x{1} <= K);
C = C + set(H * V.x{2} <= K);
C = C + set(H * V.x{3} <= K);

% Compute an explicit controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V);

% Plot the controller partition. Notice that the feasible set is constrained by
% the given polytopic constraint, since said constraint was enforced on x0
% (V.x{1}) as well
plot(ctrl);

echo off
