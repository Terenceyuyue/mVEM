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

% Add the polytopic constraint on each predicted state excluding the initial
% state x0
C = C + set(H * V.x{2} <= K);
C = C + set(H * V.x{3} <= K);

% Compute an explicit controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V);

% Plot the controller partition. Notice that the feasible set is now NOT
% constrained by H*x0<=K, since we have added said constraint on x1 and x2.
plot(ctrl);

% Simulate the open-loop system
x0 = [4.5; -1.5];
X = sim(ctrl, x0, struct('openloop', 1))

% Notice that 2nd and 3rd element of the open-loop trajectory indeed satisfy
% given polytopic constraints

echo off
