Double_Integrator
close all

echo on

% Switch off stabilizing target set
probStruct.Tconstraint = 0;

% Prediction horizon 5
probStruct.N = 5;

% Use linear norm
probStruct.norm = 1;

% Construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct, 'online');

% Logic constraints:
% we now add a constraint which will tell the controller that u_1 must be above
% 0.5 if an only if u_0 is below zero.

C = C + set(iff(V.u{1} <= 0, V.u{2} >= 0.5));

% Compute an explicit controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V, 'online');

% Simulate the open-loop system
x0 = [1; 0];
[X, U] = sim(ctrl, x0, struct('openloop', 1));
U

% Notice that since U(1) (which corresponds to u_0) is below zero, U(2) (which
% corresponds to u_1) is correcly bounded from below by 0.5.
%
% Now we change the initial state:
x0 = [-1; 0];
[X, U] = sim(ctrl, x0, struct('openloop', 1));
U

% Now since U(1) is >= 0, U(2) will be less than 0.5 (by the "if and only if"
% implication).

echo off
