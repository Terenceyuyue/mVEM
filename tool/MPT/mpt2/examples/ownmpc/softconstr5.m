clear sysStruct probStruct
Double_Integrator
close all

echo on
% Demonstration of soft constraints
%
% Turn off stabilizing target set
probStruct.Tconstraint = 0;

% Set prediction horizon
probStruct.N = 5;

% Set output constraints
sysStruct.ymax = [5; 5];
sysStruct.ymin = [-5; -5];

% Soft constraints can also be used together with mpt_ownmpc():
probStruct.Sy = 100;         % soften both output constraints
probStruct.symax = [1; 1];   % maximum allowed violation of each constraint

% Prepare constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct, 'online');

% Now we can add a logical which says that if the constraint on the first
% element of the first predicted output (y_0) is violated, u_0 must be above 0.5
% and u_1 must be between -0.7 and -0.4:
C = C + set(iff(V.sy{1}(1) > 0, V.u{1} > 0.5));
C = C + set(iff(V.sy{1}(1) > 0, -0.7 < V.u{2} < -0.4));

% Compute an online controller
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V);

% Compute open-loop optimizer
x0 = [5.5; 0];
u = ctrl(x0, struct('openloop', 1))

% Notice that the open-loop optimizer indeed satisfies given logic conditions,
% since constraints on y_0 are violated (y_0 = x_0 in this example), therefore
% the respective slack variable V.sy{1} is above zero.

echo off
