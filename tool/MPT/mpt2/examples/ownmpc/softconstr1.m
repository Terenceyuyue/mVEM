clear sysStruct probStruct
Double_Integrator
close all

echo on
% Demonstration of soft constraints
%
% Enforce stabilizing target set
probStruct.Tconstraint = 1;

% Set prediction horizon
probStruct.N = 3;

% Set output constraints
sysStruct.ymax = [5; 5];
sysStruct.ymin = [-5; -5];

% Our system has 2 outputs with constraints set on both of them. Now we would
% like to consider both output constraints as soft. To do so, we define penalty
% on slack variables, which will turn on soft constraints:
probStruct.Sy = 100;

% When using soft constraints, you should always set sysStruct.Pbnd to limit the
% exploration space.
sysStruct.Pbnd = unitbox(2, 20);

% Now calculate an explicit controller
ctrl = mpt_control(sysStruct, probStruct);

% Plot the controller partition
plot(ctrl);

% The original hard constraints on system outputs have been set to +/- 5. By
% using soft constraints, we are allowing the constraints to be violated, as can
% be seen from the figure. If you wonder why there are white spaces on the plot,
% the explanation is easy -- we have enforced a stabilizing target set
% constraint, i.e. the last predicted state must lie in the orange set in the
% center. Since this is still a hard constraint which was not softened, there
% are portions of the state-space for which no feasible control law exists,
% despite the softened output constraints.

echo off
