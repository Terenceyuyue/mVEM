% "Turbo Car" demo

% Copyright (C) 2005 by Michal Kvasnica (kvasnica@control.ee.ethz.ch)

clear sysStruct probStruct

% "turbocar" is basically a double integrator with one exception - it has one
% additional boolean input, which denotes the "turbo". if this input is 1, influence
% of u(1) in the state update equation is doubled. there is a hard limit on how
% many times the "turbo" can be used - this is given by bounds on the 3rd state
% variable. every time the "turbo" is used, x(3) is decremented by 1. lower
% limits is hence set to 0, upper limit can vary between 0 and 5 (only integer
% values are allowed).

% task of the controller is to drive the position to 0 while penalizing usage of
% the "turbo"

% import dynamical model from HYSDEL source
sysStruct = mpt_sys(which('turbo_car.hys'));

% state constraints
sysStruct.xmax = [20; 5; 5];    % limits on position, speed and "maximum turbo counter"
sysStruct.xmin = [-20; -5; 0];

% input constraints
sysStruct.umax = [1; 1];
sysStruct.umin = [-1; 0];

% output constraints
sysStruct.ymax = 20;
sysStruct.ymin = -20;

% penalty on states
probStruct.Q = eye(3);

% penalty on inputs
probStruct.R = [0.1 0; 0 3];   % note that we penalize usage of turbo quite heavily


%penalty on final state
probStruct.P_N = diag([1 1 0]);

% prediction horizon
probStruct.N = 5;

% Inf-norm type cost function
probStruct.norm = Inf;

% compute on-line constroller
onlc = mpt_control(sysStruct, probStruct, 'online');

% simulate the closed loop for 20 sampling instances starting from x0=[20;0;5]
mpt_plotTimeTrajectory(onlc,[20;0;5],20)

open_system('turbocarsim');
sim('turbocarsim', 20);