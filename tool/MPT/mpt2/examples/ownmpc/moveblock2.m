clear sysStruct probStruct
opt_sincos

echo on
% Move blocking can of course be used together with mpt_ownmpc() and/or with
% soft constraints
%
% set finite prediction horizon
probStruct.N = 2;

% use linear cost function
probStruct.norm = 1;

% only consider one free control move
probStruct.Nc = 1;

% consider output constraints as soft constraints, allow a maximum violation of
% 1.5:
probStruct.symax = [1.5; 1.5];
probStruct.Sy = diag([10 100]);   % penalize violations of hard constraints

% construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct);

% add custom term to the objective function
O = O + norm(V.sy{1}-2, 1);

% calculate an explicit controller
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V);

% plot the controller partition
plot(ctrl);

% calculate a closed-loop trajectory
x0 = [5; 5]; simsteps = 20;
figure;
simplot(ctrl, x0, simsteps);

echo off