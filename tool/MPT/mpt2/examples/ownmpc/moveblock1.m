clear sysStruct probStruct
opt_sincos

echo on
% Move blocking for control problems based on hybrid systems
%
% set finite prediction horizon
probStruct.N = 5;

% only consider two free control moves
probStruct.Nc = 2;

% compute an on-line controller
ctrl = mpt_control(sysStruct, probStruct, 'online');

% obtain open-loop trajectory
x0 = [5; 5];
[X, U] = sim(ctrl, x0, struct('openloop', 1));
U

echo off