clear sysStruct probStruct Options
close all
echo on

% import a nonlinear system
sysStruct = mpt_sys(@duffing_oscillator);

% penalize quadratic cost objective
probStruct.norm = 2;
probStruct.Q = eye(2);
probStruct.R = 1/10;

% set prediction horizon
probStruct.N = 3;

% calculate an on-line controller
ctrl = mpt_control(sysStruct, probStruct, 'online');

% simulate the closed-loop system for a period of 10 seconds
x0 = [2.5; 1];
Tfinal = 3;
simsteps = Tfinal / sysStruct.Ts;

% For safety, monitor what is happening in YALMIP
Options.verbose = 2;

% use a local nonlinear solver
Options.nlsolver = 'local';

% simulate and plot the closed-loop trajectory
simplot(ctrl, x0, simsteps, Options);

figure

% use YALMIPs global nonlinear solver
Options.nlsolver = 'global';

% use fmincon as the upper bound solver (set this option to 'none' if you don't
% have fmincon installed)
Options.uppersolver = 'fmincon';

% with a quadratic objective, a QP lower bound solver is recommended
% (testing indicates that CLP is prefered). However, since CLP only is
% available in the Windows distribution, we use the LP solver CDD instead
Options.lowersolver = 'cdd';

% Stop the branch-and-bound tree early
Options.nliter = 15;

% This takes longer time but aims at a globally optimal solution
simplot(ctrl, x0, simsteps, Options);

echo off
