function [xnext, yk] = di_sim_fun(xk, uk)
% DI_SIM_FUN Example of a function to model arbitrary dynamics
%
% If you want to use the sim() or simplot() functions to simulate controllers
% with different dynamical systems than those which have been used to compute
% the controller, you can do so by providing a handle to a special simulation
% function when calling sim() or simplot(), e.g.:
%
%   sim(ctrl, @di_sim_fun, x0)
%
% The function handle must point to a function which takes exactly two inputs
% and produces exactly two outputs. The first input is always the state at time
% "k" (xk) and input at time "k" (uk). First output argument must be the next
% state update (xnext) and output associated to state xk (yk). For example:
%
%     function [xnext, yk] = your_sim_fun(xk, uk)
%
%     A = [1 1; 0 1]; B = [1; 0.5]; C = [1 0];
%     xnext = A*xk + B*uk;
%     yk    = C*xk;
%
% You can use arbitrary dynamics, i.e. also non-linear functions.

A = [1 1; 0 1];
B = [1; 0.5];
C = [1 0];
xnext = A*xk + B*uk;
yk = C*xk;
