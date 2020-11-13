clc
echo on
yalmip('clear')
%*********************************************************
%
% Modeling and solving KYP problems
%
%*********************************************************
%
% This example requires the KYPD-solver.
%
% We will show how we can use the dedicated
% KYPD solver to solve a L2-gain analysis problem.
%
% minimize t
%  s.t [A'P+PA+C'C PB;B'P -t] < 0
%                           P > 0
% 
pause

% Define a random stable system
n = 30;
A = randn(n);A = A - max(real(eig(A)))*eye(n)*1.5;
B = randn(n,1);
C = randn(1,n);
pause

% YALMIP variables
t = sdpvar(1,1);
P = sdpvar(n,n);
pause

% KYP constraint 
% NOTE : Do not add the constraint P>0
% if you want to use the KYPD solver)
F = set(kyp(A,B,P,blkdiag(C'*C,-t)) < 0)
pause

% Solve using dedicated solver
sol1 = solvesdp(F,t,sdpsettings('solver','kypd'));
pause

% Solve using SeDuMi (or any other standard SDP solver)
sol2 = solvesdp(F,t,sdpsettings('solver','sedumi'));
pause

% Compare solution time (this can be a bit mis-leadinhg due to memory caching etc...)
% KYPD should typically be a lot faster than the straight-forward solution
sol1.solvertime/sol2.solvertime

pause
echo off