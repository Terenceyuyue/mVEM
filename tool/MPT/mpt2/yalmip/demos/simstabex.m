yalmip('clear');
echo on
clc
%*********************************************************
%
% Simultaneuous stabilization using non-convex 
% semidefinite programming.
%
%*********************************************************
%
% This example solves a simultaneuous stabilization problem
%
% The example is taken from
%  "A path-following method for solving BMI problems in control"
%  A. Hassibi, J. How, and S. Boyd. 
%  Proceedings of American Control Conference, 2:1385-1389, June 1999. 
%
%
% The goal is to find a state-feedback, stabilizing a set
% of systems, under  constraints on the feedback matrix.
% 
% The performance measure is to maximize the smallest 
% decay-rate for the closed-loop systems.
%
% The optimization problem can be stated as:
%
%  max t
%  subject to  (A1+B*K)'*P1+P1*(A1+B*K) < -2*t*P1
%              (A2+B*K)'*P2+P2*(A2+B*K) < -2*t*P2
%              (A3+B*K)'*P3+P3*(A3+B*K) < -2*t*P3
%                                   -50 < K
%                                     K < 50
%                              P1,P2,P3 > 0
%
% Clearly, this is a BMI due to the products between Pi 
% and K, and the products t*Pi
pause
clc

% Generate system data
A1 = [1 -1 0;1 1 0;0 0 -0.5];
A2 = [1.5 -7 0;7 1.5 0;0 0 1];
A3 = [-0.5 -3 0;3 -0.5 0;0 0 2];
B = [0.2477 -0.1645;0.4070 0.8115;0.6481 0.4083];
pause

% An initial guess is often important when solving BMIs
%
% The following set of LMIs gives us a stabilizing controller and
% a reasonable initial guess
% (let P1=P2=P3=P, multply BMIs with Q=P^-1, let  L=KP^-1, t=0 and 
%  neglect amplitude constraints)
L = sdpvar(2,3);
Q = sdpvar(3,3);
F = set(Q>0) + set(trace(Q)==1); % Normalize
F = F + set(Q*A1'+(B*L)'+(A1*Q+B*L) < 0);
F = F + set(Q*A2'+(B*L)'+(A2*Q+B*L) < 0);
F = F + set(Q*A3'+(B*L)'+(A3*Q+B*L) < 0);
pause
% Find a feasible solution using penbmi as convex SDP solver
solvesdp(F,-trace(Q),sdpsettings('solver','penbmi'));

% The solution Q and L can be used later to define initial guesses for P and K
P0 = inv(double(Q));
K0 = double(L)*P0;
pause
clc
% Define variables in actual BMI-problem
t = sdpvar(1,1);
K = sdpvar(2,3);
P1 = sdpvar(3,3);
P2 = sdpvar(3,3);
P3 = sdpvar(3,3);

% BMI problem (P>0 changed to P>I for numerical reasons)
F = set(P1>eye(3)) + set(P2>eye(3))+set(P3>eye(3));
F = F + set(-50 < K < 50);
F = F + set((A1+B*K)'*P1+P1*(A1+B*K) < -2*t*P1);
F = F + set((A2+B*K)'*P2+P2*(A2+B*K) < -2*t*P2);
F = F + set((A3+B*K)'*P3+P3*(A3+B*K) < -2*t*P3);
pause

% Solve with the calculated initial guesses (and t = 0 initially)
setsdpvar(P1,P0);
setsdpvar(P2,P0);
setsdpvar(P3,P0);
setsdpvar(K,K0);
setsdpvar(t,0);
solvesdp(F,-t,sdpsettings('usex0',1));

% Decay-rate obtained
double(t)
% Feedback matrix
double(K)
pause
echo off