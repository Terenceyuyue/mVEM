clc
echo on
%*********************************************************
%
% Decay-rate estimation using non-convex SDP
%
%*********************************************************
%
% The problem we solve is estimatation of decay-rate of a 
% linear system x' = Ax. This can be formulated as a 
% generalized eigenvalue problem (GEVP)
%
% max alpha
% s.t A'P+PA < -2alphaP
%          P > I
%
% This time, we solve it as a BMI using PENBMI (hence you 
% need PENBMI to run this demo). Note, this is a quasi-convex
% problem, and PENBMI is actually guaranteed to find the
% global optima.
pause
clc

% Define the variables
A = [-1 2;-3 -4];
P = sdpvar(2,2);
alpha = sdpvar(1,1);
pause

% Define the GEVP
F = set(P>eye(2))+set(A'*P+P*A < -2*alpha*P) + set(alpha > 0);
pause
% Maximize alpha (minimize -alpha)
solvesdp(F,-alpha);
double(alpha)
pause