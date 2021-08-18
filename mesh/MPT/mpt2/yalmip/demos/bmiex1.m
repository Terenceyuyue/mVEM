echo on
clc
% This example shows how to solve a simple BMI using PENBMI
%
% The problem is to find a feedback u = Lx so that the 
% L2 gain from w to y is minimized, for the system
%
% x' = Ax+Bu+Gw
% y  = Cx
%
% This can be formulated as the BMI
%
% min t
% s.t                                 P > 0
%     [(A+BL)P+P(A+BL)+C'C PG; G'P -tI] < 0
%
pause
clc

% Create system data
A = [-1    -1    -1;
     1     0     0;
     0     1     0];
B = [1;0;0];
C = [0 0 1];
G = [-1;-1;-1];

pause

% Define decision variables
P = sdpvar(3,3);
L = sdpvar(1,3);
t = sdpvar(1,1);
pause

% A reasonble initial guess is valuable
[L0,P0]=lqr(A,B,eye(3),1);
setsdpvar(P,P0);
setsdpvar(L,-L0);
setsdpvar(t,100);
pause
clc

% Define the BMI
F = set([(A+B*L)'*P+P*(A+B*L)+C'*C P*G;G'*P -t]<0) + set(P>0)
pause

% Solve, using the initial guess on P, L and t
solvesdp(F,t,sdpsettings('solver','penbmi','usex0',1))
double(L)
pause
echo off