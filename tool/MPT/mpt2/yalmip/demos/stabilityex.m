clc
echo on
%*********************************************************
%
% Lyapunov analysis using semidefinite programming
%
%*********************************************************

% Create a stable matrix
A = [-1 2;0 -2];
pause % Strike any key to continue. 

% Create symmetric matrix (full syntax)
P = sdpvar(2,2,'symmetric'); 
pause % Strike any key to continue. 

% Add SETs for stability
F = set(P>eye(2)) + set(A'*P+P*A<0)
pause % Strike any key to continue. 

% Find feasible solution, minimize trace(P)
solution = solvesdp(F,trace(P));
pause % Strike any key to continue. 

% Extract numerical solution
P_feasible = double(P);

% Check solution
eig(P_feasible)
eig(A'*P_feasible+P_feasible'*A)
pause % Strike any key to continue. 

% Checking the constraints can also be done with
checkset(F)

pause % Strike any key to continue. 
echo off