clc
echo on
clc
echo on
%*********************************************************
%
% Model predictive control example
%
%*********************************************************
%
% In this example, we solve quadratic problems using
% semidefinite programming, second order cone programming
% and quadratic programming.

% MPC settings
N = 5;
pause % Strike any key to continue. 

% Create numerics for a discretized double integrator
A = [2 -1;1 0];
B = [1;0];
C = [0.5 0.5];

[H,S] = create_CHS(A,B,C,N)
pause % Strike any key to continue. 

% Initial state
x = [2;0];

% Define free variables
t = sdpvar(1,1);
U = sdpvar(N,1);   
pause % Strike any key to continue. 

% Define the prediction vector
Y = H*x+S*U; 
pause % Strike any key to continue. 

% Control constraints 
F = set(-1 < U < 1);

% Terminal constraint
F = F+set(-1 < Y(N) < 1);  
pause % Strike any key to continue. 

% Our goal is to minimize the quadratic function Y'*Y+U'*U
%
% Performance constraint written as an SDP using Schur complement
% (very inefficient way to solve a QP...)
F = F+set([t Y' U';Y eye(N) zeros(N,N);U zeros(N,N) eye(N)]>0)
pause % Strike any key to continue. 

% Solve
sol = solvesdp(F,t);
pause % Strike any key to continue. 

% Look at solution
double(U) 
double(Y(N)) 
pause % Strike any key to continue. 
clc
% More efficient implementation using SOCP... (if SOCP solver available)
F =  set(-1 < U < 1) + set(-1 < Y(N) < 1) + set(cone([Y;U],t));
sol = solvesdp(F,t);
pause % Strike any key to continue. 
clc

% Even more efficient implementation if QP solver is available
F =  set(-1 < U < 1) + set(-1 < Y(N) < 1);
sol = solvesdp(F,Y'*Y+U'*U);
pause % Strike any key to continue. 


echo off