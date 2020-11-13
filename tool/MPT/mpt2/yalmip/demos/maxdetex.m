clc
echo on
%*********************************************************
%
% Determinant maximization problem
%
%*********************************************************
%
% YALMIP can be used to model determinant maximization 
% problems, and solve these problems using any SDP solver 
% (i.e., the dedicated solver MAXDET is not needed)
pause
clc
% Consider the discrete-time system x(k+1)=Ax(k+1)+Bu(k), u = -Lx
A = [1 0;0.4 1];
B = [0.4;0.08]; 
L = [1.9034 1.1501];

% We want to find the largest possible invariant ellipsoid x'(Y^-1)x < 1
% for which a control constraint |u|<1 is satisfied
% This can be formulated as a determinant maximization problem
pause % Strike any key to continue. 

% Define the symmetric matrix Y
Y = sdpvar(2,2);
pause % Strike any key to continue. 

% Define LMI for invariance and constraint satisfaction
F = set( [Y Y*(A-B*L)';(A-B*L)*Y Y] > 0);
F = F + set(L*Y*L' <1);
pause % Strike any key to continue. 

% Y should be positive definite and we want to maximize det(Y)
%
% In YALMIP, we can solve this by maximizing det(Y)^(1/(length(Y)))
% which can be modeled using semidefinite and second order cone
% constraints. This function, the geometric mean of the eigenvalues,
% is available in the nonlinear operator geomean.
%
% If you want to get slightly smaller models, you may want to
% replace the operator geomean with geomean2. For matrices with
% dimension not a power of 2, the associated SDP model will
% become a bit more efficiently constructed. Note though that
% geomean2 does not model the geometric mean, but a monotonic
% transformation of this function.

% NOTE : if you want to use the dedicated solver MAXDET, you must use
% the objective function -logdet(Y) instead. The command logdet may
% however become obsolete in a future version.
%
pause

solution = solvesdp(F,-geomean(Y));
pause % Strike any key to continue. 

% Get result
Y = double(Y)
pause % Strike any key to continue.
echo off