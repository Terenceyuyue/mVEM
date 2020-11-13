clc
echo on

% YALMIP supports mixed integer (LP/QP/SOCP/SDP) programming.  If no
% specialized MIP solver is installed, such as GLPK or CPLEX, an internal 
% branch-and bound scheme is applied, using the installed continuous 
% solvers for solving the node problems. 
%
% Notice that this is a very rudimentary branch-and-bound solver, so
% for performance and anything but small academic examples, a dedicated 
% integer solver is required.
%
% Type help bnb to obtain inforamtion about the solver.
%
% Let us test the internal solver on MILP, MIQP, MISOCP and MISDPP

pause
yalmip('clear')
clc

% We are given a noisy Toeplitz matrix P, and the goal is to find an 
% integer Toeplitz matrix Z minimizing sum_i sum_j |Pij-Zij| 
%
% Obviously, this is a MILP problem.
%
pause % Strike any key to continue. 

% Create a noisy Toeplitz matrix
n = 5;
P = toeplitz(randn(n,1)*100)+randn(n,n)*5
pause % Strike any key to continue. 

% Define an integer Toeplitz matrix
Z = intvar(n,n,'toeplitz');
pause % Strike any key to continue. 

% Introduce a variable to deal with the absolute values
% note : P is not symmetric so constraints are automatically
% interpreted as element-wise.
t = sdpvar(n,n,'full');
F = set(-t < P-Z < t);
pause % Strike any key to continue. 

% Solve (using the internal branch-and-bound solver
% Turn on display in branching, but turn off display in 
% the local LP-solver
solvesdp(F,sum(sum(t)),sdpsettings('solver','bnb','verbose',2));
pause % Strike any key to continue. 

% Check the result
P
double(Z)
pause

clc
% As an alternative, let us solve the corrsponding MIQP
% problem where we minimize sum_i sum_j |Pij-Zij|^2
pause

% Define the residuals
e = P(:)-Z(:)
% and optimize
solvesdp([],e'*e,sdpsettings('solver','bnb','verbose',2));
pause % Strike any key to continue. 

% Check the result
P
double(Z)
pause


clc
% Finally, let us do things the hard way, and formulate
% the problem a mixed integer second order cone problem!
pause

% Define the residuals
e = P(:)-Z(:)
% and a bound
t = sdpvar(1,1);
% and optimize the bound t subject to ||e||<t
solvesdp(set(cone(e,t)),t,sdpsettings('solver','bnb','verbose',2));
pause % Strike any key to continue. 

% Check the result
P
double(Z)
pause


clc
% ...or why not a plain stupid implementation using
% mixed integer semi-definite programming!
pause

% A Shur complement gives a LMI
F = set([t e';e eye(length(e))]>0);
solvesdp(F,t,sdpsettings('solver','bnb','verbose',2));
pause % Strike any key to continue. 

% Check the result
P
double(Z)
pause
echo off


