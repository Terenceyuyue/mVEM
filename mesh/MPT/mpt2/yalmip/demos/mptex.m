yalmip('clear')
clc
echo on
%*********************************************************
%
% Multi-parametric programing
%
%*********************************************************
%
% This example requires the MPT-Toolbox
%
% We return to the MPC-problem. 
% This time however, we seek the explicit solution as a 
% function of the current state x
pause

% Create numerics for a discretized double integrator
N = 5;
A = [2 -1;1 0];
B = [1;0];
C = [0.5 0.5];
[H,S] = create_CHS(A,B,C,N);
pause

% Variables
x = sdpvar(2,1);
U = sdpvar(N,1);
pause

% Parameterized predictions
Y = H*x+S*U; 

%Control constraints
F = set(-1 < U < 1);

% Terminal constraint
F = F+set(-1 < Y(N) < 1); 
pause

% Exploration space
F = F + set(-5 < x < 5);
pause

% Call solvemp (same syntax as solvesdp, except additional 
% fourth argument to specify the parametric variable)
mpsol = solvemp(F,Y'*Y+U'*U,[],x);
pause

% We can, e.g., plot the solution
plot(mpsol{1}.Pn)
pause

% Solve the problem with a L1 cost
mpsol = solvemp(F,Y'*Y + U'*U,[],x);
pause

% Plot the value function
clf
mpt_plotPWA(mpsol{1}.Pn,mpsol{1}.Bi,mpsol{1}.Ci)
pause


% Solve the problem with a L1 cost, but this time
% only request U(1) in the returned solution
mpsol = solvemp(F,norm([Y;U],1),[],x,U(1));
pause

% Plot the optimizer
clf
mpt_plotPWA(mpsol{1}.Pn,mpsol{1}.Fi,mpsol{1}.Gi)
pause

% Extract value of optimizer at point [0.1;0.2]
[ii,jj] = isinside(mpsol{1}.Pn,[0.1;0.3]);
mpsol{1}.Fi{jj}*[0.1;0.3] + mpsol{1}.Gi{jj}
pause

% To avoid learning MPT commands, simply use 
% some more outputs from SOLVEMP, and YALMIP PWA
% functions will be generated automatically
[mpsol,diagnost,Z,Valuefunction,Optimizer] = solvemp(F,norm([Y;U],1),[],x,U(1));
clf
plot(Valuefunction)
pause

clf
plot(Optimizer)
pause

assign(x,[0.1;0.3]);
double(Optimizer)
pause

clc
% The Valuefunction and Optmizer are standard
% variables in YALMIP (so called nonlinear operators)
% which we can use as variables when setting up other
% optmimization problem
%
% Here we solve the silly problem of finding the state
% with maximal x(1) coordinate, while having an optmial
% cost less than 10.
%
% Since the value function for a simple parametric LP 
% is convex and can be written as the maximum of a set 
% of linear functions, the problem will be an LP, and 
% YALMIP keeps track of convexity etc, since the value
% function behind the scenes is described with a so 
% called convexity-aware nonlinear operator.
pause
solvesdp(set(Valuefunction < 10),-x(1)),
double(x)
double(Valuefunction)
pause

% Why not solve a nonconvex problem and find the maximum 
% of the convex value function. Note that this will generate
% a binary LP problem, so you need to have an MILP solver
pause
solvesdp([],-Valuefunction)
double(x)
double(Valuefunction)
pause

% Note that the transparant use of the value functions and the
% optimizer as standard piecewise YALMIP variables still is 
% under development. Quadratic functions are for instance 
% not fully supported yet.
%
% Learn a lot more in the HTML manual...

pause
echo off
