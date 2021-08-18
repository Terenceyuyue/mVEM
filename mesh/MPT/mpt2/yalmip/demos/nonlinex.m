clc
echo on

% Starting from YALMIP version 3, polynomial expressions are supported.
%
% These nonlinear expressions can be used for, e.g.
% SDPs with BMI constraints, or to solve sum-of-squares problems
%
% The nonlinear expressions are built up from sdpvar objects,
% and are manipulated in same way
pause % Strike any key to continue. 
yalmip('clear')
clc

% To begin with, define a scalar variable
x = sdpvar(1,1);
pause % Strike any key to continue. 

% Nonlinear expressions are easily built
p = 1+x+x^2+x^3;
pause % Strike any key to continue. 

% Matrices can also be nonlinear
Y = sdpvar(3,3);
Z = Y*Y+Y.*Y;
pause

% Polynomials are displayed without any symbolic information...
y = sdpvar(1,1);
p = x*x+y^4+x*y
pause

% But can be displayed better with the command sdisplay
% (this only works if the involved variables are explictely
% defined as scalars)
sdisplay(p)
pause

% Polynomials can be diffrentiated
dp = jacobian(p);
sdisplay(dp)
pause

% ...w.r.t a specific variable
dp = jacobian(p,y);
sdisplay(dp)
pause


% ...why not twice
sdisplay(jacobian(jacobian(p)'))
pause

% ...or
sdisplay(hessian(p))
pause

% Of course, all standard linear operators applies
% to the nonlinear objects
x = sdpvar(3,1);
p = trace(x*x') + sum(x.^2)
pause

clc
% Finally, a word of caution. Use yalmip('clear') 
% while working with polynomial expressions.
% 
% The reason is that every time a nonlinear variable is
% defined, a description on how it is created is
% saved inside YALMIP. With many nonlinear terms
% this list grows fast, making YALMIP slower and slower
% since the list has to be searched in when polynomial
% expressions are manipulated.
%
% (after this short session, there are already 27 nonlinear terms)
pause
echo off