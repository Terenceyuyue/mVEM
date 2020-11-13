clc
echo on
%*********************************************************
%
% Geomtric programming
%
%*********************************************************
%
% Geometric programming is optimization involving nonlinear
% terms with non-integer and negative powers.
%
% MOSEK can be used to solve a special class of this, so
% called posynomial geometric programming where all coefficients
% in the objective function and the constraints are positive, and
% the decision variables are constrained to be non-negative)
pause

% To define these problems, we first note that we can define
% variables with negative and non-integer powers
%
% (note, only scalar unit sdpvar variables can be raised to negative
%  and non-integer powers. Hence (1+x)^pi is not a valid command)
x = sdpvar(1,1);
degree(x^pi)
pause


% The following example is taken from the MOSEK manual
% (note, non-negativity does not need to be specified)
t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);

obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);

F = set((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1 < 1);
pause

% Standard call to solve problem
solvesdp(F,obj)

pause
echo off
