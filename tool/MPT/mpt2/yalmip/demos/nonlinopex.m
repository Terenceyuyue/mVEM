clc
echo on

% A recent and devloping extension in YALMIP is support
% for nonlinear operators such as min, max, norm and more.
%
% Although nonlinear, and often non-differentiable, the resulting
% optmization problem is in many cases still convex, and
% can be modelled using suitable additional variables and constraits.
%
% If these operators are used, YALMIP will derive a suitable
% convex model and solve the resulting problem. It may also happen
% that YALMIP fails to build a convex model, since the rules to detect
% convexity only are sufficient but not necessary.
%
% These extended operators should only be used of you now how
% to model them manually, why it can be done and when it can be done.

pause % Strike any key to continue. 
yalmip('clear')
clc

% To begin with, define some scalars
sdpvar x y z

% Nonlinear expressions are easily defined
p = min(x+y+norm([x;y],1),x)-max(min(x,y),max(x,y));
pause

% The result is a linear variable, but it is special.
% This can be seen when displayed (note the "derived")
p
pause

% These expressions can be used as any other expression
% in YALMIP. The difference is when optmization problems
% are solved. YALMIP will start by trying to expand the
% definitions of the derived variables, and try to maintain
% convexity while doing so.
pause

% Let us solve the linear regression again (from DEMO 2)
a = [1 2 3 4 5 6];
t = (0:0.2:2*pi)';
x = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
y = x*a'+(-4+8*rand(length(x),1));
a_hat = sdpvar(1,6);
residuals = y-x*a_hat';
pause % Strike any key to continue. 


% Minimize L1 error (uses ABS operator)
solvesdp([],sum(abs(residuals)));
double(a_hat)
pause

% Minimize Linf error (uses both MAX and ABS)
solvesdp([],max(abs(residuals)));
double(a_hat)
pause

% Minimize L1 error even easier (uses NORM operator)
% NOTE : This is much faster than explicitely 
% introducing the absolute values.
solvesdp([],norm(residuals,1));
double(a_hat)
pause

% Minimize Linf error even easier (uses NORM operator)
solvesdp([],norm(residuals,inf));
double(a_hat)
pause

% Regularized solution!
solvesdp([],1e-2*norm(a_hat)+norm(residuals,inf));
double(a_hat)
pause

% Minimize Linf with performance constraint on L2
%
% First, get the best possible L2 cost
solvesdp([],norm(residuals));
optL2 = double(norm(residuals));

% Now optimize Linf with performance deteriation constraint on L2
F = set(norm(residuals) < optL2*1.2);
obj = norm(residuals,inf);
pause

solvesdp(set(norm(residuals) < optL2*1.2),norm(residuals,inf));

double(a_hat)
pause

% Well, you get the picture...
%
% Here is an example were the convexity check (correctly) fails
solvesdp(set(norm(residuals) < norm(residuals,1)),norm(residuals,inf))
pause

% ...and here is an example where YALMIP fails to prove convexity
% (this problem is convex)
solvesdp([],norm(max(0,residuals)))
pause

% The rules for convexity preserving operations are currently very simple, 
% and will most likely be improved in a future version.
%
% Still, rather complicated construction are possible.
sdpvar x y z
F = set(max(1,x)+max(y^2,z)<3)+set(max(1,-min(x,y))<5)+set(norm([x;y],2)<z);
sol = solvesdp(F,max(x,z)-min(y,z)-z);
pause

% The nonlinear operators currently supported are
%
% ABS      : Absolute value of matrix    
% MIN      : Minimum of column values    
% MAX      : Minimum of column values    
% SUMK     : Sum of k largest (eigen-) values
% SUMABSK  : Sum of k largest (by magniture) (eigen-) values 
% GEOMEAN2 : (Almost) Geometric mean of (eigen-) values (used for determinant maximization). 
%
% Adding new operators is rather straightforward and is 
% described in the HTML-manual
pause
echo off

