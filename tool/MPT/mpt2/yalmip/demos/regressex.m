clc
echo on

% This examle will show how to define LP and QP problems
%
% We will use YALMIP to formulate couple of regression problems 
% (linear L_1, L_2 and L_inf regression)
%
pause % Strike any key to continue.  % Strike any key to continue. 

% To begin with, we generate data according to a noisy linear regression
% 
% y(t) = x'(t)*a+e(t), t=0,0.01,...,2*pi

a = [1 2 3 4 5 6];
t = (0:0.02:2*pi)';
x = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
y = x*a'+(-4+8*rand(length(x),1));
pause % Strike any key to continue. 
plot(y)
pause % Strike any key to continue. 

% We now want to estimate the parameter a, given y, and we will 
% try 3 different approaches; L_1, L_2 and L_inf. 
%
% We start by defining the decision variable

a_hat = sdpvar(1,6);
pause % Strike any key to continue. 

% By using a_hat and the regressor x, we can define the residuals

residuals = y-x*a_hat';
pause % Strike any key to continue. 

% To solve the L_1 problem, we define a variable that will work
% as a bound on |y-x'*a|
bound = sdpvar(length(residuals),1);
pause % Strike any key to continue. 

% Express that bound is absolute value of the residuals (usinhg a double-sided constraint)
F = set(-bound < residuals < bound);
pause % Strike any key to continue. 

% Minimize sum of bound
% NOTE, This takes time if you only have linprog installed
solvesdp(F,sum(bound));
% Optimal L_1 estimate
a_L1 = double(a_hat)
pause % Strike any key to continue. 


clc

% The L_2 solution is easily obtained from a QP
% (note the nonlinear objetive function. More about
%  this in demo #11)
solvesdp([],residuals'*residuals);
a_L2 = double(a_hat)
pause % Strike any key to continue. 

% Finally, let us minimize the L_inf norm. This corresponds to
% minimizing the largest (absolute value) residual
% We introduce a scalar bound on the largest value
bound = sdpvar(1,1);
F = set(-bound < residuals < bound);
pause % Strike any key to continue. 

% Solve by minimizing the bound
solvesdp(F,bound);
a_Linf = double(a_hat)
pause % Strike any key to continue. 

plot([y x*a' x*a_L1' x*a_L2' x*a_Linf'])
legend('Measured','True','L1','L2','Linf')

echo off










