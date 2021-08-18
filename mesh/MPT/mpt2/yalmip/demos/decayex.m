clc
echo on
%*********************************************************
%
% Advanced YALMIP modeling, generalized eigenvalue problem
%
%*********************************************************
% 
% The problem we solve here is estimatation of decay-rate
% of a linear system x' = Ax. 
%
% This can be formulated as
%
% max alpha
% s.t A'P+PA < -2alphaP
%          P > I
%
pause

% Define the variables
A = [-1 2;-3 -4];
P = sdpvar(2,2);
pause

 %To find a lower bound on alpha, we solve a standard Lyapunov stability problem. 
F = set(P>eye(2))+set(A'*P+P*A < -eye(2));
pause
solvesdp(F,trace(P));
P0 = double(P)
pause

%For this particular solution, the decay-rate is
alpha_lower = -max(eig(inv(P0)*(A'*P0+P0*A)))/2
pause

% We now find an upper bound on the decay-rate by doubling alpha until
% the problem is infeasible. To find out if the problem is infeasible,
% we check the problem field in the solution structure. The
% meaning of this variable is explained in the help text for the command 
% yalmiperror. Infeasibility is detected if the value is 1. To
% reduce the amount of information written on the screen, we run the
% solver in a completely silent mode. This can be accomplished by
% using the verbose and warning options in sdpsettings.
%
% Note, on some solvers, just checking infeasibility flag is not enough. A
% more careful bisection code should check also for numerical problems that 
% could indicate infeasibility
pause

options = sdpsettings('verbose',0,'warning',0);
alpha_upper = alpha_lower*2;
F = set(P>eye(2))+set(A'*P+P*A < -2*alpha_upper*P);
sol = solvesdp(F,[],options);
while ~(sol.problem==1)
  alpha_upper = alpha_upper*2;
  F = set(P>eye(2))+set(A'*P+P*A < -2*alpha_upper*P);
  sol = solvesdp(F,trace(P),options);
end
alpha_upper
pause

% Having both an upper bound and a lower bound allows us to perform a
% bisection. (see code in decayex.m)
pause
echo off

tol = 0.01;
alpha_works = alpha_lower;
clc
disp('################################');
disp('#   Lower     Upper     Current');
disp('################################');
while (alpha_upper-alpha_lower)>tol 
  alpha_test = (alpha_upper+alpha_lower)/2;
  disp([alpha_lower alpha_upper alpha_test])
  F = set(P>eye(2))+set(A'*P+P*A < -2*alpha_test*P);
  sol = solvesdp(F,trace(P),options);
  if sol.problem==1
    alpha_upper = alpha_test;
  else
    alpha_lower = alpha_test;
    alpha_works = alpha_test;
  end
end


% A lower bound on the optimal value is thus
alpha_works
pause
echo off