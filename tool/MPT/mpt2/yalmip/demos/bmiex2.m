echo on
clc
% This example shows how to (locally) solve a small BMI using 
% a simple linearization-based algorithm.
%
% The main motivation of this example is not to describe a 
% particularily efficient solver, but to describe how easily
% a BMI solver can be implemented using high-level YALMIP code.
%
% The problem is to find a feedback u = Lx so that the 
% L2 gain from w to y is minimized, for the system
%
% x' = Ax+Bu+Gw
% y  = Cx
%
% This can be formulated as the BMI
%
% min t
% s.t                                 P > 0
%     [(A+BL)P+P(A+BL)+C'C PG; G'P -tI] < 0
%
pause
clc

% Create system data
A = [-1    -1    -1;
     1     0     0;
     0     1     0];
B = [1;0;0];
C = [0 0 1];
G = [-1;-1;-1];

% Define decision variables
P = sdpvar(3,3);
L = sdpvar(1,3);
t = sdpvar(1,1);
pause

% A reasonble initial guess is valuable
[L0,P0]=lqr(A,B,eye(3),1);
setsdpvar(P,P0);
setsdpvar(L,-L0);
setsdpvar(t,100);
pause
clc

% Recover all involved decision variables
% in one single variable (simplifies code later)
x = recover(getvariables([P(:);L(:);t]))
pause

% Define the nonlinear matrix (simplifies the code)
H = -[(A+B*L)'*P+P*(A+B*L)+C'*C P*G;G'*P -t];
pause

% Save old iterates and old objective function
x0 = double(x);
t0 = double(t);
pause

% Linearized constraints
F = set(linearize(H)>0) + set(P>0)    
pause

% add a trust region
F = F + set(cone(x-x0,0.5*norm(x0)));
pause

% Solve linearized problem
solvesdp(F,t)
pause

% Optimal decision variables for linearized problem
xnew = double(x);
pause

% Original problem is not guaranteed to be feasible
% Line-search for feasible (and improving) solution
pause
alpha = 1;
while (min(eig(double(H)))<0) | (min(eig(double(P)))<0) | (double(t)>t0*0.9999)
    alpha = alpha*0.5;
    setsdpvar(x,x0+alpha*(xnew-x0));
end

% Current (squared) gain
double(t)
pause

% repeat....
%
% for i = 1:15
%     
%     % Save old iterates and old objective function
%     x0 = double(x);
%     t0 = double(t);
% 
%     % Linearized constraints
%     F = set(linearize(H)>0) + set(P>0);    
%     % add a trust region
%     F = F + set(cone(x-x0,0.25*norm(x0)));
%         
%     % Solve linearized problem
%     solvesdp(F,t,sdpsettings('verbose',0));
%     
%     % Optimal decision variables for linearized problem
%     xnew = double(x);
%     
%     % Original problem is not guaranteed to be feasible
%     % Line-search for feasible (and improving) solution
%     alpha = 1;
%     while (min(eig(double(H)))<0) | (double(t)>t0*0.9999)
%         alpha = alpha*0.5;
%         setsdpvar(x,x0+alpha*(xnew-x0));
%     end
%     double(t)
% end
pause

echo off
for i = 1:15
    % Save old iterates and old objective function
    x0 = double(x);
    t0 = double(t);

    % Linearized constraints
    F = set(linearize(H)>0) + set(P>0);
    % add a trust region
    F = F + set(cone(x-x0,0.5*norm(x0)));
        
    % Solve linearized problem
    solvesdp(F,t,sdpsettings('verbose',0));
    
    % Optimal decision variables for linearized problem
    xnew = double(x);
    
    % Original problem is not guaranteed to be feasible
    % Line-search for feasible (and improving) solution
    alpha = 1;
while (min(eig(double(H)))<0) | (min(eig(double(P)))<0) | (double(t)>t0*0.9999)
        alpha = alpha*0.5;
        setsdpvar(x,x0+alpha*(xnew-x0));
    end
    disp(['#' num2str(i)   '  L2 gain : ' num2str(double(t))])
end
clc
echo on
% An alternativ is to work with the complete set of 
% constraints by linearizing the LMI object
%
% Define the constraints

F = set(H>0) + set(P>0)
pause
% Solve linearized problem
solvesdp(linearize(F),t)
while (min(eig(double(H)))<0) | (double(t)>t0*0.9999)
    alpha = alpha*0.5;
    setsdpvar(x,x0+alpha*(xnew-x0));
end
double(t)

pause
% The rest is done similarily...
pause
clc

% In fact, all this can be done automatically
% by using the solver bmilin (coded using high-level YALMIP)
pause
solvesdp(F,t,sdpsettings('solver','bmilin','usex0',1));
pause
echo off


