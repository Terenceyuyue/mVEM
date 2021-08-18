
% MPC setting
N =10;

% Create numerics for a discretized double integrator
A = [2 -1;1 0];
B = 0.1*[1;0];
C = [0.5 0.5];

[H,S] = create_CHS(A,B,C,N)

% Initial state
x = [2;0];

% Define free variables
t = sdpvar(1,1);
U = [intvar(10,1)]   

% Define the prediction vector
Y = H*x+S*U; 

% Control constraints (we add a tag for later)
control_bound = 10;
F = lmi('U<control_bound','Upper bound');
F = F+lmi('U>-control_bound','Lower bound') ;
% Terminal constraint
F = F+lmi('Y(N)>-1');  
F = F+lmi('Y(N)<1'); 
t = sdpvar(2*N,1);
%F = lmi([Y;U]<t) + lmi(t>-[Y;U]);

solvesdp(F,Y'*Y+U'*U,sdpsettings('verbose',0,'bnb.solver','sedumi'))