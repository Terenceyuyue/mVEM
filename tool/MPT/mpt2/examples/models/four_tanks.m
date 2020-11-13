clear model problem

% define parameters of the "four tanks" model
A1 = 28; % [cm^2] cross-section of tank 1
A2 = 32; % [cm^2] cross-section of tank 2
A3 = 28; % [cm^2] cross-section of tank 3
A4 = 32; % [cm^2] cross-section of tank 4
a1 = 0.071; % [cm^2] cross-section of outlet hole of tank 1
a2 = 0.057; % [cm^2] cross-section of outlet hole of tank 2
a3 = 0.071; % [cm^2] cross-section of outlet hole of tank 3
a4 = 0.057; % [cm^2] cross-section of outlet hole of tank 4
g = 981; % [cm/s^2]
h10 = 12.4; h20 = 12.7; h30 = 1.8; h40 = 1.4;
k1 = 3.33; k2 = 3.35; kc = 0.5; v10 = 3; v20 = 3;
g1 = 0.7; g2 = 0.6;

T1 = A1/a1 * sqrt(2*h10 / g); T2 = A2/a2 * sqrt(2*h20 / g);
T3 = A3/a3 * sqrt(2*h30 / g); T4 = A4/a4 * sqrt(2*h40 / g);

A = [-1/T1 0 A3/(A1*T3) 0; 0 -1/T2 0 A4/(A2*T4); 0 0 -1/T3 0; 0 0 0 -1/T4];
B = [g1*k1/A1 0; 0 g2*k2/A2; 0 (1-g2)*k2/A3; (1-g1)*k1/A4 0];
C = [kc 0 0 0; 0 kc 0 0];
D = zeros(2);

% sampling time
Ts = 2;

% define the model by discretizing a corresponding state-space representation
model = mpt_sys(ss(A, B, C, D), Ts);

% state constraints
model.xmax = [10; 10; 10; 10];
model.xmin = [-10; -10; -10; -10];

% input constraints
model.umax = [1; 1];
model.umin = [-1; -1];

% output constraints
model.ymax = [10*kc; 10*kc];
model.ymin = [-10*kc; -10*kc];



% prediction horizon
problem.N = 2;

% control horizon
problem.Nc = 2;

% quadratic performance index
problem.norm = 2;

% penalty on states
problem.Q = eye(4);

% penalty on outputs
problem.Qy = 100*eye(2);

% penalty on inputs
problem.R = eye(2);

% the controller must track a given reference
problem.tracking = 1;
