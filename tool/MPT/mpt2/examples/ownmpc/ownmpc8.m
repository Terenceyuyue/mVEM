Double_Integrator
close all

echo on

% Switch off stabilizing target set
probStruct.Tconstraint = 0;

% Set penalty on last predicted state
probStruct.P_N = 100*eye(2);

% Prediction horizon 1
probStruct.N = 1;

% Use linear norm
probStruct.norm = 1;

% Construct constraints and objective
[C, O, V] = mpt_ownmpc(sysStruct, probStruct);

% Colision avoidance constraints:
% we would like the controller to steer the system in a way such that the
% system states never enter a set of "unsafe" states.
%
% First we define the set of unsafe states as a polytope object
Punsafe = unitbox(2, 0.3) + [1; 1];

% Now compute the set of "safe" state using set difference 
Pbox = unitbox(2, 10);
Psafe = Pbox \ Punsafe;

% Now add constraints telling that each predicted state must lie inside of the
% "safe" set of states
for k = 1:length(V.x),
    C = C + set(ismember(V.x{k}, Psafe));
end

% Compute an explicit controller with custom constraints
ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V);

% Plot a closed-loop trajectory
x0 = [-1; 2.2];
simplot(ctrl, struct('x0', x0));

% Notice how system states avoid the "unsafe" (empty) box-shaped set 

echo off
