clear sysStruct probStruct
Double_Integrator

echo on

% It is now possible to use hard terminal state constraints defined in
% "probStruct.xN"
probStruct.xN = [0; 0];

% switch off stabilizing target set, since the terminal state constraint on it's
% own guarnatees stability
probStruct.Tconstraint = 0;
probStruct.P_N = zeros(2);

% calculate an explicit controller
ctrl = mpt_control(sysStruct, probStruct);

% compute an open-loop trajectory
x0 = [2; 0];
X = sim(ctrl, x0, struct('openloop', 1));
X

% notice that the final predicted state meets given constraint
%
% now plot the closed-loop trajectory
simplot(ctrl, x0);

echo on
