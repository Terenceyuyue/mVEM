% verification demo

% Copyright (C) 2005 by Michal Kvasnica (kvasnica@control.ee.ethz.ch)

% task - check if system states can enter in N steps some given set Xf starting
% from a set of initial conditions X0 assuming:
%
%  1. the system input belongs to some set of initial conditions U0
%  2. the system input is driven by an explicit control law

% here we show how to handle the second case - system input is driven by an
% explicit control law

fprintf('Verification demo - Part 2\n\n');

% load dynamical system
Double_Integrator

% compute explicit controller
expc = mpt_control(sysStruct, probStruct);

% define set of initial condintions as a polytope object
X0 = unitbox(2,1) + [3;0];

% set of final states, does not contain the origin
Xf1=unitbox(2,0.1)+[-0.2; -0.2];



% number of steps
N = 10;

% perform verification
disp('Performing verification...');
[canreach, Nf] = mpt_verify(expc, X0, Xf1, N);

if canreach
    fprintf('\nController drives system states into Xf in %d steps\n\n', Nf);
else
    fprintf('\nController DOES NOT drive system states into Xf in %d steps\n\n', N);
end


% now define different final set, this time it includes the origin
Xf2 = unitbox(2, 0.1);

% perform verification
fprintf('\nPerforming verification with different target set...\n');
[canreach, Nf] = mpt_verify(expc, X0, Xf2, N);

if canreach
    fprintf('\nController drives system states into Xf in %d steps\n\n', Nf);
else
    fprintf('\nController DOES NOT drive system states into Xf in %d steps\n\n', N);
end


% we can check the result by computing reachable sets:
fprintf('\nDouble-check using reach set computation...\n\n');
R = mpt_reachSets(expc, X0, N);

% plot controller partition, set of initial states and set of final states
plot(expc.Pn, 'y', X0, 'r', Xf1, 'g', R, 'b');
title('Xf is not reachable from X0');


% plot controller partition, set of initial states and set of final states
figure;
plot(expc.Pn, 'y', X0, 'r', Xf2, 'g', R, 'b');
title('Xf is reachable from X0');
