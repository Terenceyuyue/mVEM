% reachability computation demo

% Copyright (C) 2005 by Michal Kvasnica (kvasnica@control.ee.ethz.ch)

% we want to compute N-Step reachable sets for a linear system in the form
% x(k+1) = A*x(k) + B*u(k)
%
% assuming that
%   x(0) \in X0
%   u(k) is driven by an explicit controller
%

fprintf('Reachability analysis demo - Part 2\n\n');


% load dynamical system
Double_Integrator

probStruct.N = 3;
probStruct.norm = 1; % lower dimensional mappings happen often in 1/Inf norm problems

% compute explicit controller
expc = mpt_control(sysStruct, probStruct);


% define set of initial conditions
X0 = unitbox(2,1) + [3;0];


%=========================================================
% compute the 5-Steps reachable set
disp('Computing N-Step reachable sets');
N = 5;
R=mpt_reachSets(expc, X0, N);
%=========================================================


% plot the sets along with X0
plot(expc.Pn, 'y', X0, 'r', R, 'g');
title('Controller partition - yellow, Initial conditions (X0) - red, Reach sets - green');




fprintf('\nNow computing also lower-dimensional sets...\n\n');

% there is one very important point that has to be mentioned:
%
% by default, mpt_reachSets function returns only fully dimensional reachable
% sets. it could happen that certain mappings are lower dimensional, in which
% case you have to tell the function explicitly that you want to compute these
% lower dimensional sets as well.


% compute the 5-Steps reachable set
disp('Computing N-Step reachable sets (include lower-dimensional sets)');
N = 5;


%=========================================================
% compute lower-dimensional sets as well
Options.lowdim = 1;

[R, V]=mpt_reachSets(expc, X0, N, Options);
%=========================================================



% plot the sets along with X0
figure;
plot(expc.Pn, 'y', X0, 'r', R, 'g');
title(sprintf('Controller - yellow, X0 - red, Reach sets - green, Lower-dimensional sets - blue points'));
hold on;

% plot the lower-dimensional sets which are returned as set of vertices
for ii = 1:length(V),
    v = V{ii};
    plot(v(:,1), v(:,2), 'bx');
end