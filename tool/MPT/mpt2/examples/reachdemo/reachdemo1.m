% reachability computation demo
% Copyright (C) 2005 by Michal Kvasnica (kvasnica@control.ee.ethz.ch)

% first we want to compute N-Step reachable sets for a linear system in the form
% x(k+1) = A*x(k) + B*u(k)
% x(0) \in X0
% u(k) \in U0
%
% i.e. compute set of states which are reachable in N steps starting from set of
% initial conditions X0 assuming the system input u(k) belongs to some bounded
% set U0
fprintf('Reachability analysis demo - Part 1\n\n');

% first define matrices of the dynamical system
A=[-1 -4; 4 -1]; 
B=[1;1]; 
C=[1 0]; 
D=0; 
syst=ss(A,B,C,D); 
systd=c2d(syst,0.02); 

% import dynamical system to MPT system structure
sysStruct=mpt_sys(systd); 

% MPT requires at least input and output constraints to be defined, altough they
% have no effect on the reachability computation
sysStruct.ymax=10; 
sysStruct.ymin=-10; 
sysStruct.umax=1; 
sysStruct.umin=-1; 

% define set of initial condintions as a polytope object
X0=polytope([0.9 0.1; 0.9 -0.1; 1.1 0.1; 1.1 -0.1]); 

% set of admissible inputs as a polytope object
U0=unitbox(1,0.1);  % inputs should be such that |u| <= 0.1

% compute the 50-Steps reachable set
disp('Computing N-Step reachable sets');
N = 50;
R=mpt_reachSets(sysStruct, X0, U0, N);

% plot the sets along with X0
plot(X0, 'r', R, 'g');
title('Initial conditions - red, Reachable sets - green');
