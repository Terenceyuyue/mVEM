% verification demo

% Copyright (C) 2005 by Michal Kvasnica (kvasnica@control.ee.ethz.ch)

% task - check if system states can enter in N steps some given set Xf starting
% from a set of initial conditions X0 assuming:
%
%  1. the system input belongs to some set of initial conditions U0
%  2. the system input is driven by an explicit control law

% in this example we focus on the first case - system input can be anything
% inside of a set of admissible inputs given by a polytope U0

fprintf('Verification demo - Part 1\n\n');

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

% set of final states
Xf1=unitbox(2,0.1)+[-0.2; -0.2]; 

% number of steps
N = 50;


% start verification
disp('Performing verification...');
[canreach, Nf] = mpt_verify(sysStruct, X0, Xf1, N, U0);

if canreach
    fprintf('\nTarget set CAN be reached in %d steps\n\n', Nf);
else
    fprintf('\nTarget set CAN NOT be reached in %d steps\n\n', N);
end


% now define different final set
Xf2 = unitbox(2, 0.1);

disp('Performing verification with different target set...');
[canreach, Nf] = mpt_verify(sysStruct, X0, Xf2, N, U0);

if canreach
    fprintf('\nTarget set CAN be reached in %d steps\n\n', Nf);
else
    fprintf('\nTarget set CAN NOT be reached in %d steps\n\n', N);
end



% we can double-check the results by computing reachable sets:
fprintf('\nDouble-checking results...\n');
R = mpt_reachSets(sysStruct, X0, U0, N);

% plot the results
plot(X0, 'r', Xf1, 'g', R, 'b');
title('Xf1 (green) is reachable from X0 (red)');
figure;
plot(X0, 'r', Xf2, 'g', R, 'b');
title('Xf2 (green) is not reachable from X0 (red)');
