%RUNMOVEBLOCKING Demonstrates MPT Move Blocking routines
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Demonstrates on the Double Integrator example how the move blocking 
% routines work.
%
% Open-loop solutions are calculated for the following three cases: 
%   - no blocking on inputs u_k or differences of inputs u_{k+1} - u_k.
%   - fix inputs u_k to be constant during the prediction horizon.
%   - fix differences of consecutive inputs u_{k+1} - u_k to be constant 
%     during the prediction horizon.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% none
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% none
%

% Copyright is with the following author(s):
%
% (C) 2004 Raphael Cagienard, Automatic Control Laboratory, ETH Zurich,
%          cagienard@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
% 
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

mptOptions.verbose = 0;

% DOUBLE_INTEGRATOR 2nd order LTI example
Double_Integrator;

% Calculate open-loop solution
Options.openloop=1;

% Initial state
x0=[1; -1];

fprintf('\n\n')
disp('---------------------------------------------------------------------------')
disp('Example Move Blocking Routines')
disp('---------------------------------------------------------------------------')
fprintf('\n')
disp('---------------------------------------------------------------------------')
disp('Double Integrator')
disp('Quadratic cost objective (norm=2)')
disp('Finite horizon solution (N = 5)')
disp(['Initial state x0 = [' num2str(x0') ']'])
disp('---------------------------------------------------------------------------')
fprintf('\n')


%no blocking
disp('---------------------------------------------------------------------------')
fprintf('\n')
disp('NO BLOCKING: calculate open-loop solution for prediction horizon N = 5.');
p=input('... (press enter to continue)');
fprintf('\n')

ctrlStruct=mpt_control(sysStruct,probStruct);
figure;
[X,U,Y,cost,trajectory]=mpt_plotTimeTrajectory(ctrlStruct,x0,[],Options);
subplot(223);
h=title(['No blocking'],'FontWeight','bold');

fprintf('\n')
disp('Open-loop solution for prediction horizon N = 5 and no blocking.')
disp('Notice that the inputs are changing at each time step during the prediction')
disp('horizon (bottom, left figure).')
p=input('... (press enter to continue)');
fprintf('\n')


%inputblocking
disp('---------------------------------------------------------------------------')
fprintf('\n')
disp('INPUTBLOCKING: fix consecutive inputs u_k to be constant during the');
disp('prediction horizon.'); 
disp('Example: u_1, u_2=u_3=u_4=u_5 --> probStruct.inputblocking=[1 4] (entries');
disp('define how many consecutive inputs are constant).');
p=input('... (press enter to continue)');
fprintf('\n')

probStruct.inputblocking=[1 4];
probStruct.deltablocking=[];
ctrlStruct=mpt_control(sysStruct,probStruct);
figure
[X,U,Y,cost,trajectory]=mpt_plotTimeTrajectory(ctrlStruct,x0,[],Options);
subplot(223);
h=title(['Inputblocking=[1 4]'],'FontWeight','bold');

fprintf('\n')
disp('Open-loop solution for N = 5 and inputblocking = [1 4].')
disp('Notice that the predicted inputs are constant for the last 4 steps during')
disp('the prediction horizon (bottom, left figure).')
p=input('... (press enter to continue)');
fprintf('\n')


%deltablocking
disp('---------------------------------------------------------------------------')
fprintf('\n')
disp('DELTABLOCKING: fix differences u_k - u_{k+1} of consecutive inputs to be ');
disp('constant during the prediction horizon.'); 
disp('Example: u_1 and u_5 independent, interpolate u_2, u_3 and u_4 such that');
disp('(u_1-u_2)=(u_2-u_3)=(u_3-u_4)=(u_4-u_5) --> probStruct.deltablocking=[1 5] ');
disp('(entries define wich inputs u_k are independent).');
p=input('... (press enter to continue)');
fprintf('\n')

probStruct.inputblocking=[];
probStruct.deltablocking=[1 5];
ctrlStruct=mpt_control(sysStruct,probStruct);
figure
[X,U,Y,cost,trajectory]=mpt_plotTimeTrajectory(ctrlStruct,x0,[],Options);
subplot(223);
h=title(['Deltablocking=[1 5]'],'FontWeight','bold');

fprintf('\n')

disp('Open-loop solution for N = 5 and deltablocking = [1 5].')
disp('Notice the differences between the inputs are constent during the.')
disp('prediction horizon (bottom, left figure).')

fprintf('\n')
disp('---------------------------------------------------------------------------')
input('end ... (press enter to finish)');
fprintf('\n\n')