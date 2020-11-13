%MPT_DEMO5 Tracking functionality
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Explains tracking functionality
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
%(C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch

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

clc
disp('This tutorial introduces how to design controllers which guarantee tracking properties.');
fprintf('\n\n');
disp('There are two types of tracking supported in MPT:');
disp('');
disp(' - Free state tracking where the reference point can be specified after the controller is computed');
disp(' - Fixed reference tracking where the reference point has to be specified before computation');
fprintf('\n');
disp('Free state tracking is supported only for LTI systems, whereas the fixed tracking works also for PWA systems.');
fprintf('\n\n');
disp('Press any key to continue...');
pause
clc

disp('To specify free state tracking, set the following field in your problem definition structure:');
fprintf('\n');
disp('probStruct.tracking = 1;');
fprintf('\n');
disp('The controller is then obtained by:');
fprintf('\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n\n');
disp('You can visualize closed-loop trajectories either with:');
disp('>> simplot(ctrl)');
fprintf('\n');
disp('or by:');
disp('>> simplot(ctrl, x0)');
fprintf('\n');
disp('In the first case you will be prompted to enter the reference point for the simulation.');
disp('Enter the reference point when asked, e.g. [-3; 0]');
fprintf('\n');
pause
load ltiDIft
simplot(ctrlStruct);
fprintf('\n\n');
disp('You can also provide the reference point directly:');
disp('>> Options.reference = [-4; 0];');
disp('>> simplot(ctrl, [2;2], [], Options);');
Options.reference=[-4;0];
simplot(ctrlStruct, [2;2], [], Options);
fprintf('\n\n');
disp('The resulting trajectories are depicted on the figure');
fprintf('\n\n');
disp('Press any key to continue...');
pause

% clc
% close all
% disp('An alternative to free state tracking is output regulation. In this case, a reference point')
% disp('for the output is provided in probStruct.yref:');
% fprintf('\n');
% disp('probStruct.tracking = 0;');
% disp('probStruct.yref = 2;');
% disp('probStruct.Qy = 100;');
% fprintf('\n');
% disp('Notice that an additional parameter - ''Qy'' has to be given as well for output tracking.');


clc
close all
disp('For fixed reference tracking, reference point for outputs has to be given');
disp('in the problem structure:');
fprintf('\n');
disp('probStruct.tracking = 0;');
disp('probStruct.yref = [2; 0];');
fprintf('\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n\n');
disp('Closed loop trajectories can now be plotted either by the point-and-click interface...');
disp('>> simplot(ctrl)');
fprintf('\n');
disp('or by:');
disp('>> simplot(ctrl, x0)');
fprintf('\n');
disp('In both cases, the output (in this case a one-to-one matching of the states,');
disp('is steered towards the reference point from any feasible location.');
load ltiDIrt
simplot(ctrlStruct);
simplot(ctrlStruct, [0;0], []);
fprintf('\n\n');
disp('Resulting trajectories are depicted on the figure');
fprintf('\n\n');
disp('Press any key to continue...');
pause

fprintf('\n\nThe end.');
