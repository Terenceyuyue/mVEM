%MPT_DEMO6 Illustrates implementation and visualization of the control law
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Illustrates implementation and visualization of the control law
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
%(C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%         grieder@control.ee.ethz.ch

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
disp('This demonstration will lead you through routines capable to implement and visualize');
disp('control laws produced by the MPT toolbox.');
fprintf('\n\n');
disp('This topic involves:');
fprintf('\n');
disp(' - Computation of open-loop and closed-loop trajectories');
disp(' - Plotting of polyhedral partitions');
disp(' - Plotting of open-loop and closed-loop trajectories');
disp(' - Visualization of the control action');
disp(' - Visualization of the value function');
fprintf('\n\n');
disp('We will present the functionality on an LTI Double Integrator example.');
disp('To speed up the explanation, the data will be loaded from a saved file.');
load ltiDIfh
fprintf('\n\n');
disp('Press any key to continue...');
pause
clc

Options.newfigure=1;
disp('Plotting of polyhedral partitions:');
fprintf('\n');
disp('when one calls mpt_control to compute the control law, the output arguments are as follows:')
fprintf('\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n');
disp('where "ctrl" contains the description of the polyhedral partition over which the PWA control law,');
disp('specified by cell arrays ctrl.Fi and ctrl.Gi, is defined. The partition can be');
disp('visualized by typing');
fprintf('\n');
disp('>> plot(ctrl);');
fprintf('\n\n');
plot(ctrlStruct,Options);
disp('Press any key to continue...');
pause

clc
disp('As mentioned before, the control law is Piece-wise affine and unique in each');
disp('segment of the polyhedral partition. When implementing this kind of control,');
disp('one first needs to determine in which polytope a given state lies in. Then,');
disp('the control action is given by:');
fprintf('\n');
disp('u(k) = Fi{n}*x(k) + Gi{n}');
fprintf('\n');
disp('where the index "n" corresponds to the "active" region, i.e. the region');
disp('which x(k) belongs to.')
fprintf('\n');
disp('Moreover, when implementing the control law in a receding horizon manner, one');
disp('needs to repeat this procedure at every time instance, i.e. identify an active');
disp('region, compute the control input, apply it to the system, measure the new');
disp('state, and start again from the beginning.');
fprintf('\n\n');
disp('The MPT toolbox provides a user-friendly implementation of the receding-horizon');
disp('strategy. This feature is handled by the following function:');
fprintf('\n');
disp('[X,U,Y,D,cost] = sim(ctrl,x0,Options)');
fprintf('\n');
disp('You need to provide the following arguments:');
disp('  ctrl    - controller object obtained by mpt_control');
disp('  x0      - initial state for which you want to compute a control action');
fprintf('\n');
disp('And one optional argument:');
disp('  Options    - user parameters (see help mpt_computeTrajectory for details)');
fprintf('\n');
disp('Outputs are:');
disp('  X    - a matrix containing the time-evolution of the states subject to control');
disp('  U    - a matrix in which the control moves with respect to time are stored');
disp('  Y    - evolution of the system outputs');
disp('  D    - evolution of the additive disturbances (if present)');
disp('  cost - value of the cost function');
fprintf('\n\n');
disp('Assume now we want to control the initial state x=[-2; 1]');
fprintf('\n');
disp('>> [X,U,Y,D,cost] = sim(ctrl,[-2;1]);');
fprintf('\n');
disp('Press any key to start the computation...');
pause
[X,U,Y,D,cost] = sim(ctrlStruct,[-2;1]);
fprintf('\n');
X,U
fprintf('\n');
disp('You can see how the state converges to zero, which testifies that the control');
disp('is stabilizing. Moreover, all system constraints were satisfied.');
fprintf('\n');
disp('Press any key to continue...');
pause

clc
disp('The previous method is useful if you intend to implement the explicit control law');
disp('to some physical device. For visualization purposes, the MPT toolbox provides');
disp('the following function which plots all trajectories:');
fprintf('\n');
disp('simplot(ctrl,x0)');
fprintf('\n');
disp('Mandatory input arguments are:');
disp('  ctrl    - contoller object obtained by mpt_control');
disp('  x0      - initial state for which you want to compute a control action');
fprintf('\n\n');
disp('To visualize the control trajectories for the given initial state x0=[-6; 3], you call:');
fprintf('\n');
disp('>> simplot(ctrl,[-2; 1]);');
fprintf('\n');
Options.newfigure=1;
simplot(ctrlStruct,[-2;1],[],Options);
disp('Press any key to continue...');
pause

clc
disp('If your system is 2-dimensional, you can use a mouse point-and-click interface to');
disp('plot system trajectory:');
fprintf('\n\n');
disp('>> simplot(ctrl);');
fprintf('\n\n');
disp('On the figure which will open, use your mouse to pick up an initial state.');
disp('For that state, closed-loop trajectory will be computed and displayed on the figure.');
disp('You can choose several points. Break the function using right mouse button');
fprintf('\n');
disp('Press any key to continue...');
pause
Options.newfigure=1;
simplot(ctrlStruct,Options);
fprintf('\n');
disp('Press any key to continue...');
pause


clc
disp('If your system is 2-dimensional, you can also visualize the control action by');
disp('calling the following function:');
fprintf('\n\n');
disp('>> plotu(ctrl);');
fprintf('\n');
disp('Input arguments have the same meaning as in previous functions.');
fprintf('\n');
disp('Rotate the newly opened figure to see it from different perspectives.');
fprintf('\n');
plotu(ctrlStruct,Options);
disp('Press any key to continue...');
pause
fprintf('\n\n');

clc
disp('Each controller object contains also value function which may be PWA or PWQ');
disp('over polyhedra, depending on norm of the cost function (probStruct.norm) for');
disp('which the controller was obtained. To visualize the value function, call:');
fprintf('\n\n');
disp('>> plotj(ctrl);');
fprintf('\n');
disp('Rotate the newly opened figure to see it from different perspectives.');
fprintf('\n');
plotj(ctrlStruct);
disp('Press any key to continue...');
pause
fprintf('\n\n');

disp('The end');
