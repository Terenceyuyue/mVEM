%MPT_DEMO3 Explains control routines for LTI systems
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Explains control routines of the MPT toolbox for LTI systems
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
disp('This tutorial is devoted to introduce the functionality of');
disp('the MPT toolbox in terms of control algorithms.');
fprintf('\n\n');
disp('You can solve the following optimization problems with MPT:');
disp('');
disp(' - Constrained Finite Time Optimal Control Problem for constrained linear and PWA systems');
disp(' - Constrained Infinite Time Optimal Control Problem for constrained linear and PWA systems');
disp(' - Time-Optimal Control Problem for constrained linear and PWA systems');
disp(' - Low-complexity controllers for constrained linear and PWA systems');
fprintf('\n');
disp('In the linear case, the toolbox supports also additive and parametric uncertainty.');
disp('For Piecewise-affine (PWA) systems, additive disturbance can be specified.');
fprintf('\n\n');
disp('Press any key to continue...');
pause
clc

disp('In order to solve a specific problem in a multi-parametric fashion, the user needs');
disp('to specify both the SYSTEM and the PROBLEM he wants to solve.');
fprintf('\n');
disp('Let us first focus on Linear Time-Invariant (LTI) systems:');
disp('x(k+1)=Ax(k)+Bu(k)');
disp('y(k)=Cx(k)+Du(k)');
fprintf('\n');
disp('subject to constraints:');
disp('ymin <= y(k) <= ymax');
disp('umin <= u(k) <= umax');
fprintf('\n\n');
disp('To specify the dynamical system, the user has to fill out the sysStruct structure:');
fprintf('\n');
disp('>> sysStruct.A = [1 1; 0 1];');
disp('>> sysStruct.B = [1; 0.5];');
disp('>> sysStruct.C = [1 0; 0 1];');
disp('>> sysStruct.D = [0; 0];');
fprintf('\n');
disp('>> sysStruct.ymin = [-5; -5];');
disp('>> sysStruct.ymax = [5; 5];');
disp('>> sysStruct.umin = -1;');
disp('>> sysStruct.umax = 1;');
fprintf('\n\n');
disp('The description of the problem the user wants to solve is provided by means of the');
disp('problem structure probStruct:');
fprintf('\n');
disp('Prediction horizon:');
disp('>> probStruct.N = 5;');
fprintf('\n');
disp('Norm of the state in the cost function (allowed values are 1, 2, and Inf):');
disp('>> probStruct.norm = 2;');
fprintf('\n');
disp('Weight on states in the cost function:');
disp('>> probStruct.Q = eye(2);');
fprintf('\n');
disp('Weight on inputs in the cost function:');
disp('>> probStruct.R = 1;');
fprintf('\n\n');
disp('Press any key to continue...');
pause

Options.newfigure=1;
try
    evalin('base','load ltiDIfh');
catch
    error('Demos data not accessible! Please add mpt/demos to your path.');
end
fprintf('\n');
disp('Now, having both the system as well as the problem specified, you can start the computation:');
fprintf('\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
disp('>> plot(ctrl);');
fprintf('\n');
fprintf('The solution consists of %d regions and is depicted on your screen.\n', length(ctrlStruct.Pn));
mpt_plotPartition(ctrlStruct,Options);
fprintf('\n\n');
disp('Press any key to continue...');
pause
clc
try
    evalin('base','load ltiDIih');
catch
    error('Demos data not accessible! Please add mpt/demos to your path.');
end
disp('MPT_CONTROL is the main control routine of the toolbox. Depending on the values');
disp('  of sysStruct and probStruct it solves the particular problem by calling one of');
disp('  the auxiliary functions.');
fprintf('\n\n');
disp('Let us now compute an infinite-time solution for the same LTI system as before.');
fprintf('\n');
disp('Only change that is necessary is to set the prediction horizon to infinity:');
fprintf('\n');
disp('>> probStruct.N = Inf;');
fprintf('\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
disp('>> plot(ctrl);');
fprintf('\n');
fprintf('The controller, which consists of %d regions, is plotted on your screen.\n', length(ctrlStruct.Pn));
plot(ctrlStruct,Options);
fprintf('\n\n');    
disp('Press any key to continue...');
pause

clc
try
    evalin('base','load ltiDIlc');
catch
    error('Demos data not accessible! Please add mpt/demos to your path.');
end
disp('Using similar modification, the user can request the computation of low-complexity');
disp('  controllers by providing a level of sub-optimality.');
fprintf('\n');
disp('Please consult the manual for more details.');
fprintf('\n\n');
disp('>> probStruct.N = 1;');
disp('>> probStruct.subopt_lev = 2;');
fprintf('\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
disp('>> plot(ctrl);');
fprintf('\n');
fprintf('The resulting reduced-complexity explicit controller is defined over %d regions.\n',length(ctrlStruct.Pn));
plot(ctrlStruct,Options);
fprintf('\n\n');
disp('Press any key to continue...');
pause

clc
try
    evalin('base','load ltiDIadd');
catch
    error('Demos data not accessible! Please add mpt/demos to your path.');
end

disp('The MPT toolbox comes along with a set of algorithms that are able to calculate explicit');
disp('  control laws that are robust against additive and polytopic uncertainties.');
fprintf('\n\n');
disp('To specify an additive disturbance, introduce the following field to your sysStruct:');
fprintf('\n');
disp('>> probStruct.N = 3;');
disp('>> sysStruct.noise = polytope([eye(2); -eye(2)], [1;1;1;1]*0.1);');
fprintf('\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n\n');
fprintf('The multi-parametric program results in %d regions for prediction horizon N=3.\n',length(ctrlStruct.Pn));
plot(ctrlStruct,Options);
fprintf('\n\n');
disp('The end');

