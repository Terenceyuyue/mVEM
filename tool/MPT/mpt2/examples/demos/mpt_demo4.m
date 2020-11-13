%MPT_DEMO4 Control of PWA systems
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Explains control routines for PWA systems
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
disp('The MPT toolbox contains routines able to compute the explicit solution for PWA systems.');
disp('Let''s focus on a sample PWA system consisting of 2 dynamics.');
fprintf('\n\n');
disp('The affine state and output update functions have the form:');
disp('  x(k+1) = A{i} x(k) + B{i} u(k) + f{i}');
disp('  y(k)   = C{i} x(k) + D{i} u(k) + g{i}');
fprintf('\n');
disp('Dynamic no. 1 is defined as follows:');
fprintf('\n');
disp('>> sysStruct.A{1} = 0.8*[cos(-pi/3) -sin(-pi/3); sin(-pi/3) cos(-pi/3)];');
disp('>> sysStruct.B{1} = [0;1];');
disp('>> sysStruct.f{1} = [0;0];');
disp('>> sysStruct.C{1} = [1 0;0 1];');
disp('>> sysStruct.D{1} = [0;0];');
disp('>> sysStruct.g{1} = [0;0];');
fprintf('\n');
disp('The dynamic is active in a polyhedral set given by the following equation:');
disp('  guardX{1}*x(k) <= guardC{1}');
fprintf('\n');
disp('>> sysStruct.guardX{1} = [1 0];');
disp('>> sysStruct.guardC{1} = [  0];');
fprintf('\n\n');    
disp('And then dynamic no. 2:');
fprintf('\n');
disp('>> sysStruct.A{2} = 0.8*[cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)];');
disp('>> sysStruct.B{2} = [0;1];');
disp('>> sysStruct.f{2} = [0;0];');
disp('>> sysStruct.C{2} = [1 0;0 1];');
disp('>> sysStruct.D{2} = [0;0];');
disp('>> sysStruct.g{2} = [0;0];');
fprintf('\n');
disp('The dynamic is active in a polyhedral set given by the following equation:');
disp('  guardX{2}*x(k) <= guardC{2}');
fprintf('\n');
disp('>> sysStruct.guardX{2} = [-1 0];');
disp('>> sysStruct.guardC{2} = [   0];');
fprintf('\n');
disp('Along with the system constraints:');
disp('>> sysStruct.ymin = [-10;-10];');
disp('>> sysStruct.ymax = [10;10];');
disp('>> sysStruct.umin = -1;');
disp('>> sysStruct.umax = 1;');
fprintf('\n\n');
disp('Press any key to continue...');
pause

clc
disp('Various optimization objectives can be applied to PWA systems:');
disp(' - Cost optimal solution for finite horizon (probStruct.subopt_lev=0)');
disp(' - Infinite time optimal solution (probStruct.subopt_lev=0)');
disp(' - Minimum time approach (probStruct.subopt_lev=1)')
disp(' - Low complexity strategy (probStruct.subopt_lev=2)');
fprintf('\n');
disp('(the first two options are available only for linear performance index, i.e. probStruct.norm has to be either 1 or Inf )');
fprintf('\n\n');
disp('>> probStruct.Q = eye(2);');
disp('>> probStruct.R = 1;');
disp('>> probStruct.norm = Inf;');
fprintf('\n');
disp('Cost optimal solution for finite prediction horizon is enforced by:');
fprintf('\n');
disp('>> probStruct.N = 5;');
disp('>> probStruct.subopt_lev = 0;');
fprintf('\n\n');
disp('The controller is then obtained by');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n');
disp('and is depicted on the figure');
try
    evalin('base','load pwa_sincos_n5');
catch
    error('Demos data not accessible! Please add mpt/demos to your path.');
end
Options.newfigure = 1;
plot(ctrlStruct,Options);
fprintf('\n\n');
disp('Press any key to continue...');
pause

clc
disp('Infinite time optimal solution will be computed if probStruct is set as follows:');
fprintf('\n');
disp('>> probStruct.N = Inf;');
disp('>> probStruct.subopt_lev = 0;');
fprintf('\n\n');
disp('The controller is calculated with');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n');
disp('and depicted on the figure');
try
    evalin('base','load pwa_sincos_it');
catch
    error('Demos data not accessible! Please add mpt/demos to your path.');
end

Options.newfigure = 1;
plot(ctrlStruct,Options);
fprintf('\n\n');
disp('Press any key to continue...');
pause

clc
disp('Minimum time controller is specified by the following combination of fields:');
fprintf('\n');
disp('>> probStruct.norm = 2; % can be any norm, i.e. 1, 2, Inf');
disp('>> probStruct.N = Inf;');
disp('>> probStruct.subopt_lev = 1;');
fprintf('\n\n');
disp('and can be obtained in the usual way:');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n');
disp('Partitions are shown on the figure.');
try
    evalin('base','load pwa_sincos_mt');
catch
    error('Demos data not accessible! Please add mpt/demos to your path.');
end
Options.newfigure = 1;
plot(ctrlStruct,Options);
fprintf('\n\n');
disp('Press any key to continue...');
pause


clc
disp('Low complexity controller is defined by:');
fprintf('\n');
disp('>> probStruct.norm = 1; % can be either 1-norm of Inf-norm');
disp('>> probStruct.N = Inf;');
disp('>> probStruct.subopt_lev = 2;');
fprintf('\n\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n');
disp('>> plot(ctrl)');
try
    evalin('base','load pwa_sincos_lc');
catch
    error('Demos data not accessible! Please add mpt/demos to your path.');
end
Options.newfigure = 1;
plot(ctrlStruct,Options);
fprintf('\n\n');
disp('Press any key to continue...');
pause


clc
disp('Similarly as in the linear case, additive disturbance can be specified for PWA systems:');
fprintf('\n\n');
disp('>> sysStruct.noise = polytope([eye(2); -eye(2)], [1;1;1;1]*0.5);');
fprintf('\n');
disp('A call to:');
fprintf('\n');
disp('>> ctrl = mpt_control(sysStruct, probStruct);');
fprintf('\n\n');
disp('will then lead to a robustly stabilizing controller for a given PWA system.');
fprintf('\n\n');
disp('The end');
fprintf('\n\n');
disp('We strongly recommend consulting the manual and help files');
disp('to find out more about control features of the MPT toolbox.');
disp('You can also try the "runExample" demo to choose different');
disp('sample systems and problems to solve.');
