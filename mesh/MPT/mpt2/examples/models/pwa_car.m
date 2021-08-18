%PWA_CAR 2nd order PWA model of a car moving on road with different slopes
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Model of a friction-less car moving in a horizontal plane with different
% slopes. The controller is required to steer the car to the origin (=zero
% position, zero acceleration) as fast as possible. Dynamical behavior of the
% system is defined by Newton's laws as follows:
%
%    dp/dt = v
%  m dv/dt = u - m g sin(alpha)
%
% discretization of the above system with sampling time of 0.1 seconds yelds the
% following affine system:
%
% [ p(k+1) ]   [1 0.1] [p(k)]   [0.005]        [c              ]
% [ v(k+1) ] = [0  1 ] [v(k)] + [0.1  ] u(k) + [-g sin(alpha_i)]
%
% The track has 4 segments (x1 denotes the horizontal position):
%
% x1 >= -0.1          ->    alpha = 0 degrees
% -4 <= x1 <= -3     ->    alpha = 0 degrees
% -3 <= x1 <= -0.1    ->    alpha = 10 degrees
% x1 <= -4           ->    alpha = -5 degrees
%
% The car is assumed to start from the 2nd region. Because of the constraints on
% control input, it is not possible to climb up the hill towards the origin
% directly, additional speed has to be gained by first moving the opposite
% direction and even then fully accelerate to reach the origin.
%
% USAGE:
%   pwa_car
%   ctrlStruct = mpt_control(sysStruct,probStruct)'
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% none
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sysStruct, probStruct - system and problem definition structures stores
%                         in the workspace
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

clear sysStruct probStruct

A = [1 0.1; 0 1];
B = [0.005; 0.1];
C = eye(2);
D = [0;0];
slope1 = pi*10/180;
slope2 = -pi*5/180;

slope3 = 0;
g = 9.81;


% PWA description 1: x(k+1)=A{1} x(k) + B{1} u(k) + f{1}
sysStruct.A{1} = A;
sysStruct.B{1} = B;
sysStruct.f{1} = [0; 0];
% y(k) = C{1} x(k) + D{1} u(k) + g{1}
sysStruct.C{1}=[1 0;0 1];
sysStruct.D{1}=[0;0];
% guardX{1}*x(k)+guardU{1}*u(k)<=guardC{1}
sysStruct.guardX{1}=[-1 0];
sysStruct.guardC{1}=[0.1];

% PWA description 1: x(k+1)=A{1} x(k) + B{1} u(k) + f{1}
sysStruct.A{2} = A;
sysStruct.B{2} = B;
sysStruct.f{2} = [0; 0];
% y(k) = C{1} x(k) + D{1} u(k) + g{1}
sysStruct.C{2}=[1 0;0 1];
sysStruct.D{2}=[0;0];
% guardX{1}*x(k)+guardU{1}*u(k)<=guardC{1}
sysStruct.guardX{2}=[1 0; -1 0];
sysStruct.guardC{2}=[-3; 4];

% PWA description 1: x(k+1)=A{1} x(k) + B{1} u(k) + f{1}
sysStruct.A{3} = A;
sysStruct.B{3} = B;
sysStruct.f{3} = [0.01541; -g*sin(slope2)];
% y(k) = C{1} x(k) + D{1} u(k) + g{1}
sysStruct.C{3}=[1 0;0 1];
sysStruct.D{3}=[0;0];
% guardX{1}*x(k)+guardU{1}*u(k)<=guardC{1}
sysStruct.guardX{3}=[1 0];
sysStruct.guardC{3}=[-4];

% PWA description 1: x(k+1)=A{1} x(k) + B{1} u(k) + f{1}
sysStruct.A{4} = A;
sysStruct.B{4} = B;
sysStruct.f{4} = [-0.02568; -g*sin(slope1)];
% y(k) = C{1} x(k) + D{1} u(k) + g{1}
sysStruct.C{4}=[1 0;0 1];
sysStruct.D{4}=[0;0];
% guardX{1}*x(k)+guardU{1}*u(k)<=guardC{1}
sysStruct.guardX{4}=[1 0; -1 0];
sysStruct.guardC{4}=[-0.1; 3];




%set constraints on output
sysStruct.ymin    =   [-7;-40];
sysStruct.ymax    =   [1;40];
%set constraints on input
sysStruct.umin    =   -5;
sysStruct.umax    =   5;
%set constraints on input slew rate
sysStruct.dumin   =   -inf;
sysStruct.dumax   =   inf;


% string labels for states, inputs and outputs
sysStruct.xlabels = { 'Position', 'Speed' }; 
sysStruct.ulabels = 'Acceleration';
sysStruct.ylabels = { 'Position', 'Speed' };


% problem definition:
probStruct.Q  = diag([100 100]);  % weight on the state (nx x nx)
probStruct.R  = 0.1;         % weight on the input (nu x nu)
probStruct.N = Inf;        % prediction horizon
probStruct.norm = 2;
probStruct.y0bounds = 1;
probStruct.subopt_lev=1;
sysStruct.Pbnd = polytope([-10 -40; -10 40; 1 40; 1 -40]);