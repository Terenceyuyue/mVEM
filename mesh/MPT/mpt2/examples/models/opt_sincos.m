%PWA_SINCOS 2nd order PWA example with 2 PWA dynamics
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Example taken from:
%
%  author = {A. Bemporad and M. Morari},
%  month = March,
%  year = 1999,
%  title = {Control of Systems Integrating Logic, Dynamics, and Constraints},
%  journal = Automatica,
%  volume = 35,
%  number = 3,
%  pages = {407--427},
%  annote = {\protect{\rm Special issue on hybrid systems}},
%
%
% Default values are:
%   * Prediction horizon N = Inf
%   * Weight in the cost function Q = I, R = 1
%   * Norm = 2
%   * No bounds on x0
%   * Level of suboptimality = 1   (needed for PWA systems!)
%   * No uncertainty
%
% Note: with probStruct.subopt_lev=1 the solution consists of 174 regions
%       with probStruct.subopt_lev=2 the solution consists of 174 regions
%
%
% USAGE:
%   pwa_sincos
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
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

alpha= -pi/3;

% PWA description 1: x(k+1)=A{1} x(k) + B{1} u(k) + f{1}
sysStruct.A{1} = 0.8*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
sysStruct.B{1}=[0;1];
sysStruct.f{1}=[0;0];
% y(k) = C{1} x(k) + D{1} u(k) + g{1}
sysStruct.C{1}=[1 0;0 1];
sysStruct.D{1}=[0;0];
% guardX{1}*x(k)+guardU{1}*u(k)<=guardC{1}
sysStruct.guardX{1}=[1 0];
sysStruct.guardC{1}=[  0];

alpha= pi/3;
% PWA description 1: x(k+1)=A{2} x(k) + B{2} u(k) + f{2}
sysStruct.A{2} = 0.8*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];
sysStruct.B{2}=[0;1];

sysStruct.f{2}=[0;0];
% y(k) = C{2} x(k) + D{2} u(k) + g{2}
sysStruct.C{2}=[1 0;0 1];
sysStruct.D{2}=[0;0];
% guardX{2}*x(k)+guardU{2}*u(k)<=guardC{2}
sysStruct.guardX{2}=[-1 0];
sysStruct.guardC{2}=[  0];

%set constraints on output
sysStruct.ymin    =   [-10;-10];
sysStruct.ymax    =   [10;10];
%set constraints on input
sysStruct.umin    =   [-1];
sysStruct.umax    =   [1];
%set constraints on input slew rate
sysStruct.dumin   =   -inf;
sysStruct.dumax   =   inf;
%state constraints
sysStruct.xmax    = [20; 20];
sysStruct.xmin    = [-20; -20];


% problem definition:
probStruct.Q  = 1*eye(2);  % weight on the state (nx x nx)
probStruct.R  = 1;         % weight on the input (nu x nu)
probStruct.N = Inf;        % prediction horizon
probStruct.norm = Inf;     % infinity norm in cost function
probStruct.subopt_lev=0;   % cost-optimal solution
probStruct.y0bounds=1;     % impose constraints also on x(0)
probStruct.P_N = zeros(2); % cost-to-go