%PWA_DI 2nd order PWA example with 4 PWA dynamics
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Sample PWA system with 2 states and 1 control input.
% The affine dynamics is defined over 4 different regions.
% 
% Note: with probStruct.subopt_lev=1 the solution consists of 220 regions
%       with probStruct.subopt_lev=2 the solution consists of 144 regions
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
%
% USAGE:
%   pwa_DI
%   ctrlStruct = mpt_control(sysStruct,probStruct);
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

%   This file contains the following data:
%
%   Q,R                     Cost objective  J=sum x'Qx + u'Ru
%   A,B,f                   x(k+1)=Ax(k)+Bu(k)+f
%   C,D,g                   y(k)=Cx(k)+Du(k)+g
%   guardX, guardU, guardC: The dynamics are active if the following holds guardX*x(k)+guardU*u(k)<=guardC
%   bndA,bndb               Defines the area of interest. Can also be used to impose constraints on x(0).
%
%   All this data is stored in the struct "sysStruct" which is subsequently stored in the cell array "pwa"
%   which contains an entry for each dynamics.

clear sysStruct probStruct
%++++++++++++++++++++++++++++++
%   PARAMETERS FOR ALL DYNAMICS
%++++++++++++++++++++++++++++++

%set constraints on output
sysStruct.ymin    =   [-5;-5];
sysStruct.ymax    =   [5;5];

%set constraints on input
sysStruct.umin    =   -1;
sysStruct.umax    =   1;

%set constraints on input slew rate
sysStruct.dumin   =   -inf;
sysStruct.dumax   =   inf;


%define the area of state space which is of interest (can be VERY large)
%this is merely a box which defines where to compute the controller
sysStruct.Pbnd=polytope([eye(2); -eye(2)], [1;1;1;1]*50);

%++++++++++++++++++++++++++
% DYNAMICS 1
%++++++++++++++++++++++++++
%x(k+1)=Ax(k)+Bu(k)+f
sysStruct.A{1}=[1 1;0 1];
sysStruct.B{1}=[1;0.5];
sysStruct.f{1}=[0;0];
%set where dynamics are active guardX*x+guardU*U<=guardC
sysStruct.guardX{1} = [-1 0;0 -1];
sysStruct.guardC{1} = [0;0];
%y(k)=Cx(k)+Du(k)+g
sysStruct.C{1}=[1 0;0 1];
sysStruct.D{1}=[0;0];
sysStruct.g{1}=[0;0];


%++++++++++++++++++++++++++
% DYNAMICS 2
%++++++++++++++++++++++++++
%x(k+1)=Ax(k)+Bu(k)+f
sysStruct.A{2} = [1 1;0 1];
sysStruct.B{2}=[-1;-0.5];
sysStruct.f{2}=[0;0];
%set where dynamics are active guardX*x+guardU*U<=guardC
sysStruct.guardX{2} = [1 0; 0 1];
sysStruct.guardC{2} = [0;0];
%y(k)=Cx(k)+Du(k)+g
sysStruct.C{2}=[1 0;0 1];
sysStruct.D{2}=[0;0];
sysStruct.g{2}=[0;0];


%++++++++++++++++++++++++++
% DYNAMICS 3
%++++++++++++++++++++++++++
%x(k+1)=Ax(k)+Bu(k)+f
sysStruct.A{3} = [1 -1;0 1];
sysStruct.B{3}=[-1;0.5];
sysStruct.f{3}=[0;0];
%set where dynamics are active guardX*x+guardU*U<=guardC
sysStruct.guardX{3} = [1 0;0 -1];
sysStruct.guardC{3} = [0;0];
%y(k)=Cx(k)+Du(k)+g
sysStruct.C{3}=[1 0;0 1];
sysStruct.D{3}=[0;0];
sysStruct.g{3}=[0;0];


%++++++++++++++++++++++++++
% DYNAMICS 4
%++++++++++++++++++++++++++
%x(k+1)=Ax(k)+Bu(k)+f
sysStruct.A{4} = [1 -1;0 1];
sysStruct.B{4}=[1;-0.5];
sysStruct.f{4}=[0;0];
%set where dynamics are active guardX*x+guardU*U<=guardC
sysStruct.guardX{4} = [-1 0;0 1];
sysStruct.guardC{4} = [0;0];
%y(k)=Cx(k)+Du(k)+g
sysStruct.C{4}=[1 0;0 1];
sysStruct.D{4}=[0;0];
sysStruct.g{4}=[0;0];


sysStruct.StateName = {'Position', 'Speed'};
sysStruct.InputName = 'Acceleration';
sysStruct.OutputName = {'Position', 'Speed'};


probStruct.Q  = 1*eye(2);  % norm of the state (nx x nx)
probStruct.R  = 1;         % norm of the input (nu x nu)
probStruct.N = Inf;
probStruct.norm = 2;
probStruct.subopt_lev=1;
