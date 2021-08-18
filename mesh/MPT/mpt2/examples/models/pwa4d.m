%PWA3D 3rd order PWA example with 2 PWA dynamics
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% The example is taken from:
% 
% author = {D.Q. Mayne and S. Rakovi\'{c}},
% year = 2003,
% title = {Model predicitve control of constrained piecewise affine discrete-time systems},
% journal = {Int. J. of Robust and Nonlinear Control},
% volume = 13,
% number=3,
% month=apr,
% pages = {261--279}
%
% 
% Default values are:
%   * Norm = 2
%   * Prediction horizon N = Inf
%   * Weight in the cost function Q = I, R = 0.1
%   * No bounds on x0
%   * Level of suboptimality = 1   (needed for PWA systems!)
%   * No uncertainty 
%
% Note: with probStruct.subopt_lev=1 the solution consists of 549 regions
%       with probStruct.subopt_lev=2 the solution consists of 298 regions
%
%
% USAGE:
%   pwa3d
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


clear sysStruct probStruct
H=[eye(4);-eye(4)];
K=[10;5;10;10;10;5;10;10];

% PWA description 1: x(k+1)=Adyn{1} x(k) + Bdyn{1} u(k) + fdyn{1}
sysStruct.A{1} = [1 0.5 0.3 0.5; 0 1 1 1;0 0 1 1;0 0 0 1];
sysStruct.B{1} = [0; 0; 0; 1];
sysStruct.f{1} = [0; 0; 0; 0];
sysStruct.C{1} = [1, 0, 0, 0];
sysStruct.D{1} = [0];
sysStruct.g{1} = [0];

% guardX*x+guardU*u<=guardC
sysStruct.guardX{1}=[0 1 0 0; H];
sysStruct.guardC{1}=[  1; K];
sysStruct.guardU{1}=[  0; zeros(size(H,1),size(sysStruct.B{1},2))];

sysStruct.A{2} = [1 0.2 0.3 0.5; 0 0.5 1 1; 0 0 1 1; 0 0 0 1];
sysStruct.B{2} = [0; 0; 0; 1];
sysStruct.f{2} = [0.3; 0.5; 0; 0];
sysStruct.C{2} = [1, 0, 0, 0];
sysStruct.D{2} = [0];
sysStruct.g{2} = [0];

% guardX*x+guardU*u<=guardC
sysStruct.guardX{2}=[0 -1 0 0; H];
sysStruct.guardC{2}=[  -1; K];
sysStruct.guardU{2}=[  0; zeros(size(H,1),size(sysStruct.B{2},2))];

sysStruct.ymin=-10;
sysStruct.ymax=10;
sysStruct.umin = -1;
sysStruct.umax = 1;
sysStruct.dumin=-inf;
sysStruct.dumax=inf;

probStruct.Q  = eye(4);  % norm of the state (nx x nx)
probStruct.R  = 0.1;         % norm of the input (nu x nu)
probStruct.subopt_lev = 1;
probStruct.norm = 1;
probStruct.N=Inf;
probStruct.y0bounds=1;