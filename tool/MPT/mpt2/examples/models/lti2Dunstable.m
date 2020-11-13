%LTI2DSTABLE 2nd order stable LTI example
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
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
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

% x(k+1)=Ax(k)+Bu(k)
sysStruct.A = [5.049 1.17; 9.362 2.708];
sysStruct.B = [1.17; 1.708];

% y(k)=Cx(k)+Du(k)
sysStruct.C = [0 0.5];
sysStruct.D = 0;

% set constraints on inputs
sysStruct.umax = 2;
sysStruct.umin = -1;

% constraints on output
sysStruct.ymax = 30;
sysStruct.ymin = -30;

% constraints on states
sysStruct.xmax = [10; 10];
sysStruct.xmin = [-10; -10];

% use quadratic cost
probStruct.norm = 1;

% penalty on states
probStruct.Q = eye(2);

% penalty on inputs
probStruct.R = 1;

% prediction horizon
probStruct.N = 3;
