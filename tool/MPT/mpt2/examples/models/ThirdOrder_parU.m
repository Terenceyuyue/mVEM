%THIRDORDER_PARU 3rd order LTI example with polytopic uncertainty
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% 3rd order LTI system with 3 states and 2 control inputs
% 
% Default values are:
%   * Prediction horizon N = 5
%   * Weight in the cost function Q = I, R = 0.1I
%   * Norm = 2
%   * No bounds on x0
%   * Level of suboptimality = 0
%   * Parametric uncertainty
%
% USAGE:
%   ThirdOrder_parU
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

%definition of parametric uncertainty
pct=0.1;

%x(k+1)=Ax(k)+Bu(k)
sysStruct.A=[0.5 0.5 -0.2;0 0.7 -0.2; -0.2 0 0.5];
sysStruct.B=[0 1;1 0; 0 1];

%y(k)=Cx(k)+Du(k)
sysStruct.C=eye(3);
sysStruct.D=zeros(3,2);
	
%set constraints on output
sysStruct.ymin    =    -5*[1 1 1]';
sysStruct.ymax    =   5*[1 1 1]';

%set constraints on input
sysStruct.umin    =   -2*ones(2,1);
sysStruct.umax    =   2*ones(2,1);
 
%set constraints on input slew rate
sysStruct.dumin   =   -inf*ones(2,1);
sysStruct.dumax   =   inf*ones(2,1);


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%define uncertainty
sysStruct.Acell{1}=sysStruct.A*(1-pct);
sysStruct.Acell{2}=sysStruct.A*(1+pct);

sysStruct.Bcell{1}=sysStruct.B;
sysStruct.Bcell{2}=sysStruct.B;


probStruct.norm=2;
probStruct.Q=eye(3);   
probStruct.R=eye(2);
probStruct.N=5;
probStruct.x0bound=0;
probStruct.subopt_level=0;