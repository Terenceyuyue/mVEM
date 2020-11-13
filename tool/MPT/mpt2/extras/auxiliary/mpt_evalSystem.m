function [A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct)
%MPT_EVALSYSTEM Extracts data from sysStruct and probStruct structures
%
% [A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Extracts informations from the system and problem definition structures
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct      - System structure in the sysStruct format
% probStruct     - Problem structure in the probStruct format
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% A,B,C,D        - System dynamics matrices
% Q,R            - Weighting matrices in the cost index
% ymin, ymax     - min/max constraints on the output
% umin, umax     - min/max constraints on the control input
% dumin, dumax   - min/max constraints on the slew rate of control input
% bndA, bndb     - region of exploration bndA*x<=bndb
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

%
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

error(nargchk(2,2,nargin));

if ~isfield(sysStruct,'verified'),
    sysStruct=mpt_verifySysStruct(sysStruct);
end

if ~isfield(probStruct,'verified'),
    probStruct=mpt_verifyProbStruct(probStruct);
end

if iscell(sysStruct.A),
    A=sysStruct.A{1};
    B=sysStruct.B{1};
    C=sysStruct.C{1};
    D=sysStruct.D{1};
else
    A=sysStruct.A;
    B=sysStruct.B;
    C=sysStruct.C;
    D=sysStruct.D;
end
    
ymin=sysStruct.ymin;
ymax=sysStruct.ymax;
umin=sysStruct.umin;
umax=sysStruct.umax;
dumin=sysStruct.dumin;
dumax=sysStruct.dumax;

if iscell(probStruct.Q),
    Q = probStruct.Q{end};
else
    Q=probStruct.Q;
end
if iscell(probStruct.R),
    R = probStruct.R{end};
else
    R=probStruct.R;
end

if(isfield(sysStruct,'Pbnd'))
    [bndA, bndb]=double(sysStruct.Pbnd);
    Pbnd = sysStruct.Pbnd;
elseif isfield(sysStruct,'bndA'),
    bndA=sysStruct.bndA;
    bndb=sysStruct.bndb;
else
    bndA=[];
    bndb=[];
end

if(~isfield(sysStruct,'polyUncert'))
    polyUncert=0;
else
    polyUncert=sysStruct.polyUncert;
end
if(~isfield(sysStruct,'addUncert'))
    addUncert=0;
else
    addUncert=sysStruct.addUncert;
end