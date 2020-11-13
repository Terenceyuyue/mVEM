function [sysStructCell]=mpt_randLTISys(maxSys,nx,nu,Options)
%MPT_RANDLTISYS Generates random LTI systems
%
%   [sysStructCell]=mpt_randLTISys(maxSys,nx,nu,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% 
% Generates random LTI systems
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
%   maxSys      -   Number of systems to be generated
%   nx          -   number of states of each system; scalar
%   nu          -   number of inputs of each system; scalar
%   Options
%       .ymin / .ymax       -   Output constraints for each system (Default: -/+ 10)
%       .umin / .umax       -   Input constraints for each system (Default: -/+ 1)
%       .dumin / .dumax     -   Input constraints for each system (Default: -/+ Inf)
%       .enforceStability   -   Consider stable dynamics only, if set to 1
%                               (Default: 0)
%       .maxVal             -   Maximum absolute value of elements in dynamic matrices
%                               (Default: 2)
%       .justUnstable       -   If true, only unstable systems will be returned
%                               (Default: false)
%       
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
%   sysStructCell           -   cell array of size "maxSys" containing different sysStructs
%

% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich
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

if(nargin<3)
    error('mpt_randPWASys: insufficient number of inputs')
end
if(nargin<4)
    Options=[];
end
if(isfield(Options,'enforceStability'))
    enforceStability=Options.enforceStability;
else
    enforceStability=0; %consider stable dynamics only
end
if(isfield(Options,'maxVal'))
    maxVal=Options.maxVal;
else
    maxVal=2;  %Maximum absolute value of elements in dynamic matrices
end
if(isfield(Options,'umax'))
    umax=Options.umax;
else
    umax=1*ones(nu,1);  
end
if(isfield(Options,'umin'))
    umin=Options.umin;
else
    umin=-1*ones(nu,1);  
end
if(isfield(Options,'dumax'))
    dumax=Options.dumax;
else
    dumax=Inf*ones(nu,1);  
end
if(isfield(Options,'dumin'))
    dumin=Options.dumin;
else
    dumin=-Inf*ones(nu,1); 
end
if(isfield(Options,'ymax'))
    ymax=Options.ymax;
else
    ymax=10*ones(nx,1); 
end
if(isfield(Options,'ymin'))
    ymin=Options.ymin;
else
    ymin=-10*ones(nx,1); 
end
if ~isfield(Options, 'justUnstable')
    Options.justUnstable = 0;
end






%%%%%%%%%   SYSTEM GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rand('state',sum(100*clock));

for sys_ctr=1:maxSys
    notFinished=1;
    while(notFinished)
        sysStruct.A = (rand(nx,nx)-0.5)*maxVal;
        sysStruct.B = (rand(nx,nu)-0.5)*maxVal;
        if Options.justUnstable & all(abs(real(eig(sysStruct.A)))<1)
            notFinished = 1;
        elseif(enforceStability==1 & (rank(ctrb(sysStruct.A,sysStruct.B))~=nx | ...
                any(abs(real(eig(sysStruct.A)))>1)))
            notFinished=1;
        else
            notFinished=0;
        end
    end
    
    %y(k)=Cx(k)+Du(k)
    sysStruct.C = eye(nx);
    sysStruct.D = zeros(nx,nu);
    
    %set constraints on output
    sysStruct.ymin    =   ymin;
    sysStruct.ymax    =   ymax;
    
    %set constraints on input
    sysStruct.umin    =   umin;
    sysStruct.umax    =   umax;
    sysStruct.dumin   =   dumin;
    sysStruct.dumax   =   dumax;
    
    sysStructCell{sys_ctr}=sysStruct;
end

