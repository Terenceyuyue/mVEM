function [sysStructCell,sys_dyn]=mpt_randPWASys(maxSys,nx,nu,noDyn,Options)
%MPT_RANDPWASYS generates random PWA systems
%
%   [sysStructCell,sys_dyn]=mpt_randPWASys(maxSys,nx,nu,noDyn,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% 
% Generates random PWA systems
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
%   maxSys      -   Number of systems to be generated
%   nx          -   Number of states of each system; scalar
%   nu          -   Number of inputs of each system; scalar
%   noDyn       -   Number of dynamics of the PWA system (default=4)
%                   (Note: the larger this number, the less 'random' the partitions
%                          because of the implemented construction method)
%
%   Options
%       .xmin / .xmax       -   State constraints for each system (Default: -/+ 10)
%       .umin / .umax       -   Input constraints for each system (Default: -/+ 1)
%       .dumin / .dumax     -   Input constraints for each system (Default: -/+ Inf)
%       .enforceStability   -   Consider stable dynamics only, if set to 1
%                               (Default: 0)
%       .originInDyn        -   Origin is on the boundary of the different systems
%                               (Default: 0)
%       .maxVal             -   Maximum absolute value of elements in dynamic matrices
%                               (Default: 2)
%       
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
%   sysStructCell           -   cell array of size "maxSys" containing different sysStructs
%   sys_dyn                 -   cell array of polytopes which define the dynamic partition for
%                               each sysStruct
%

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
    noDyn=3;
end
if(nargin<5)
    Options=[];
end
if(isfield(Options,'enforceStability'))
    enforceStability=Options.enforceStability;
else
    enforceStability=0; %consider stable dynamics only
end
if(isfield(Options,'originInDyn'))
    originInDyn=Options.originInDyn;
else
    originInDyn=0;  %only 1 dynamic contains the origin in its interior
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
if(isfield(Options,'xmax'))
    xmax=Options.xmax;
else
    xmax=10*ones(nx,1); 
end
if(isfield(Options,'xmin'))
    xmin=Options.xmin;
else
    xmin=-10*ones(nx,1); 
end
if(noDyn<4)
    error('The number of dynamics must be at least 4 for random system generation')
end

%%%%%%%%%   SYSTEM GENERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noDyn=noDyn-1;                  %subtract initial set which is computed anyway
minVal=min(abs([xmin;xmax]));   %initialize

for sys_ctr=1:maxSys
    
    %%%%%%%%%  dynamic definition
    rand('state',sum(100*clock));
    
    %random polytope for dynamic number 1
    if(originInDyn)
        %only 1 dynamics contains the origin in the interior
        vert=[rand(nx+1,nx)*minVal]; %add more points for more dynamics    
        %make sure there is one vertex in each quadrant
        for i=1:nx
            tmp=zeros(1,nx);
            tmp(i)=minVal;
            vert(i,:)=vert(i,:)-tmp;
        end
        
        dyn1=hull(vert);
        
        while(nconstr(dyn1)<noDyn)
            i=mod(i,nx);     
            i=i+1;
            tmp=zeros(1,nx);
            tmp(i)=minVal;
            
            vert=[vert;rand(1,nx)*minVal-tmp];
            dyn1=hull(vert); 
        end
    else
        %several dynamics touch the origin
        vert=[rand(nx,nx)*minVal; zeros(1,nx)]; %add more points for more dynamics, but keep origin
        dyn1=hull(vert);
        while(nconstr(dyn1)<noDyn)
             vert=[vert; rand(1,nx)*minVal];
             dyn1=hull(vert);
        end
    end
    
    if(minVal==0)
        error('Random system generation does not work for constraints set to zero... sorry:)')
    end
    
    dynamics=polytope([eye(nx);-eye(nx)],[xmax;-xmin]);        %enforce constraint satisfaction for x0
    dynamics=[dynamics\dyn1 dyn1];
    sys_dyn{sys_ctr}=dynamics;
    
    for i=1:length(dynamics)
        [sysStruct.guardX{i},sysStruct.guardC{i}]=double(dynamics(i));
    end
    %%%%%%%%%   generate system dynamics
    for dyn_ctr=1:length(dynamics)
        notFinished=1;
        while(notFinished)
            sysStruct.A{dyn_ctr}= (rand(nx,nx)-0.5)*maxVal;
            sysStruct.B{dyn_ctr}= (rand(nx,nu)-0.5)*maxVal;
            if(enforceStability==1 & rank(ctrb(sysStruct.A{dyn_ctr},sysStruct.B{dyn_ctr}))~=nx)
                notFinished=1;
            else
                notFinished=0;
            end
        end
        
        %y(k)=Cx(k)+Du(k)
        sysStruct.C{dyn_ctr} = eye(nx);
        sysStruct.D{dyn_ctr} = zeros(nx,nu);
        
        %set constraints on output
        sysStruct.ymin    =   xmin;
        sysStruct.ymax    =   xmax;
        
        %set constraints on input
        sysStruct.umin    =   umin;
        sysStruct.umax    =   umax;
        sysStruct.dumin   =   dumin;
        sysStruct.dumax   =   dumax;
    end
    
    sysStructCell{sys_ctr}=sysStruct;
end

