function [Pn,Fi,Gi,dynamics,probStruct]=mpt_computePWATset(sysStruct,probStruct,originindynamics,Options)
%MPT_COMPUTEPWATSET Computes a stabilizing control invariant set (+ controllers) around the origin
%
% [Pn,Fi,Gi,dynamics,probStruct]=mpt_computePWATset(sysStruct,probStruct,originindynamics,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Function computes a stabilizing piecewise linear feedback law for the PWA system
% sysStruct, i.e., for all dynamics which contain the origin. Susequently the maximal 
% invariant subset of the resulting partition is computed. Furthermore, the probstruct is 
% updated with the cost-to-go required to obtain closed-loop stability if the set is used
% as a target set in MPC algorithms.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct           - System structure in the sysStruct format
% probStruct          - Problem structure in the probStruct format
%
% originindynamics    - Optional field: if dynamics which contain the origin are already known,
%                       the user can pass this vector of integers. We recommend to leave this 
%                       field empty. 
%
% Options.verbose     - Level of verbosity (see help mpt_init for more details)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% Pn             - polytope array which describes the control invariant set
% Fi{},Gi{}      - cell arrays containing the feedback law for each dynamic i
% dynamics       - vector which associates each region in Pn to a certain dynamic i
% probStruct     - updated probStruct; cost-to-go P was updated.
%
% see also MPT_INFSETPWA, MPT_GETSTABFEEDBACK

% Copyright is with the following author(s):
%
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

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end
if nargin<4,
    Options=[];
end
if ~isfield(Options,'verbose'),
    Options.verbose = mptOptions.verbose;
end
if(nargin<2)
    error('Not enough input arguments')
end

nPWA=length(sysStruct.A);
nx=length(sysStruct.A{1});
nu=size(sysStruct.B{1},2);


if nargin<3 | isempty(originindynamics)
        %obtain dynamics which contain origin
        originindynamics=[];
        if isfield(probStruct,'xref'),
            origin = probStruct.xref;
        else
            origin = zeros(nx,1);
        end
        ctr=0;
        for ii=1:nPWA
            if all(sysStruct.f{ii}==0) | isfield(probStruct,'xref') | isfield(probStruct,'uref')
                %.... othewise the origin would not be an equilibrium point
                if all(sysStruct.guardU{ii}==0) & max(sysStruct.guardX{ii}*origin - sysStruct.guardC{ii}) <= 0,
                    originindynamics = [originindynamics ii];
                elseif(any(sysStruct.guardU{ii}~=0))
                    tempP=polytope(sysStruct.guardU{ii},-sysStruct.guardX{ii}*origin + sysStruct.guardC{ii});
                    if(isfulldim(tempP))
                        originindynamics = [originindynamics ii];
                    end
                end
            end
        end
        
        if(isempty(originindynamics))
            error('mpt_computePWATset: No dynamic has the origin as an equilibrium point !! Aborting computation...');
        end
        if Options.verbose>=1,
            disp(['origin included in: ' num2str(originindynamics)]);
        end
end



% here we compute the initial invariant target set for all dynamics which contain the origin in their interior
if Options.verbose>=1,
    disp('Computing target set');
end
userTset=0; % true if user provided the terminal set, false otherwise
%----------------------------------------------------------------------------
%Compute invariant target set for PWA system
%----------------------------------------------------------------------------

%%%COMPUTE FEASIBLE SET FOR LQR CONTROLLERS AT TIME 0
P_CL=polytope;
Ain={};
Bin={};
for jj=1:length(originindynamics),
    ii=originindynamics(jj);
    Ain{jj}=sysStruct.A{ii};
    Bin{jj}=sysStruct.B{ii};
end
gfOptions = Options;
gfOptions.verbose = (Options.verbose>1);
[F,lP,feasible] = mpt_getStabFeedback(Ain,Bin,probStruct.Q,probStruct.R,gfOptions);
probStruct.P_N=lP;

if ~feasible | any(isnan(F{1}))
    error('No stabilizing feedback around origin found ...')
end  

len_od = length(originindynamics);
A_CL = cell(1,len_od);
f = cell(1,len_od);
Hx = cell(1,len_od);
Kx = cell(1,len_od);
for ii=1:len_od
    dyn = originindynamics(ii);
    A_CL{ii}=sysStruct.A{dyn}+sysStruct.B{dyn}*F{ii};
    f{ii}=sysStruct.f{dyn};
    %%CONSTRUCT CONSTRAINT MATRICES FOR t=0
    Hx{ii}=[F{ii};-F{ii}];                         %input constraints
    Kx{ii}=[sysStruct.umax;-sysStruct.umin];     %input constraints
    if isfield(sysStruct, 'xmax'),
        Hx{ii}=[Hx{ii}; eye(nx); -eye(nx)];
        Kx{ii}=[Kx{ii}; sysStruct.xmax; -sysStruct.xmin];
    end
    if isfield(sysStruct, 'ymax'),
        Hx{ii}=[Hx{ii}; sysStruct.C{dyn}*eye(nx); -sysStruct.C{dyn}*eye(nx)];
        Kx{ii}=[Kx{ii}; sysStruct.ymax; -sysStruct.ymin];
    end
    Hx{ii}=[Hx{ii}; sysStruct.guardX{dyn}+sysStruct.guardU{dyn}*F{ii}];
    Kx{ii}=[Kx{ii}; sysStruct.guardC{dyn}];
    P_CL=[P_CL polytope(Hx{ii},Kx{ii})];
end

%COMPUTE POSITIVE INVARIANT SUBSET OF PARTITION
[Pn,dynamics]=mpt_infsetPWA(P_CL,A_CL,f,sysStruct.noise,Options); %compute invariant target set
%Associate controllers to regions
if ~isfulldim(Pn),
    error('mpt_computePWATset: Invariant set is empty! Check your system definition.');
end
len_Pn = length(Pn);
Fi = cell(1,len_Pn);
Gi = cell(1,len_Pn);
for ii=1:len_Pn
    Fi{ii}=F{dynamics(ii)};
    Gi{ii}=zeros(nu,1);
end
dynamics=originindynamics(dynamics);


if(all(~isfulldim(Pn)))
    error('mpt_computePWATset: Invariant target set is empty')
end
