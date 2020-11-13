function ctrlStruct=mpt_optControl(sysStruct,probStruct,Options)
%MPT_OPTCONTROL Solves the CFTOC problem for a given LTI system
%
% ctrlStruct = mpt_optControl(sysStruct,probStruct)
% ctrlStruct = mpt_optControl(sysStruct,probStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Solves a finite horizon optimal control problem for a given LTI problem:
%       x(k+1)=Ax(k)+Bu(k)
%       y(k)=Cx(k)+Du(k)
%
%   With objective function:
%    min_u ||P_N x_n||_p + \sum_{i=}^horizon (||Qx||_p + ||Ru||_p)
%   s.t.
%        ymin<= y <=ymax, umin<=u<=umax, dumin<=u(t)-u(t-1)<=dumax
%        bndA*x<=bndb
%  
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct            - System structure in the sysStruct format
% probStruct           - Problem structure in the probStruct format
% Options.verbose      - level of verbosity
% Options.lpsolver     - which LP solver to use (help mpt_solveLP)
% Options.qpsolver     - which QP solver to use (help mpt_solveQP)
% Options.step_size    - length of step over a facet
% Options.debug_level  
%         Due to numerical problems tiny regions are sometimes difficult to
%         calculate, i.e. are not identified at all. This may create "gaps"
%         in the computed control law. For the exploration, these will be
%         jumped over and the exploration in the state space will continue.
%         "debug_level" can have three values:
%      
%         0: No debug done
%         1: A tolerance is given to find gap in the region partition,
%            small empty regions inside the region partition will be discarded.
%            Note that this is generally not a problem, since the feedback law 
%            is continuous and can therefore be interpolated easily.
%            Correction to the calculation of the outer hull.
%         2: Zero tolerance to find gap in the region partition, empty regions
%            if they exist, will be detected, i.e. the user will be notified.
%            Correction to the calculation of the outer hull.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrlStruct    - Controller structure with following fields:
%   Pn,Fi,Gi    - for region Pn(i).H*x <= Pn.K(i) computed input is U=Fi{i}*x+Gi{i}   
%   Ai,Bi,Ci    - cost associated to each region (x'Aix + x'Bi + Ci)
%   Pfinal      - Defines the feasible state space partition (i.e. union of
%                 all regions) as Phard.H*x<=Phard.K
%   dynamics    - dynamics active in region Pn(i)
%   details     - contains additional information:
%      activeConstraints  - Cell Array which stores the active constraints 
%                           of the optimizer in each region.
%
% see also MPT_CONTROL, MPT_OPTINFCONTROL
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
% (C) 2003 Kari Unneland, Automatic Control Laboratory, ETH Zurich

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

error(nargchk(2,3,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end


if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end

if probStruct.subopt_lev > 0,
    error('mpt_optControl: This is an OPTIMAL control function, probStruct.subopt_lev must be 0 !');
end
if isinf(probStruct.N),
    error('Prediction horizon must be finite for this function.');
end

horizon=probStruct.N;

if(nargin<3)
    Options=[];
end
if(~isfield(Options,'lpsolver'))
    Options.lpsolver=mptOptions.lpsolver;
end
if(~isfield(Options,'setHorizon'))
    Options.setHorizon = horizon; %add set constraint as terminal set constraint
end
if(~isfield(Options,'debug_level'))
    Options.debug_level=mptOptions.debug_level;
end
if(~isfield(Options,'abs_tol'))
    Options.abs_tol=mptOptions.abs_tol;
end
if(~isfield(Options,'ispwa'))
    % even though this function does not provide a full support for PWA systems,
    % it can be used to compute iterative solution for one fiex dynamics. Since
    % sysStruct will be in PWA format, we need this flag to suppress error
    % messages and to enforce slightly different convergence checks
    Options.ispwa=0;
end

if Options.ispwa==0,
    if ~(sysStruct.type==0 | strcmp(sysStruct.type,'LTI'))
        error('mpt_optControl can only handle LTI systems! Please use mpt_optPWAcontrol for PWA systems');
    end
end

starttime = cputime;

lpsolver=Options.lpsolver;       
useSymmetry=probStruct.useSymmetry;
setHorizon=Options.setHorizon; 

if(probStruct.useSymmetry)
    [constraintsAreSymmetric]=sub1_checkSymmetry(sysStruct,probStruct);
    useSymmetry=constraintsAreSymmetric;
    if(useSymmetry==1 & constraintsAreSymmetric==0)
        disp('##############################################################################');
        disp(' mpt_mpqpControl: Constraints are not symmetric')
        disp('                  It is not possible to use symmetric features')
        disp('                  FOR MORE INFORMATION TYPE: "help mpt_mpqpControl"')
        disp('##############################################################################');
        useSymmetry=0;
    end
end


%------------------------------------------
% Construct matrices for the given horizon
%-------------------------------------------

if (isfield(probStruct,'inputblocking') | isfield(probStruct,'deltablocking'))
    opt = Options;
    opt.noConstraintReduction = 1;
    % do note reduce constraints, we postpone that to mpt_blockingMatrices
    % (faster runtime)
    [Matrices]=mpt_constructMatrices(sysStruct,probStruct,opt,setHorizon);
    [Matrices]=mpt_blockingMatrices(Matrices,sysStruct,probStruct,Options);
else
    [Matrices]=mpt_constructMatrices(sysStruct,probStruct,Options,setHorizon);
end

if isinf(-Matrices.W),
    % problem is infeasible
    error('Problem is infeasible!');
end


%compute matrices for move blocking

Pbnd = sysStruct.Pbnd;

ctrlStruct.sysStruct = sysStruct;
ctrlStruct.probStruct = probStruct;

nx = size(sysStruct.A,2);

%solve multiparametric program
if probStruct.norm==2,
    [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mpqp(Matrices,Options);
    details.activeConstraints = activeConstraints;
    ctrlStruct.Pn = Pn;
    ctrlStruct.Fi = Fi;
    ctrlStruct.Gi = Gi;
    ctrlStruct.Ai = details.Ai;
    ctrlStruct.Bi = details.Bi;
    ctrlStruct.Ci = details.Ci;
    details = rmfield(details,'Ai');
    details = rmfield(details,'Bi');
    details = rmfield(details,'Ci');
    ctrlStruct.Pfinal = Phard;
    ctrlStruct.dynamics = ones(1,length(Pn));
    ctrlStruct.details = details;
else
    if strcmp(sysStruct.type,'PWA') | sysStruct.type==1,
        nu=size(sysStruct.B{1},2);
    else
        nu=size(sysStruct.B,2);
    end
    Options.nu = nu*probStruct.N;     % length of optimizer vector without slacks
    [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp(Matrices, Options);
    ctrlStruct.Pn = Pn;
    ctrlStruct.Fi = Fi;
    ctrlStruct.Gi = Gi;
    ctrlStruct.Ai = {};
    for ii=1:length(Pn),
        ctrlStruct.Ai{end+1} = zeros(nx);
    end
    ctrlStruct.Bi = details.Bi;
    ctrlStruct.Ci = details.Ci;
    ctrlStruct.Pfinal = Phard;
    ctrlStruct.dynamics = ones(1,length(Pn));
    ctrlStruct.details = details;
end

ctrlStruct.overlaps = 0;
endtime = cputime;
ctrlStruct.details.runTime = endtime - starttime;

return



%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [constraintsAreSymmetric]=sub1_checkSymmetry(sysStruct,probStruct)
constraintsAreSymmetric=1;

[A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct);

if(isempty(bndA) & (any(isinf(ymin)) | any(isinf(ymax))))
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    disp('mpt_mpqp: Warning: no boundaries specified, regions may be open or unbounded')
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
end
if(all(umin==-umax) & all(ymin==-ymax) &all(dumin==-dumax))
   if(isempty(bndA))
       % nothing happens, constraints are symmetric
   elseif(mod(length(bndA),2)==0)
       %reorder bndA
       bndA2=[];
       bndb2=[];
       k=0;
       limit1=size(bndA,1)/2;
       for i=1:limit1
           oppositefound=0;
           for j=2:size(bndA,1)    
               if(all(bndA(1,:)==-bndA(j,:))&all(bndb(1,:)==bndb(j,:)))
                   remove=j;
                   k=k+1;
                   bndA2(k,:)=bndA(1,:);
                   bndb2(k,:)=bndb(1,:);
                   k=k+1;
                   bndA2(k,:)=bndA(j,:);
                   bndb2(k,:)=bndb(j,:);
                   oppositefound=1;
               end
           end
           if(oppositefound==1)
               bndA(1,:)=[];
               bndb(1,:)=[];
               bndA(remove-1,:)=[];
               bndb(remove-1,:)=[];
           else
               constraintsAreSymmetric=0;
           end
       end
   else
       constraintsAreSymmetric=0;
   end
else
    constraintsAreSymmetric=0;
end
if(constraintsAreSymmetric & ~isempty(bndA))
    sysStruct.bndA=bndA2;
    sysStruct.bndb=bndb2;
end
