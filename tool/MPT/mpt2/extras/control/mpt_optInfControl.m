function ctrlStruct=mpt_optInfControl(sysStruct,probStruct,Options)
%MPT_OPTINFCONTROL Solves the infinite-time constrained optimal control problem for LTI systems
%
% ctrlStruct=mpt_optInfControl(sysStruct,probStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Explicitly solves the problem...
%         
% min_{u(0),...} [sum_{i=0}^infty   x(i)'Qx(i) + u(i)'Ru(i)]
% subject to  u(i)  \in U   for i=0,1,....
%             x(i)  \in X   for i=1,2,....
%   
% The resulting input sequence is still of finite dimension. As soon as the 
% state enters a region of the state space where the Riccati LQR feedback law
% satisfies the system constraints for all time, no further inputs are computed.
% This script combines multiparametric-programming techniques with reachability 
% analysis.
%  
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct            - System structure in the sysStruct format
% probStruct           - Problem structure in the probStruct format
%
%	See the MPT Manual for additional details on the structure format or   
%       consult one of the example systems (e.g. Double_Integator) which were  
%	provided with this package.                                            
%
% Options.verbose      - level of verbosity
% Options.lpsolver     - which LP solver to use (help mpt_solveLP)
% Options.qpsolver     - which QP solver to use (help mpt_solveQP)
% Options.step_size    - length of step over a facet
% Options.maxHorizon   - The maximum horizon which is used for computations;
%                        Leave empty to obtain the whole set.
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
% ctrlStruct    - controller structure with the following fields:
%   Pn,Fi,Gi    - for region Pn(i).H*x <= Pn.K(i) computed input is U=Fi{i}*x+Gi{i}   
%   Ai,Bi,Ci    - cost associated to each region (x'Aix + x'Bi + Ci)
%   Pfinal      - Defines the feasible state space partition (i.e. union of
%                 all regions) as Phard.H*x<=Phard.K
%   dynamics    - dynamics active in region Pn(i)
%   details     - contains additional information:
%      regionHorizon      - Vector containing the number of steps required for each
%                           region to reach the control invariant set. 
%		                    (i.e. dimension of the associated input sequence)
%      activeConstraints  - Cell Array which stores the active constraints 
%                           of the optimizer in each region.
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
%
% "Computation of the Constrained Infinite Time Linear Quadratic Regulator",
% P. Grieder, F. Borrelli, F. Torrisi, M. Morari; In the proceedings of the
% American Control Conference (ACC) 2003, Denver, Colorado
%
% see also MPT_CONTROL, MPT_OPTCONTROL, MPT_ITERATIVE, MPT_ITERATIVEPWA

% Copyright is with the following author(s):
%
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
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

error(nargchk(2,3,nargin));

zero_tol = 1e-14;     % everything below this value is considered zero

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if ~isfield(sysStruct,'verified'),
    verOpt.verbose=1;
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end

if ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    probStruct=mpt_verifyProbStruct(probStruct,verOpt);
end

if ~(sysStruct.type==0 | strcmp(sysStruct.type,'LTI'))
    error('mpt_optInfControl can only handle LTI systems! Please use mpt_optInfControlPWA for PWA systems');
end

if(nargin<3)
    Options = [];
end

if(~isfield(Options,'lpsolver'))
    Options.lpsolver=mptOptions.lpsolver;
end
if(~isfield(Options,'qpsolver'))
    Options.qpsolver=mptOptions.qpsolver;
end
if(~isfield(Options,'maxHorizon'))
    Options.maxHorizon = Inf; %add set constraint as terminal set constraint
end
if(~isfield(Options,'debug_level'))
    Options.debug_level=mptOptions.debug_level;
end
if ~isfield(Options,'step_size'),
    Options.step_size=mptOptions.step_size;
end
if(~isfield(Options,'verbose'))
    Options.verbose=mptOptions.verbose;
end

if ~isfield(Options, 'statusbar'),
    Options.statusbar = 0;
end
if ~isfield(Options, 'closestatbar')
    Options.closestatbar = 1;
end
closestatbar = Options.closestatbar;
Options.closestatbar = 0;
statusbar = Options.statusbar;
if statusbar,
    Options.statusbar = 0;
end

if probStruct.N<Inf
    disp('mpt_optInfControl: probStruct.N should be Inf in this function');
    probStruct.N=Inf;
end
Options.noConstraintReduction=1;    % needed for mpt_constructMatrices since N=inf

lpsolver=Options.lpsolver;       
useSymmetry=probStruct.useSymmetry;
maxHorizon=Options.maxHorizon; 
DEBUG=Options.debug_level;

if isfield(sysStruct,'noise'),
    if mpt_isnoise(sysStruct.noise),
        error('mpt_optInfControl: Additive disturbances not supported in infinite-time control!');
    end
end

if isfield(sysStruct,'Aunc'),
    error('mpt_optInfControl: Parametric uncertainty not supported in infinite-time control!');
end

starttime = cputime;

[A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct);
nu=size(B,2);
nx=length(A);

EMPTY_ROW_TOL=mptOptions.abs_tol;     % Tolerance for declaring that row is empty
CONSTR_TOL = mptOptions.rel_tol;      % Tolerance for declaring that constr. is redundant

stepSize = Options.step_size;         % How far to go from border of two regions
%----------Debug tolerance-------------------------------------------------------------
if(DEBUG==2)
    TOLERANCE=0;
else
    TOLERANCE=stepSize*10;
end
%-------------------------------------------------------------------------------------------------
if Options.qpsolver==1,
    options=optimset(optimset('quadprog'),'Display','off','LargeScale','off');  %set quadprog options
else
    options = [];
end
horizon=1;       %initial horizon (reachability) is 1
tmpProbStruct = probStruct;
tmpProbStruct.N=1;
Options.includeLQRset = 0;
% transform problem into the form: min (U'HU + xFU + x'Yx)  s.t. G*U<=W+E*x
[G,W,E,H,F,Y,Cf,Cx,Cc,probStruct.useSymmetry,bndA,bndb,Pinvset]=mpt_constructMatrices(sysStruct,tmpProbStruct,Options);

Matrices.G=G;
Matrices.W=W;
Matrices.E=E;
H=(H+H')/2;
Matrices.S=E+G*inv(H)*F';
S = Matrices.S;
Matrices.H=H;
Matrices.F=F;
Matrices.Y=Y;
Matrices.Cx=Cx;
Matrices.Cf=Cf;
Matrices.Cc=Cc;
Matrices.bndA=bndA;
Matrices.bndb=bndb;

nR = 1;                                            % number of regions => initialize to one
activeConstraints{nR}  = [];           % no active constraints in the control-invariant set

constraintStorage      = [];              % initializes the storage structure for all active constraints                                                      

Fi={};      % control law                                       => optimal control in region k is defined by u=Fi{k}*x+Gi{k}
Gi={};      % control law
hardA=[];   % outer hull
hardb=[];   % outer hull
nHard=0;
xB=[]; %structure for storing all the points xBeyond
xBRegion=[]; %structure for storing region associated to point
Pn = polytope;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   Compute control invariant set for the unconstrained LQR feedback law
%   ATTENTION: G,W must be previously computed using horizon 1 !!
%1)
%get LQR feedback law
[Fi{nR},PLQR] = mpt_dlqr(A,B,Q,R);
Fi{nR}=-Fi{nR};
Gi{nR}=zeros(length(umax),1);
%2)
%compute unconstrained target region region for reachability (constraints must hold at t=0)

Hinv = inv(H);
uA = (G*-Hinv*F')-E;  %U=-Hinv*F' is feedback law  ; ie G*U-E*x<=W ; bndA*x<=bndb
ub = W;
totA  = [uA ; bndA;C;-C];
totb  = [ub ; bndb*1e6;ymax;-ymin];  %10000 is an arbitrary value; it is chosen to encompass the entire unconstrained region
X = polytope(totA, totb);
localOptions = Options;
localOptions.abs_tol = CONSTR_TOL;
[PFinal,tstar,fd,isemptypoly] = mpt_infset(A+B*Fi{nR},X,100,sysStruct.noise,localOptions);
if(isemptypoly)
    error('mpt_optInfControl: There is no positive invariant region for the LQR controller')
end


%3)
%compute unconstrained region (constraints must hold at t'=t+1)
totA  = [uA ; bndA;];
totb  = [ub ; bndb];  
X = polytope(totA, totb);

[Pcr,tstar,fd,isemptypoly] = mpt_infset(A+B*Fi{nR},X,100,sysStruct.noise,Options);
if(isemptypoly)
    error('mpt_optInfControl: There is no positive invariant region for the LQR controller')
end

% Pcr = Pinvset;
% PFinal = Pinvset;

Pn = [Pn Pcr];

clear keptrows isemptypoly
%initialize xBeyond=> stores points which are only just outside the region

for i=1:nconstr(Pn(nR)),
    xBeyond{nR}{i}=[];                    
end

%initialize value function J(x)=x'Ai{i}x+x'Fi{i}+Ci{i}
Ai{1} = PLQR;
Bi{1} = zeros(1,nx);
Ci{1} = 0;

%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%%%%%%%%%%%%%%%%%%%%%%%
%%                    %%
%% DO THE EXPLORATION %%
%%                    %%
%%%%%%%%%%%%%%%%%%%%%%%%

region=1;                                                 % start exploration from region 1
horizon=1;                                                 % first check reachability for horizon 1
keptrows{region}=2*ones(1,nconstr(Pn(nR)));         % set to 2 so that all borders are examined
regionHorizon(region)=0;                                                 % region 1 needs 0 steps to enter region 1 s
oldNrOfRegions=1;                                                    % only one region so far
maxRegion=1;                                                     % only one region so far
completedRegions=0;                                                     %no regions completed so far
checkedFacet{region}=[];                                                     %no facets checked so far

if statusbar,
    statbar.handle = mpt_statusbar('Computing...');
    progress = 0;
    Options.verbose = -1;
end


while (region<=nR & horizon<=maxHorizon),  
    borderVec =1: nconstr(Pn(region));
    if(region>maxRegion)
        checkedFacet{region}     =   [];   %initialize neighbor storage
        maxRegion                =   region;
    end 

    progress = mod(horizon, 20) / 20;
    if statusbar,
        if isempty(mpt_statusbar(statbar.handle, progress, 0, 0.9)),
            mpt_statusbar;
            error('Break...');
        end
    end

    %----------------------------------------------------------------------------------
    % Remove Facets from vector that do not need to be checked...
    %----------------------------------------------------------------------------------
    if(horizon==regionHorizon(region))
        borderVec(find(keptrows{region}==0))=[]; %do not check facets which originated from reachability restrictions 
    end   
    if(horizon>regionHorizon(region))
        borderVec(find(keptrows{region}==1))=[]; %do not check facets which originated from constraints restrictions 
    end
    if(~isempty(checkedFacet{region}) & ~isempty(borderVec))
        borderVec(find(ismember(borderVec,checkedFacet{region})==1))=[];   %only check borders with no adjacent neighbour
    end
    
    
    for border=borderVec,    %only search borders which originate from input/output constraints (instead of reachability)
        [Hp,Kp] = double(Pn(region));
        if(isempty(xBeyond{region}{border}))
            % Create a point close to the border  => do not compute this point more than once
            [xBorder,RBorder]=facetcircle(Pn(region),border,Options);  
            
            xBeyond{region}{border}=xBorder+stepSize*Hp(border,:)';
            
            insideBound = 1;
            %check if point is still within the bounded state-space
            if(~isempty(bndA)&~all(bndA*xBeyond{region}{border}-bndb<=0))
                insideBound   =   0;
            end
            
            nMatch=0;
            if(insideBound)
                isinOpt.abs_tol = 0;
                isinOpt.fastbreak = 1;
                [isthere,iwhere] = isinside(Pn, xBeyond{region}{border}, isinOpt);
                for iii=1:length(iwhere),
                    if iwhere(iii)~=region,
                        nMatch = nMatch+1;
                        % neighbor=iii;
                    end
                end
            end
            
        else
            %Point close to border was alerady computed
            insideBound = 1;
            nMatch      = 0;
        end
        
        
        if (nMatch>=1 | ~insideBound),
            if nMatch>1,
                disp('mpt_optInfControl:  Polyhedra should not overlap');
                %this sometimes happens due to numerical inconsistencies and can "tolerated" 
                %it's best to check for overlapping regions anyway once the algorithm has completed
            end
            checkedFacet{region} =   [checkedFacet{region} border];  %add current border to vector which stores the explored adjacent regions
        else  
            
            %==============================================%
            % Solve optimization problem for point xBeyond %
            %==============================================%
            [zopt,lambda,how,exitflag,Vz]=mpt_solveQP(H,zeros(horizon,length(umax)),G,W+S*xBeyond{region}{border},[],[],[],Options.qpsolver,options);
            ii=find(abs(G*zopt-W-S*xBeyond{region}{border})<=zero_tol);   % active constraints
            kk=1;
            while(mpt_iscombequal(ii,activeConstraints{region}) & kk<=100 & exitflag>0 & (isempty(bndA)|all(bndA*xBeyond{region}{border}-bndb<=0)))
                %this is not possible, unless numerical errors occur
                %remedy this by increasing the step size until a new region is reached
                
                xBeyond{region}{border}=xBorder+stepSize*10^(kk/4)*Hp(border,:)';
                if(kk<5 & DEBUG==2)
                    xB(size(xB,1)+1,:)= xBeyond{region}{border}';
                    xBRegion(size(xB,1),:)=region;
                end
                if(~isempty(bndA)&~all(bndA*xBeyond{region}{border}-bndb<=0))
                    %outside bounds
                    kk=10000;
                else
                    [zopt,lambda,how,exitflag,Vz]=mpt_solveQP(H,zeros(horizon,length(umax)),G,W+S*xBeyond{region}{border},[],[],[],Options.qpsolver,options);
                    ii=find(abs(G*zopt-W-S*xBeyond{region}{border})<=zero_tol);   % active constraints
                end
                kk=kk+1;
            end
            %reset xBeyond to initial value
            
            xBeyond{region}{border}=xBorder+stepSize*Hp(border,:)';
            
            if(exitflag>0)
                [constraintStorage,noNewFacetFound]=mpt_checkStorage(constraintStorage,ii);      %function that checks whether the active constraints were already covered
            end
            if((exitflag<=0) | noNewFacetFound)
                %set was already covered 
                if(exitflag==0)
                    error('mpt_optInfControl: Maximum number of iterations in QP solver exceeded... results could be false.')
                end
                if(exitflag<=0 & kk==1)
                    checkedFacet{region} =   [checkedFacet{region} border]; 
                    nHard=nHard+1;
                    
                    hardA(nHard,:)=Hp(border,:);
                    hardb(nHard,:)=Kp(border,:);
                end
            else
                nR=nR+1; % new region
                activeConstraints{nR}=ii;   %store active constraints

                %get active subset of region for the current horizon settings  
                [Pcr,Fi{nR},Gi{nR}] = sub6_computelaw(ii,nu*horizon,nx,Matrices);
                opt.lpsolver=Options.lpsolver;
                opt.abs_tol=CONSTR_TOL;

                [Pret,keptrows{nR},feasible] = mpt_getReachSubset(Pcr,PFinal,A,B,Fi{nR},Gi{nR},horizon,opt);
                
                
                if(feasible)
                    if (mod(region,20)==0)
                        if Options.verbose>0,
                            disp(sprintf('Region : %d/%d      ',region,nR));
                        end
                    end
                    checkedFacet{region} =   [checkedFacet{region} border];  %add current border to vector which stores the explored adjacent regions
                    regionHorizon(nR)    =   horizon;
                    Pn(nR) = Pret;
                    
                    %initialize xBeyond
                    
                    for i=1:nconstr(Pret),
                        xBeyond{nR}{i}=[];
                    end
                    constraintStorage=mpt_addToStorage(activeConstraints{nR},constraintStorage,horizon);
                    % Add Symmetric Region
                    if(useSymmetry)
                        nR=nR+1; % new region
                        
                        %compute active constraints for mirror region
                        clear ii2;
                        for j=1:length(ii)
                            if(mod(ii(j),2)==1)
                                ii2(j) = ii(j)+1;
                            else
                                ii2(j) = ii(j)-1;
                            end
                        end
                        
                        activeConstraints{nR}=ii2'; 
                        keptrows{nR}=keptrows{nR-1};
                        Pn(nR) = -Pn(nR-1);
                        Fi{nR} = Fi{nR-1};
                        Gi{nR} = -Gi{nR-1};
                        
                        if(exist('regionHorizon','var'))
                            regionHorizon(nR)=horizon;
                        end
                        
                        
                        if(iscell(xBeyond))
                            %initialize xBeyond in the infinite horizon case
                            for i=1:nconstr(Pret),
                                xBeyond{nR}{i}=[];
                            end
                        end
                        
                        constraintStorage=mpt_addToStorage(activeConstraints{nR},constraintStorage,horizon);
                    end
                else
                    %empty region =>  reduce counter so that the region which was just updated will be overwritten
                    checkedFacet{nR} = [];
                    activeConstraints{nR}=[];
                    nR=nR-1;
                end
            end
        end % nMatch==0
        if(DEBUG==2)
            xB(size(xB,1)+1,:)= xBeyond{region}{border}'; % For each border, store explored point
            xBRegion(size(xB,1),:)=nR;
        end    
    end % border
    
    if(length(checkedFacet{region})==nconstr(Pn(region)))
        %all facets of this region have been successfully explored
        if(isempty(find(completedRegions==region)))
            %do not check this region again
            completedRegions = [completedRegions region];
        end
    end
    
    region=region+1;
    while(~isempty(find(completedRegions==region)))
        region=region+1; %increase start region counter
    end
    
    
    if(region>nR)
        %store value function
         Matrices.H=Matrices.H/2; %to extract proper cost without 1/2
         for i=(oldNrOfRegions+1):nR,
             Ai{i}= Fi{i}'*Matrices.H*Fi{i}+ F*Fi{i} + Matrices.Y;
             Bi{i}= (2*Gi{i}'*Matrices.H+Matrices.Cf)*Fi{i}+Gi{i}'*Matrices.F'+Matrices.Cx;
             Ci{i}= Gi{i}'*Matrices.H*Gi{i}+Matrices.Cf*Gi{i}+Matrices.Cc;
         end

        if(oldNrOfRegions==nR)
            %number of regions hasn't changed with increased horizon => end reached
            if Options.verbose>1,
                disp('Number of regions hasn''t changed from N to N+1');
            end
            [Pn,Fi,Gi,activeConstraints]=mpt_adjustCellSize(Pn,Fi,Gi,activeConstraints);
            if ~exist('PHard','var')
                % no such varaiable
                PHard = polytope(hardA, hardb);
            end
            ctrlStruct.sysStruct = sysStruct;
            ctrlStruct.probStruct = probStruct;
            ctrlStruct.Pn = Pn;
            ctrlStruct.Fi = Fi;
            ctrlStruct.Gi = Gi;

            %The resulting value function is V(x)=x'*Ai{i}*x+Bi{i}*x+Ci{i} if x \in P(i)
            ctrlStruct.Ai = Ai;
            ctrlStruct.Bi = Bi;
            ctrlStruct.Ci = Ci;
         
            ctrlStruct.Pfinal = PHard;
            ctrlStruct.dynamics = ones(1,length(Pn));
            details.regionHorizon = regionHorizon;
            details.activeConstraints = activeConstraints;
            ctrlStruct.details = details;
            ctrlStruct.overlaps = 0;
            endtime = cputime;
            ctrlStruct.details.runTime = endtime - starttime;
            
            if statusbar,
                mpt_statusbar(statbar.handle, 1);
                mpt_statusbar;
            end
            
            return
        end
        %increase output horizon 
        horizon = horizon + 1;
        if Options.verbose>0,
            disp(sprintf('------------     New Horizon %d     ------------',horizon));
        end
        %construct new weight (H) and constraint (G,W,E) matrices for the increased horizon
        tmpProbStruct=probStruct; 
        tmpProbStruct.N=horizon;
		[G,W,E,H,F,Y,Cf,Cx,Cc,probStruct.useSymmetry,bndA,bndb,Pinvset]=mpt_constructMatrices(sysStruct,tmpProbStruct,Options);
        Matrices.G=G;
        Matrices.W=W;
        Matrices.E=E;
        H=(H+H')/2;
        Matrices.S=E+G*inv(H)*F';
        S = Matrices.S;
		Matrices.H=H;
		Matrices.F=F;
		Matrices.Y=Y;
		Matrices.Cx=Cx;
		Matrices.Cf=Cf;
		Matrices.Cc=Cc;
        
        %start again from the beginning
        region=sub_mpt_getActiveMinimum(completedRegions)+1;
        oldNrOfRegions=nR;
    end
end % END EXPLORATION

[Pn,Fi,Gi,activeConstraints]=mpt_adjustCellSize(Pn,Fi,Gi,activeConstraints);
%-----------------------------------------------------------------
% -Remove false outer bounds
% -Check if the outer hull contains empty regions
%-----------------------------------------------------------------

if statusbar,
    mpt_statusbar(statbar.handle, 1);
end

if(DEBUG==1|DEBUG==2)
    [hardA,hardb,NoSolutionToPoint]=sub1_fixouterhull(hardA,hardb,xB,Pn,lpsolver,TOLERANCE,stepSize,xBRegion);
end
PHard = polytope(hardA, hardb);

ctrlStruct.sysStruct = sysStruct;
ctrlStruct.probStruct = probStruct;
ctrlStruct.Pn = Pn;
ctrlStruct.Pn = Pn;
ctrlStruct.Fi = Fi;
ctrlStruct.Gi = Gi;

 %store value function
 Matrices.H=Matrices.H/2; %to extract proper cost without 1/2
 for i=(oldNrOfRegions+1):nR,
     Ai{i}= Fi{i}'*Matrices.H*Fi{i}+ F*Fi{i} + Matrices.Y;
     Bi{i}= (2*Gi{i}'*Matrices.H+Matrices.Cf)*Fi{i}+Gi{i}'*Matrices.F'+Matrices.Cx;
     Ci{i}= Gi{i}'*Matrices.H*Gi{i}+Matrices.Cf*Gi{i}+Matrices.Cc;
 end
         
%The resulting value function is V(x)=x'*Ai{i}*x+Bi{i}*x+Ci{i} if x \in P(i)
ctrlStruct.Ai = Ai;
ctrlStruct.Bi = Bi;
ctrlStruct.Ci = Ci;

ctrlStruct.Pfinal = PHard;
ctrlStruct.dynamics = ones(1,length(Pn));
details.regionHorizon = regionHorizon;
details.activeConstraints = activeConstraints;
ctrlStruct.details = details;
ctrlStruct.overlaps = 0;
endtime = cputime;
ctrlStruct.details.runTime = endtime - starttime;

if statusbar,
    mpt_statusbar;
end

return



%---------------------------------------------------------------------
% SUBFUNCTION 1
%-----------------------------------------------------------------------
function [hardA,hardb,NoSolutionToPoint2]=sub1_fixouterhull(hardA,hardb,xB,Pn,lpsolver,TOLERANCE,stepSize,xBRegion);
%------------------------------------------------------------------------------------
% Function which:
% -checks if all the regions lies inside the computed outer hull,
%  if not, the bound which one are outside will be removed
%- "TOLERANCE" is decided from the debug features 
% (DEBUG = 1, tolerance is equal to stepsize
%  DEBUG = 2, tolerance is strictly zero     )
%------------------------------------------------------------------------------------
for region=1:length(Pn)
    [x{region},R{region}]=chebyball(Pn(region));
    remove=find(hardA*x{region}-hardb>0);
    if(~isempty(remove))
        %-----------------------------------------------------------
        % Remove the false outer bounds from the boundary structure
        %------------------------------------------------------------
        hardA(remove,:)=[];
        hardb(remove,:)=[];
    end
end

%------------------------------------------------------------------
% Check if there exist points which are not associated with a region
% inside the outer hull
%------------------------------------------------------------------
RemoveRegion=find(xBRegion>length(Pn));
xBRegion(RemoveRegion)=[]; %remove points not associated to region anymore
xB(RemoveRegion,:)=[];
NoSolutionToPoint=[];
isinOpt.abs_tol = 2*TOLERANCE;
isinOpt.fastbreak = 1;
if(~isempty(xB))
    for i=1:size(xB,1)
        xBeyond=xB(i,:)';
        index=xBRegion(i);
        FoundInsideRegion=0;
        if isinside(Pn(index),xBeyond,isinOpt),
            FoundInsideRegion=1;
            if(~all(hardA*xBeyond-hardb<=TOLERANCE) & TOLERANCE==0)
                remove=find(hardA*xBeyond-hardb>0);
                if(~isempty(remove))
                    %-----------------------------------------------------------
                    % Remove the false outer bounds from the boundary structure
                    %------------------------------------------------------------
                    hardA(remove,:)=[];
                    hardb(remove,:)=[];
                end
            end
        else
            for regions=1:length(Pn)
                if isinside(Pn(regions), xBeyond, isinOpt),
                    FoundInsideRegion=1;
                   if(~all(hardA*xBeyond-hardb<=TOLERANCE) & TOLERANCE==0)
                        remove=find(hardA*xBeyond-hardb>0);
                        if(~isempty(remove))
                            %-----------------------------------------------------------
                            % Remove the false outer bounds from the boundary structure
                            %------------------------------------------------------------
                            hardA(remove,:)=[];
                            hardb(remove,:)=[];
                        end
                    end
                end
            end
        end
        if(FoundInsideRegion==0)
            if(all(hardA*xBeyond-hardb<=TOLERANCE))
                NoSolutionToPoint(size(NoSolutionToPoint,1)+1,:)=xBeyond';
            end
        end
    end
end
%---------------------------------------------------------------
% If a point without solution is found, check if the gap is big,
% otherwise it is ok
%-------------------------------------------------------------
NoSolutionToPoint2=[];
isinOpt.abs_tol = 0;
isinOpt.fastbreak = 1;
if(~isempty(NoSolutionToPoint)&TOLERANCE==0)
    nx=size(xBeyond,1);
    for i=1:size(NoSolutionToPoint,1)
        nMatch=0;
        for k=1:nx
            Border=zeros(1,nx);
            Border(k)=1;
            xBeyond1=NoSolutionToPoint(i,:)'+stepSize*100*Border';
            xBeyond2=NoSolutionToPoint(i,:)'-stepSize*100*Border';
            %---------------------------------------------------
            % Check if point is inside another polyhedra
            %---------------------------------------------------
            [dummy1, inwhich] = isinside(Pn, xBeyond1, isinOpt);
            nMatch1 = length(inwhich);
            [dummy1, inwhich] = isinside(Pn, xBeyond2, isinOpt);
            nMatch2 = length(inwhich);
            if(nMatch1==0|nMatch2==0)
                nMatch=0;
            end
        end
        if(nMatch==0)
            NoSolutionToPoint2=[NoSolutionToPoint(i,:);NoSolutionToPoint2];
        end
    end
end


if(~isempty(NoSolutionToPoint2))        
    disp('****************************************************************************************')
    disp('mpt_optInfControl: Points where QP have no feasible solution inside the hull:')
    disp(num2str(NoSolutionToPoint2));
    disp('****************************************************************************************')
end



%---------------------------------------------------------------------
% SUBFUNCTION 2
%-----------------------------------------------------------------------
function returnInt  =   sub_mpt_getActiveMinimum(inputVec)
%
%+++++++++++++++++++++++++++++++++++++++++++++++
% function returns smallest number directly connected to 1
% eg function returns 5 for [1 2 3 4 5 9 81]
% (C) Pascal Grieder grieder@aut.ee.ethz.ch
%+++++++++++++++++++++++++++++++++++++++++++++++
if(isempty(inputVec))
    returnInt  =   0;
    return;
end

returnInt   =   min(inputVec);

for i=1:length(inputVec)
    if(~isempty(find(inputVec==returnInt+1)))
        returnInt=returnInt+1;
    else
        break
    end
end

return

%---------------------------------------------------------------------
% SUBFUNCTION 3
%-----------------------------------------------------------------------
function [constraintStorage,noNewFacetFound]=mpt_checkStorage(constraintStorage,activeConstraint)   
%
%-----------------------------------
%   DESCRIPTION:
%-----------------------------------
% Check if region (ie active constraints) already exists
% Internal routine
%
%-----------------------------------
%    INPUT:
%-----------------------------------   
%  constraintStorage:           storage structure
%  activeConstraint:            constraint to be added to structure
%
%-----------------------------------
%   OUTPUT:
%-----------------------------------
%  constraintStorage:           storage structure
%  noNewFacetFound:             combination/activeConstraints already stored
%
%(C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich, grieder@control.ee.ethz.ch
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


if(isempty(activeConstraint))
     noNewFacetFound = 1;
     return
end

search=1;
noNewFacetFound = 0;

try, if(isempty(constraintStorage{length(activeConstraint)}{min(activeConstraint)}))
        search=0;
      end
catch,
        search=0;
end

if(search)
    storageSegment =  constraintStorage{length(activeConstraint)}{min(activeConstraint)};
    stLength       =  length(storageSegment);
    i=1;
    while(i<=stLength)
       if(length(activeConstraint)==length(storageSegment{i}))
        if(ismember(activeConstraint, storageSegment{i}))
            %combination already stored
             noNewFacetFound = 1;
             return
        end
       end
        i=i+1;
    end
end%if search





%---------------------------------------------------------------------
% SUBFUNCTION 4
%-----------------------------------------------------------------------
function [Pn,Fi,Gi,activeConstraints] = mpt_adjustCellSize(Pn,Fi,Gi,activeConstraints)
%
%-----------------------------------
%   DESCRIPTION:
%-----------------------------------
% this function removes "empty" entries from the cell arrays
% Internal routine.
%
%-----------------------------------
%    INPUT:
%-----------------------------------   
% Pn                 - polytope array defining polyhedral partition
% Fi,Gi              - cell arrays defining the PWA control law
% activeConstraints  - cell array containing active constraints for each region
%
%-----------------------------------
%   OUTPUT:
%-----------------------------------
% Pn                 - polytope array defining polyhedral partition
% Fi,Gi              - cell arrays defining the PWA control law
% activeConstraints  - cell array containing active constraints for each region
%
%(C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich, grieder@control.ee.ethz.ch
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

ctr=1;
PPn = polytope;
for i=1:length(Pn)
    if isfulldim(Pn(i)),
        PPn = [PPn Pn(i)];
        GGn{ctr}=Gi{i};
        FFn{ctr}=Fi{i};
        constr{ctr}=activeConstraints{i};
        ctr=ctr+1;
        end
end
clear Pn Fi Gi activeConstraints
Pn=PPn;
Fi=FFn;
Gi=GGn;
activeConstraints=constr;



%---------------------------------------------------------------------
% SUBFUNCTION 5
%-----------------------------------------------------------------------
function [constraintStorage]=mpt_addToStorage(activeConstraint,constraintStorage,horizon)   
%
%-----------------------------------
%   DESCRIPTION:
%-----------------------------------
% adds the set of active constraints to the storage structure
% Internal routine
%
%-----------------------------------
%    INPUT:
%-----------------------------------   
% activeConstraints
% constraintStorage
% horizon
%
%-----------------------------------
%   OUTPUT:
%-----------------------------------
% constraintStorage
%
%(C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich, grieder@control.ee.ethz.ch
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

try, 
    if(isempty(constraintStorage{length(activeConstraint)}{min(activeConstraint)}))
        segmentLength=0;
    else
        segmentLength =  length(constraintStorage{length(activeConstraint )}{min(activeConstraint)});
    end
catch,
    segmentLength=0;
end


if(length(activeConstraint )>0)
    constraintStorage{length(activeConstraint)}{min(activeConstraint)}{segmentLength+1}= [activeConstraint'];
end




%---------------------------------------------------------------------
% SUBFUNCTION 6
%-----------------------------------------------------------------------
function [Pcr,Fcontrol,Gcontrol] = sub6_computelaw(ii,nu,nx,Matrices)
%----------------------------------------------------------------
% Function which compute control law and region definition
%---------------------------------------------------------------
nii=length(ii);
if nii==0,
    %---------------------------------
    % Control law
    %----------------------------------
    Fcontrol=zeros(nu,nx); 
    Gcontrol=zeros(nu,1);  
    %----------------------------------- 
    % Region definiton
    %-----------------------------------
    crA=[Matrices.bndA; -Matrices.S];      
    crb=[Matrices.bndb;  Matrices.W];      
else
    Gt=Matrices.G(ii,:);
    if nii>nu,
        disp('mpt_mpqp: Degeneracy')
        %----------------------------
        %Deal with primal degeneracy
        %---------------------------------
        [Gt,keptrows]=mpt_getFullRankSubset(Gt,1);
        if(iscell(Gt))
            Pcr=polytope;
            ctr=0;
            for k=1:length(keptrows)
                [Ptemp,Ft,Gt] = sub6_computelaw(ii(keptrows{k}),nu,nx,Matrices);
                if(isfulldim(Ptemp))
                   if(Pcr>=Ptemp)
                       %region is already covered
                   else
                       ctr=ctr+1;
                       Pcr=[Pcr Ptemp];
                       Fc{ctr}=Ft;
                       Gc{ctr}=Gt;
                   end
                end
            end
            Pcr = envelope(Pcr); %merge pieces
            %%Pcr = set(Pcr,'keptrows',length(Matrices.bndb)+(1:nconstr(Pcr)));
            Fcontrol=Fc{1};
            Gcontrol=Gc{1};
            return
        else
           ii=ii(keptrows);
        end
    end
    
    Wt=Matrices.W(ii,:);
    St=Matrices.S(ii,:);
    Matrices.Hinv = inv(Matrices.H);
    GHG=inv(Gt*Matrices.Hinv*Gt');
    GHGinv=inv(Gt*Matrices.Hinv*Gt');
    tmat=Matrices.Hinv*Gt'*GHGinv;
    %-----------------------------------
    % Control law
    %-----------------------------------
    Fcontrol=tmat*St;      
    Gcontrol=tmat*Wt;      
    %-----------------------------------
    % Region definition
    %----------------------------------
    crA=[Matrices.bndA;  Matrices.G*Fcontrol-Matrices.S;  GHGinv*St];   
    crb=[Matrices.bndb; -Matrices.G*Gcontrol+Matrices.W; -GHGinv*Wt];   
    crA(length(Matrices.bndb)+ii,:)=[];       %remove active constraints from Gz<=W+Sx, since the inequality holds by definition
    crb(length(Matrices.bndb)+ii,:)=[];       %remove active constraints from Gz<=W+Sx, since the inequality holds by definition
end

% Switch to the uopt description: U=Z-Hinv*F'*x; %
Fcontrol=Fcontrol-Matrices.Hinv*Matrices.F';
Pcr = polytope(crA, crb);
