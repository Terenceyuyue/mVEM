function ctrlStruct = mpt_iterative(sysStruct,probStruct,Options)
%MPT_ITERATIVE Computes a time-optimal or low-complexity explicit controller for LTI systems
%
% ctrlStruct = mpt_iterative(sysStruct,probStruct)
% ctrlStruct = mpt_iterative(sysStruct,probStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function computes a minimum-time controller for the system defined in 
% "sysStruct" and problem defined in "probStruct".
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------  
% sysStruct        - System structure in the sysStruct format
% probStruct       - Problem structure in the probStruct format
%
%	See the MPT Manual for additional details on the structure format or   
%   consult one of the example systems (e.g. Double_Integator) which were  
%	provided with this package.                                            
%
% Options.lpsolver - Optional: Solver for LPs (see help mpt_solveLP for details)
% Options.abs_tol  - Optional: absolute tolerance
% Options.verbose  - Optional: level of verbosity
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used (see help mpt_init)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% ctrlStruct    - Controller structure with the following fields:
%   Pn,Fi,Gi    - for region Pn(i).H*x <= Pn(i).K control input is U=Fi{i}*x+Gi{i}   
%   Ai,Bi,Ci    - cost associated to each region (x'Aix + Bi*x + Ci)
%                 Note that Ai and Bi are zero matrices, Ci contains the
%                 set distance to the origin
%   Pfinal      - The maximum control invariant set as a polytope object
%   dynamics    - Dynamics active in region Pn(i)
%   details     - Structure with more details about the solution:
%     finalPn   \ 
%     finalFi   - Partition and PWA control law obtained in the final iteration
%     finalGi   /
%     loopCtr   - Number of iterations needed to converge
%     IterStore - Vector which stores the iteration at which each region was
%                 computed, i.e. region j was computed at iteration IterStore(j).
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% "Complexity Reduction of Receding Horizon Control", P. Grieder and M. Morari;
% In the Proceedings of the IEEE Conference on Decision and Control 2003, Maui, Hawaii
%
%
% see also MPT_ITERATIVEPWA
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
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

userSysStruct = sysStruct;
userProbStruct = probStruct;

%emptypoly=polytope;
emptypoly = mptOptions.emptypoly;
if(nargin<3)
    % default options
    Options = [];
end

if(~isfield(Options,'lpsolver'))
    Options.lpsolver=mptOptions.lpsolver;
end
if(~isfield(Options,'debug_level'))
    Options.debug_level=mptOptions.debug_level;
end
if(~isfield(Options,'maxCtr'))
    Options.maxCtr=Inf;
end
if(isfield(Options,'step_size'))
    SUBSET_TOLERANCE=min(Options.step_size*100,mptOptions.step_size);
else
    SUBSET_TOLERANCE=mptOptions.step_size;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'ispwa'),
    % even though this function does not provide a full support for PWA systems,
    % it can be used to compute iterative solution for one fiex dynamics. Since
    % sysStruct will be in PWA format, we need this flag to suppress error
    % messages and to enforce slightly different convergence checks
    Options.ispwa=0;                   
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


lpsolver=Options.lpsolver;
userTset = isfulldim(probStruct.Tset);     % true if user provided the terminal set, false otherwise

if Options.ispwa==0,
    % this routine does not work for "true" PWA systems and for 1-norm case (see the note above)
    if iscell(sysStruct.A)
        error('Sorry, this function handles LTI systems only!');
    end
end

starttime = cputime;

%extract elements from struct
[A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct);

horizon=1;                          % solving 1-step problem
emptypoly = polytope;
probStruct.N = 1;
nR = 0;
if ~isfulldim(probStruct.Tset) & probStruct.tracking==0
    if Options.ispwa,
        % limited support for PWA systems requires the Terminal set to be given
        error('mpt_iterative: Polytopic target set MUST be given for PWA system!');
    end
    nR=1;
    if probStruct.norm~=2 & ~Options.ispwa,
        
        % check if weights are symmetrical and positive definite
        if all(eig(Q)>=0) & all(all(Q==Q')),
            QQ = Q;
        else
            % otherwise use unity weights
            QQ = eye(size(A,1));
        end
        if all(eig(R)>0) & all(all(R==R')),
            RR = R;
        else
            % otherwise use unity weights
            RR = eye(size(B,2));
        end
        [Fi{nR},PLQR] = mpt_dlqr(A,B,QQ,RR);
    else
        [Fi{nR},PLQR] = mpt_dlqr(A,B,Q,R);
    end
    Fi{nR}=-Fi{nR};
    Gi{nR}=zeros(length(umax),1);
end

%construct matrices for the given horizon; matrices are needed to correctly formulate the mpQP problem
if ~Options.ispwa,
    if probStruct.norm~=2,
        origProbStruct = probStruct;
        probStruct.norm=2;
        disp('mpt_iterative: switching to quadratic cost objective (does not affect minimum-time solution properties)');
        if all(eig(Q)>=0) & all(all(Q==Q')),
        else
            % otherwise use unity weights
            if Options.verbose > 0,
                disp('mpt_iterative: 1/Inf problem given, Q not symmetrical and/or < 0, setting them to identity...');
            end
            probStruct.Q = eye(size(A,1));
        end
        if all(eig(R)>0) & all(all(R==R')),
        else
            % otherwise use unity weights
            if Options.verbose > 0,
                disp('mpt_iterative: 1/Inf problem given, R not symmetrical and/or <= 0, setting them to identity...');
            end
            RR = eye(size(B,2));
        end
        [Matrices]=mpt_constructMatrices(sysStruct,probStruct,Options);
        Pinvset = Matrices.Pinvset;
        probStruct = origProbStruct;
    else
        [Matrices]=mpt_constructMatrices(sysStruct,probStruct,Options);
        Pinvset = Matrices.Pinvset;
    end
end

if nR==1,
    Pfinal = Pinvset;
    if ~isfulldim(Pfinal),
        error('Bailing out since no fully dimensional invariant set found.');
    end
    if isempty(bndA),
        Pbnd = emptypoly;
    else
        Pbnd = polytope(bndA, bndb);
    end
    if isfulldim(Pbnd),
        Pfinal = Pfinal & Pbnd;
    end
    Pstore = emptypoly;
    Pstore(1) = Pfinal;
    Fstore{1}=Fi{1};    %store feedback law of region 1
    if isfield(probStruct,'FBgain'),
        Fstore{1} = Fstore{1} - probStruct.FBgain;
    end
    Gstore{1}=Gi{1};    %store feedback law of region 1
    
    invariantTSet=1;
    Pn=Pstore;
else
    Pfinal = probStruct.Tset;
    if isempty(bndA),
        Pbnd = emptypoly;
    else
        Pbnd = polytope(bndA, bndb);
    end
    if Options.ispwa,
        Pstore = emptypoly;
    else
        Pstore = Pinvset;
    end
    Fstore={};
    Gstore={};
    invariantTSet=0;
    Pn=Pstore;
    Fi={};
    Gi={};
end
    
    
maxIterations=100;

%initialize
IterStore(1)=0;     % iteration zero
nx=size(A,1);       % number of states
notConverged=1;
firststep=1;
loopCtr=0;
R_old=0;

if Options.ispwa,
    tmpProbStruct = probStruct;
    tmpProbStruct.Tset = polytope;
    tmpProbStruct.Tconstraint = 2;
    tmpProbStruct.subopt_lev = 0;
    tmpProbStruct.N = 1;
    BMatrices = mpt_constructMatrices(sysStruct, tmpProbStruct, Options);
    if isfield(Options,'ssInfo'),
        ssInfo = Options.ssInfo;
        nx = ssInfo.nx;
        nu = ssInfo.nu;
        ny = ssInfo.ny;
        nPWA = ssInfo.nPWA;
        nbool = ssInfo.nbool;
        ubool = ssInfo.ubool;
        intInfo = ssInfo.intInfo;
    else
        [nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);
    end
end

solveroptions = Options;
solveroptions.verbose = 0;

if statusbar,
    statbar.handle = mpt_statusbar('Computing...');
    statbar.progress = 0;
    Options.verbose = -1;
end

while(notConverged & loopCtr<Options.maxCtr)
    recomputed=0;
    PfOld=Pfinal;
    PnOld=Pn;
    FiOld=Fi;
    GiOld=Gi;
    if Options.verbose>1, 
        disp(['Loop Couter is ' num2str(loopCtr)]);
    end

    if statusbar,
        progress = mod(loopCtr, 50) / 50;
        if isempty(mpt_statusbar(statbar.handle, progress, 0, 0.9)),
            mpt_statusbar;
            error('Break...');
        end     
    end
    
    loopCtr=loopCtr+1;

    if Options.verbose > 0,
        fprintf('Iteration %d\n', loopCtr);
    end
    if ~Options.ispwa,
        tmpProbStruct = probStruct;
        tmpProbStruct.Tset = Pfinal;
        tmpProbStruct.Tconstraint = 2;
        tmpProbStruct.subopt_lev = 0;
        tmpProbStruct.N = 1;
        tmpCtrlStruct = mpt_optControl(sysStruct,tmpProbStruct,solveroptions);
        Pn = tmpCtrlStruct.Pn;
        Fi = tmpCtrlStruct.Fi;
        Gi = tmpCtrlStruct.Gi;
        if probStruct.norm==2,
            activeConstraints = tmpCtrlStruct.details.activeConstraints;
        end
        Pfinal = tmpCtrlStruct.Pfinal;
    else
        [Matrices,mfeas] = mpt_addTset(sysStruct, BMatrices, Pfinal,nx,nu,Options.pwa_index);
        if ~mfeas,
            disp('mpt_iterative: infeasible problem');
            break
        end
        if probStruct.norm==2,
            [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mpqp(Matrices,solveroptions);
        else
            [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,solveroptions);
        end
    end

    if statusbar,
        progress = mod(loopCtr, 50) / 50;
        if isempty(mpt_statusbar(statbar.handle, progress, 0, 0.9)),
            mpt_statusbar;
            error('Break...');
        end     
    end

    Piter{loopCtr} = Pfinal;
    
    if Options.ispwa==0,                    % for pwa case, this step is not needed
        % check if hull is correct
        if(invariantTSet & Options.debug_level>0 & loopCtr>1)
            issubset=le(Piter{loopCtr-1},Piter{loopCtr},Options);
            recomputed=0;
            if (~issubset)
                if Options.verbose>-1,
                    disp('mpt_iterative: ERROR: Terminal Set has shrunk since previous iteration. Recomputing Terminal Set...')
                    disp('It may be possible to fix this by modifying Options.step_size for the mpQP.')
                end
                recomputed=1;
                
                Pfinal2 = hull(Pn);
                
                issubset=le(Piter{loopCtr-1},Pfinal2,Options);
            end
            if(~issubset)
                for ii=1:length(Piter)-1,
                    PiterN{ii}=Piter{ii};
                end
                Piter=PiterN;
                Pfinal=PfOld;
                Pn=PnOld;
                Fi=FiOld;
                Gi=GiOld;
                disp('mpt_iterative ERROR: Terminal Set has shrunk since previous iteration.')
                error('mpt_iterative: Feasible Set is not enlarged !!');
            end
            if(recomputed)
                %recomputed hull and no error was found => store new hull
                if ~exist('Pfinal2','var'),
                    Pfinal2 = hull(Pn);
                end
                Piter{loopCtr}=Pfinal2;
            end
        end
    end
    
    %check whether the current outer hull is identical to the outer hull of the previous iteration
    
    identical = eq(Pfinal,tmpProbStruct.Tset,Options);
    notConverged=~identical; %convergence as soon as the convex hull does not change
    
    if(isfulldim(Pbnd) & notConverged)
        [issubset]=le(Pbnd, Pfinal, Options);
        notConverged=~issubset;
    end
    if Options.ispwa,
        % special handling of the PWA case
        if le(Pfinal, PfOld, Options) & ~firststep,
            if Options.verbose>=1,
                disp('feasible set starts to shrink, aborting');
            end
            notConverged=0;
            Pfinal=PfOld;
            Pn=PnOld;
            Fi=FiOld;
            Gi=GiOld;
        elseif ~le(PfOld,Pfinal,Options)
            if Options.verbose>=1,
                disp('feasible set did not expand, aborting');
            end
            notConverged=0;
            for i=1:length(Pn)
                if(~le(Pn(i),PfOld,Options))
                    IterStore(length(IterStore)+1)=loopCtr;
                    Pstore = [Pstore Pn(i)];
                    Fstore{length(Fstore)+1}=Fi{i};
                    Gstore{length(Gstore)+1}=Gi{i};
                end
            end
            
        else
            for i=1:length(Pn)
                if(~le(Pn(i),PfOld,Options))
                    IterStore(length(IterStore)+1)=loopCtr;
                    Pstore = [Pstore Pn(i)];
                    Fstore{length(Fstore)+1}=Fi{i};
                    Gstore{length(Gstore)+1}=Gi{i};
                end
            end
        end
        %notConverged=0;
        firststep=0;
    else   % LTI system
        %store regions and associated feedback law
        regions_added = 0;
        for i=1:length(Pn)
            if(loopCtr>1 & le(Pn(i),Piter{loopCtr-1},Options)) 
                %controller will never be applied    
            elseif(loopCtr==1 & i==1 & invariantTSet & le(Pn(i),Pstore(1),Options))
            else
                regions_added=1;
                IterStore(length(IterStore)+1)=loopCtr;
                Pstore = [Pstore Pn(i)];
                Fstore{length(Fstore)+1}=Fi{i};
                Gstore{length(Gstore)+1}=Gi{i};
            end
        end
        if probStruct.tracking & ~regions_added
            % no regions contribute to the solution, abort
            break
        end
    end
end

if statusbar,
    mpt_statusbar(statbar.handle, 1);
end

loopCtr=loopCtr-1; %the last iteration step does not "really" count
if userTset & ~Options.ispwa,
    Pstore=Pstore(2:end);
end
if userTset | probStruct.tracking,
    IterStore=IterStore(2:end);
end

ctrlStruct.sysStruct = userSysStruct;
ctrlStruct.probStruct = userProbStruct;
ctrlStruct.Pn = Pstore;
ctrlStruct.Fi = Fstore;
ctrlStruct.Gi = Gstore;
ctrlStruct.Ai = {};
ctrlStruct.Bi = {};
ctrlStruct.Ci = num2cell(IterStore(1:length(Pstore)));
for ii=1:length(Pstore);
    ctrlStruct.Ai{end+1} = zeros(nx);
    ctrlStruct.Bi{end+1} = zeros(1,nx);
end
ctrlStruct.Pfinal = Pfinal;
ctrlStruct.dynamics = ones(1,length(Pstore));
details.finalPn = Pn;
details.finalFi = Fi;
details.finalGi = Gi;
details.loopCtr = loopCtr;
details.IterStore = IterStore;
details.Piter = Piter;
ctrlStruct.details = details;
ctrlStruct.overlaps = 1;

endtime = cputime;
ctrlStruct.details.runTime = endtime-starttime;

if closestatbar,
    mpt_statusbar;
end