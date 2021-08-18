function ctrlStruct=mpt_iterativePWA(sysStruct,probStruct,Options)
%MPT_ITERATIVEPWA Computes a time-optimal or low-complexity explicit controller for PWA systems
%
% ctrlStruct = mpt_iterativePWA(sysStruct,probStruct),
% ctrlStruct = mpt_iterativePWA(sysStruct,probStruct,Options),
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the solution of a minimum-time problem for a given PWA system
%       x(k+1) = A_i x(k) + B_i u(k) + f_i
%       y(k)   = C_i x(k) + D_i u(k) + g_i
%       for i such that guardX(i) x(k) + guardU(i) u(k) <= guardC(i)  
%   s.t.
%       (ymin, ymax, umin, umax, dumin, dumax)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct           - System structure in the sysStruct format
% probStruct          - Problem structure in the probStruct format
%
% Options.iterative=0    
%       Use reduced-switching policy (0=no, 1=yes)
% Options.maxiterations=100
%       Maximum number of iterations
% Options.lyapunov_type
%       a string which denotes which Lyapunov function should be computed to
%       testify stability of one-step controllers:
%         'any'  - first a PWQ function is computed, if it fails, a PWA
%                  Lyapunov function will be computed as well
%         'pwq'  - only PWQ Lyapunov function
%         'pwa'  - only PWA Lyapunov function
%         'none' - do not compute any Lyapunov function; no closed-loop
%                  stability is guaranteed
% Options.verbose
%       Level of verbosity (see help mpt_init for more details)
% Options.PWA_savemode=0
%       Saving intermediate results in fesibility iterations.
% Options.PWA_warmend=0
%       Loads an intermediate result and finishes the controller calculations.
%           (without trying to extend the feasible region).
% Options.PWA_loadn=0
%       Loads intermediate result from iteration n. (0=Uses last iteration.)
% Options.PWA_savefile='PWA_save'
%       Specifies filename prefix for result files.
%           Filename is appended with iteration number.
% Options.PWA_savefilelast='PWA_lastsave'
%       Specifies filename for info about last sucessful iteration.
% Options.PWA_warmstart=0 (Not implemented yet)
%       Loads intermediate result (latest or n), and continues feasibility iteration. 
% Options.PWA_maxTsetTime=Inf
%       Interrupts the feasibility iterations after specified time, and continues with 
%           current feasible set.
% Options.details
%       If set to 1, provides some details about the solution in
%       truct.details. (default value controlled by mptOptions.details)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrlStruct    - Controller structure with following fields:
%   Pn,Fi,Gi    - for region Pn(i).H*x <= Pn.K(i) computed input is U=Fi{i}*x+Gi{i}   
%   Ai,Bi,Ci    - cost associated to each region (x'Aix + Bi*x + Ci)
%                 Note that Ai and Bi are zero matrices, Ci contains the
%                 step distance to the origin
%   Pfinal      - The maximum control invariant set as a polytope or a polyarray
%   dynamics    - Dynamics active in region Pn(i)
%   details     - Details about the solution
%
% see also MPT_CONTROL, MPT_ITERATIVE
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
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

if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end
origProbStruct = probStruct;
fprintf('\n');

if nargin<3,
    Options = [];
end
if ~isfield(Options,'maxiterations'),
    Options.maxiterations=100;
end
if ~isfield(Options,'details'),
    Options.details = mptOptions.details;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'feasset'),
    % if set to 1, function returns maximum controllable set instead of a
    % controller
    Options.feasset=0;
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options, 'lyapunov_type'),
    Options.lyapunov_type = 'any';
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


%%AL added options
if ~isfield(Options,'PWA_savemode'),
    Options.PWA_savemode=0;
end
if ~isfield(Options,'PWA_warmend'),
    Options.PWA_warmend=0;
end
if ~isfield(Options,'PWA_loadn'),
    Options.PWA_loadn=0;
end
if ~isfield(Options,'PWA_maxTsetTime'), 
    Options.PWA_maxTsetTime=Inf;
end
if ~isfield(Options,'PWA_savefile'), 
    Options.PWA_savefile='PWA_save';
end
if ~isfield(Options,'PWA_savefilelast'), 
    Options.PWA_savefilelast='PWA_lastsave';
end
if ~isfield(Options,'iterative'),
    Options.iterative=0;
end
if ~isfield(Options,'oldTset')
    Options.oldTset=0;
end
if ~isfield(Options,'ifa_special'),
    Options.ifa_special = 0;
end

if probStruct.subopt_lev==2,
    Options.finalOneStep = 1;
else
    Options.finalOneStep = 0;
end

if ~iscell(sysStruct.A),
    % LTI system passed, convert it to PWA
    sysStruct = mpt_lti2pwa(sysStruct);
end
[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);

ssInfo.nx = nx;
ssInfo.nu = nu;
ssInfo.ny = ny;
ssInfo.nPWA = nPWA;
ssInfo.nbool = nbool;
ssInfo.ubool = ubool;
ssInfo.intInfo = intInfo;

emptypoly=polytope;

originindynamics=[];
if isfield(probStruct,'xref'),
    origin = probStruct.xref;
else
    origin = zeros(nx,1);
end
ctr=0;
for ii=1:nPWA
    if(all(sysStruct.f{ii}==0)) | isfield(probStruct,'xref') | isfield(probStruct,'uref')
        %.... othewise the origin would not be an equilibrium point
        if all(sysStruct.guardU{ii}==0) & max(sysStruct.guardX{ii}*origin - sysStruct.guardC{ii}) <= Options.abs_tol,
            originindynamics = [originindynamics ii];
        elseif(any(sysStruct.guardU{ii}~=0))
            tempP=polytope(sysStruct.guardU{ii},-sysStruct.guardX{ii}*origin + sysStruct.guardC{ii});
            if(isfulldim(tempP))
                originindynamics = [originindynamics ii];
            end
        end
    end
end

%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%  COMPUTE TARGET SET
%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(~isfield(probStruct,'Tset') | ~isfulldim(probStruct.Tset))
    
    if(isempty(originindynamics))
        error('mpt_iterativePWA: No dynamic has the origin as an equilibrium point! You must define your own target set...');
    end
    if Options.verbose>=1,
        disp(['origin included in: ' num2str(originindynamics)]);
    end
    
    userTset=0; % true if user provided the terminal set, false otherwise
    % compute target set
    [Pn,Fi,Gi,dynamics,probStruct]=mpt_computePWATset(sysStruct,probStruct,originindynamics,Options);
else
    userTset=1; % true if user provided the terminal set, false otherwise
    Pn=probStruct.Tset;
    if Options.verbose > -1,
        fprintf('\n')
        disp('Using User-Defined Target Set...')
        disp('Warning: this may not guarantee stability')
    end
    if ~isfield(probStruct, 'P_N'),
        if Options.verbose > 0,
            fprintf('Warning: setting probStruct.P_N to stage cost.\n');
        end
        if iscell(probStruct.Q),
            probStruct.P_N = probStruct.Q{1};
        else
            probStruct.P_N = probStruct.Q;
        end
    end
    if Options.verbose > -1,
        fprintf('\n')
    end
    dynamics = [];
    for ii=1:length(Pn)
        indyn = [];
        for jj=1:nPWA,
            if isfulldim(Pn(ii) & intInfo.PdynX{jj}),
                indyn = [indyn jj];
            end
        end
        dynamics = [dynamics indyn];
        Fi{ii}=zeros(nu,nx);
        Gi{ii}=zeros(nu,1);
    end 
end

invTset = Pn;

for ii=1:nPWA,
    if ~isempty(intInfo.Xintersect{ii}),
        if length(intInfo.Xintersect{ii})>1,
            disp(sprintf('dynamics %d overlaps with dynamics %s in the X space',ii,mat2str(setdiff(intInfo.Xintersect{ii},ii))));
        end
    end
end

dynamicsorig=dynamics;
[Pu,how]=union(Pn);
TsetMerged=0;
if how,
    % union of target sets is convex, we will merge them
    Pnorig=Pn;
    Fiorig=Fi;
    Giorig=Gi;
    
    clear Pn Fi Gi
    Pn=Pu;
    Fi{1}=Fiorig{1};
    Gi{1}=Giorig{1};
    if Options.verbose>=1,
        disp('mpt_iterativePWA: union of terminal sets is convex, merging...');
    end
    dynamics=1;
    TsetMerged=1;   % set a flag so that we can restore the original solution later
end

nR=length(Pn);
nHulls=nR;

% contains indices of dynamics which do not contain the origin in their interior
originnotindyn = setdiff(1:nPWA,originindynamics);
% if the 'dynamics' array does not contain all the dynamics, we add a dummy controller to be able to initialize
% all the lists
availabledyns = [originindynamics originnotindyn];
for ii=1:nPWA,
    if ~any(dynamics==availabledyns(ii)),
        dynamics = [dynamics availabledyns(ii)];
        Pn = [Pn Pn(1)];
        Fi{length(Fi)+1}=Fi{1};
        Gi{length(Gi)+1}=Gi{1};
    end
end

dynexplore = cell(1,nPWA);
for ii=1:nPWA,
    dynexplore{ii}=[ii setdiff(availabledyns,ii)];    % gives order in which target sets will be explored
    % i.e. dynamics 1 -> target 1, dynamics 1 -> target 2, dynamics 2 -> target 2, dynamics 2 -> target 1
end

starttime=cputime;
convergedN = -1;

CSstorage = cell(1,intInfo.stacks);
maxPfinal = cell(1,intInfo.stacks);
for ctr=1:intInfo.stacks,
    CSstorage{ctr}.Pfinals = emptypoly;
    CSstorage{ctr}.dynamics = [];
    CSstorage{ctr}.Matrices = {};
    CSstorage{ctr}.setdist = [];
    maxPfinal{ctr} = emptypoly;
    maxPfinalStep{1}{ctr} = emptypoly;
end
emptyCS = CSstorage;
for ii=1:length(Pn),
    dyn = dynamics(ii);
    dyn_stack = intInfo.dyns_links(dyn,2);
    CSstorage{dyn_stack}.Pfinals = [CSstorage{dyn_stack}.Pfinals Pn(ii)];
    CSstorage{dyn_stack}.dynamics = [CSstorage{dyn_stack}.dynamics dyn];
    CSstorage{dyn_stack}.Matrices{end+1} = [];
    CSstorage{dyn_stack}.setdist = [CSstorage{dyn_stack}.setdist 0];
end
Step{1} = CSstorage;
for ii=1:intInfo.stacks,
    maxPfinalStep{1}{ii} = CSstorage{ii}.Pfinals;
end

isaddnoise=mpt_isnoise(sysStruct.noise);
Options.noNoiseOnTset=1;
Options.ispwa=1;

targetMatrices = cell(1,nPWA);
for ii=1:nPWA,
    targetMatrices{ii} = {};
end

keepM = [];
Pbnd = sysStruct.Pbnd;
[ObndA, Obndb] = double(Pbnd);

nTargets = 0;
if Options.verbose > -1,
    disp('Generating feasible set');
end
startTsetTime = cputime;

OptionsXU=Options;
OptionsXU.reduce=0;
OptionsXU.constructset=0;
OptionsXU.reduce_output=0;

startTime = cputime;
feasible = 0;

if probStruct.subopt_lev==2,
    Options.iterative=1;
end

if Options.ifa_special,
    mpStorage = {};
    dummyStruct.Pn = invTset;
    dummyStruct.Pfinal = Pn;
    mpStorage{1}{1} = dummyStruct;
    pontryagin_complexity = {};
end

maxxHull = emptypoly;

if statusbar,
    statbar.handle = mpt_statusbar('Computing...');
    statbar.progress = 0;
    Options.verbose = -1;
end

for horizon = 2:Options.maxiterations+1,
    if Options.ifa_special,
        mpStorage{horizon} = {};
    end
    if Options.verbose > -1,
        fprintf('new horizon: %d        \r', horizon-1);
    end
    expand=zeros(1,nPWA);
    maxhull{horizon} = emptypoly;
    
    for ii=1:intInfo.stacks,
        PFstorage{ii}.stack = {};
    end
    
    min_progress = 0;
    max_progress = 0.4;
    if statusbar,
        progress = mod(horizon-2, 50) / 50;
        if isempty(mpt_statusbar(statbar.handle, progress, min_progress, max_progress)),
            mpt_statusbar;
            error('Break...');
        end     
    end
    
    
    if isaddnoise,
        maxSet=minus([maxPfinalStep{horizon-1}{:}],sysStruct.noise,Options);      %take minkowski difference of the target set
        
        if Options.ifa_special,
            pontryagin_complexity{horizon-1}.before = length([maxPfinalStep{horizon-1}{:}]);
            pontryagin_complexity{horizon-1}.after = length(maxSet);
        end
        
        % at first iteration we can afford to compute convex hull and to try to reduce the complexity
        if horizon==2,
            maxSet=union(maxSet,Options);
        end
        if Options.verbose>1,
            disp(sprintf('maxSet consists of %d regions',length(maxSet)));
        end
        if Options.verbose<2,
            disp(sprintf('exploring %d target sets in this iteration',length(maxSet)));
        end
    end
    
    CSstorage = Step{horizon-1};
    newCSstorage = emptyCS;
    
    for source_dyn = 1:nPWA,
        AllTsets = [];
        Pstorage = emptypoly;
        source_dyn_stack = intInfo.dyns_links(source_dyn,2);
        
        if ~Options.iterative | 1,
            tmpProbStruct = probStruct;
            tmpProbStruct.N = 1;
            tmpProbStruct.Tset = emptypoly;
            tmpProbStruct.Tconstraint = 0;
            
            % even if 1/Inf norm is requested, we use 2-norm here
            % to get faster projection. New matrices with 1/Inf norm are
            % constructed afterwards.
            tmpProbStruct.norm=2;
            
            cmOptions = Options;
            cmOptions.pwa_index = source_dyn;
            cmOptions.includeLQRset = 0;
            cmOptions.verbose = 0;
            cmOptions.noReduce = 0;
            BMatrices=mpt_constructMatrices(sysStruct,tmpProbStruct,cmOptions); 
            W = BMatrices.W;
            if isinf(-BMatrices.W),
                % problem is infeasible
                continue
            end
        end
        
        if isaddnoise,
            PF_prev = maxSet;
        else
            PF_prev = emptypoly;
            
            % first extract targets regions from dynamics source_dyn
            if ~isempty(CSstorage{source_dyn_stack}.dynamics),
                source_dyn_ind = find(CSstorage{source_dyn_stack}.dynamics==source_dyn);
            else
                source_dyn_ind = [];
            end
            PF = CSstorage{source_dyn_stack}.Pfinals(source_dyn_ind);
            PF_prev = PF;
            % now add the rest
            other_dyn_ind = setdiff(1:length(CSstorage{source_dyn_stack}.dynamics),source_dyn_ind);
            PF = CSstorage{source_dyn_stack}.Pfinals(other_dyn_ind);
            PF_prev = [PF_prev PF];
            
            for ii=setdiff(1:intInfo.stacks,source_dyn_stack)
                PF_prev = [PF_prev CSstorage{ii}.Pfinals];
            end
            
            % PF_prev now contains all target sets to explore in this iteration
        end
        
        if statusbar,
            progress = mod(horizon-2, 50) / 50;
            if isempty(mpt_statusbar(statbar.handle, progress, 0, 0.4)),
                mpt_statusbar;
                error('Break...');
            end     
        end
        
        for reg=1:length(PF_prev)
            Tset = PF_prev(reg);
            
            [Matrices,mfeas] = mpt_addTset(sysStruct, BMatrices, Tset,nx,nu,source_dyn);
            if ~mfeas,
                continue
            end
            
            feasible = 1;
            W = Matrices.W;
            G = Matrices.G;
            E = Matrices.E;
            expands = 0;
            if isfulldim(maxPfinal{source_dyn_stack}),
                % do the regiondiff in the lifted XU space to find out, if the currently explored transition will expand our maximum feasible set
                PP=regiondiffXU(Pbnd,maxPfinal{source_dyn_stack},G,W,E,OptionsXU);
                if isfulldim(PP),
                    expands=1;
                end
            else
                % this is the first exploration, force exploration
                expands=1;
            end
            if expands,
                if Options.iterative==0,
                    P=polytope([-E G],W,0,2);
                    if isfulldim(P),
                        pOptions = Options;
                        pOptions.noReduce = 1;
                        Pfinal=projection(P,1:size(E,2),pOptions);
                    else
                        Pfinal = emptypoly;
                    end
                else
                    pOptions = Options;
                    pOptions.noReduce = 1;
                    [expanded,Pfinal] = sub_minSwitch(sysStruct, probStruct, Tset, BMatrices, nx,nu,source_dyn,cmOptions, pOptions);
                end
                if isfulldim(Pfinal),
                    if ~isminrep(Pfinal),
                        Pfinal = reduce(Pfinal);
                    end
                    
                    if Options.ifa_special,
                        mpqpOptions = Options;
                        mpqpOptions.verbose = 0;
                        [Pn,Fi,Gi,activeConstraints,PPfinal]=mpt_mpqp(Matrices,mpqpOptions);
                        dummystruct.Pn = Pn;
                        dummystruct.Pfinal = PPfinal;
                        mpStorage{horizon}{end+1} = dummystruct;
                    end
                    
                    if probStruct.norm~=2 & probStruct.subopt_lev~=2,
                        % construct new matrices if 1/Inf norm is
                        % required (remember earlier we've set the norm
                        % to 2 to get faster projection computation
                        tmpProbStruct = probStruct;
                        tmpProbStruct.N = 1;
                        tmpProbStruct.Tset = Tset;
                        tmpProbStruct.Tconstraint = 2;
                        cmOptions = Options;
                        cmOptions.pwa_index = source_dyn;
                        cmOptions.includeLQRset = 0;
                        cmOptions.verbose = 0;
                        cmOptions.noReduce = 0;
                        Matrices=mpt_constructMatrices(sysStruct,tmpProbStruct,cmOptions); 
                    end
                    
                    Matrices.setdist = horizon - 1;
                    Matrices.Tset = Tset;
                    
                    newCSstorage{source_dyn_stack}.dynamics = [newCSstorage{source_dyn_stack}.dynamics source_dyn];
                    newCSstorage{source_dyn_stack}.setdist = [newCSstorage{source_dyn_stack}.setdist horizon-1];
                    newCSstorage{source_dyn_stack}.Matrices{end+1} = Matrices;
                    [newCSstorage{source_dyn_stack}.Pfinals,keep] = sub_mpt_expandlist(Pfinal, newCSstorage{source_dyn_stack}.Pfinals);
                    notkept = find(keep==0);
                    if ~isempty(notkept),
                        kept = setdiff(1:length(keep),notkept);
                        newCSstorage{source_dyn_stack}.dynamics(notkept) = [];
                        newCSstorage{source_dyn_stack}.setdist(notkept) = [];
                        newCSstorage{source_dyn_stack}.Matrices = {newCSstorage{source_dyn_stack}.Matrices{kept}};
                    end
                end
                expand(source_dyn) = 1;
            end % expands 
        end % region
    end %source_dyn
    
    if statusbar,
        progress = mod(horizon-2, 50) / 50;
        if isempty(mpt_statusbar(statbar.handle, progress, 0, 0.4)),
            mpt_statusbar;
            error('Break...');
        end     
    end
    
    for ii=1:intInfo.stacks,
        [newCSstorage{ii}.Pfinals,keep] = reduceunion(newCSstorage{ii}.Pfinals,Options);
        notkept = find(keep==0);
        if ~isempty(notkept),
            kept = setdiff(1:length(keep),notkept);
            newCSstorage{ii}.dynamics(notkept) = [];
            newCSstorage{ii}.setdist(notkept) = [];
            newCSstorage{ii}.Matrices = {newCSstorage{ii}.Matrices{kept}};
        end
        nTargets = nTargets + length(newCSstorage{ii}.setdist);
        for jj=1:length(newCSstorage{ii}.Pfinals),
            maxPfinal{ii} = sub_mpt_expandlist(newCSstorage{ii}.Pfinals(jj),maxPfinal{ii});
        end
        maxPfinal{ii} = reduceunion(maxPfinal{ii},Options);
        
        maxPfinalStep{horizon}{ii} = maxPfinal{ii};
        maxTargetDyn{horizon}{ii} = newCSstorage{ii}.dynamics;
    end
    
    Step{horizon} = newCSstorage;
    
    if all(expand==0),
        break
    end
    
    if statusbar,
        progress = mod(horizon-1, 50) / 50;
        if isempty(mpt_statusbar(statbar.handle, progress, 0, 0.4)),
            mpt_statusbar;
            error('Break...');
        end     
    end
    
end % horizon

if ~feasible,
    if closestatbar,
        mpt_statusbar;
    end
    error('Problem is infeasible!');
end

if horizon>=Options.maxiterations & ~all(expand==0),
    disp('Maximum number of iterations reached without convergence...');
    disp('Please increase the value of Options.maxiterations and run the procedure again.');
    %error('Maximum number of iterations reached!');
    TsetResult='Feasible set search: maxiterations exceeded';
    convergedN = horizon;
end

TsetTime = cputime - startTsetTime;

if horizon>=Options.maxiterations & ~all(expand==0),
    disp('Maximum number of iterations reached without convergence...');
    disp('Please increase the value of Options.maxiterations and run the procedure again.');
    %error('Maximum number of iterations reached!');
    TsetResult='Feasible set search: maxiterations exceeded';
    convergedN = horizon;
end

maxxHull = [maxPfinal{:}];

if Options.feasset,
    if closestatbar,
        mpt_statusbar;
    end
    maxxHull = reduceunion(maxxHull);
    ctrlStruct = maxxHull;
    return
end

if Options.finalOneStep,
    
    min_progress = max_progress;
    max_progress = 0.45;
    mergeOpt = Options;
    mergeOpt.statusbar = 0;
    mergeOpt.verbose = -1;
    maxxHull = reduceunion(maxxHull);
    maxxHull = [invTset maxxHull];
    maxxHull = merge(maxxHull, mergeOpt);
    if Options.verbose > -1,
        disp('Solving one step problem...');
    end
    onestep = {};
    nR = 0;
    for dyn=1:nPWA,
        
        prog_min = (max_progress - min_progress) * (dyn-1) / nPWA + min_progress;
        prog_max = (max_progress - min_progress) * dyn / nPWA + min_progress;

        tmpProbStruct = probStruct;
        tmpProbStruct.N = 1;
        tmpProbStruct.Tset = emptypoly;
        tmpProbStruct.Tconstraint = 2;
        if origProbStruct.norm==2,
            tmpProbStruct.norm = Inf;
        end
        cmOptions = Options;
        cmOptions.pwa_index = dyn;
        cmOptions.includeLQRset = 0;
        cmOptions.verbose = 0;
        cmOptions.noReduce = 0;
        BMatrices=mpt_constructMatrices(sysStruct,tmpProbStruct,cmOptions); 
        if isinf(-BMatrices.W),
            % problem is infeasible
            continue
        end
        
        for reg = 1:length(maxxHull),
            
            if statusbar,
                if mod(reg, 5)==0,
                    if isempty(mpt_statusbar(statbar.handle, (reg-1)/length(maxxHull), prog_min, prog_max)),
                        mpt_statusbar;
                        error('Break...');
                    end     
                end
            end
            
            [Matrices,mfeas] = mpt_addTset(sysStruct, BMatrices, maxxHull(reg),nx,nu,dyn);
            if ~mfeas,
                continue
            end
            
            mplpOptions = Options;
            mplpOptions.verbose = 0;
            mplpOptions.nu = size(sysStruct.B{1},2);
            [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
            if ~isfulldim(Pn),   % mplp is numerically VERY sensitive! in case it did not find a solution, we recompute it, since it has to exist!
                Fi={}; Gi={}; details.Fi={}; details.Gi={};
                [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
            end
            if ~isempty(Fi) & isfulldim(Pn(1)),
                cs = {};
                cs.Pfinal = Pfinal;
                cs.sysStruct = sysStruct;
                cs.probStruct = tmpProbStruct;
                cs.overlaps = 0;
                cs.Pn = Pn;                       % transitional controller partition
                cs.Fi = Fi;                       % control law
                cs.Gi = Gi;                       % control law
                for dummy=1:length(Fi),
                    cs.Ai{dummy} = zeros(nx);               % value function, quadratic term
                end
                cs.Bi = details.Bi;               % value function, linear term
                cs.Ci = details.Ci;               % value function, constant term
                cs.dynamics = dyn*ones(1,length(Fi)); % dynamics
                cs.details = details;
                nR = nR + length(cs.Bi);
                onestep{end+1} = cs;
            end
        end % maxxHull
    end % nPWA
    %%nR
    
    min_progress = max_progress;
    max_progress = 0.6;

    roOptions = Options;
    if statusbar,
        roOptions.statusbar = 1;
        roOptions.status_min = min_progress;
        roOptions.status_max = max_progress;
        roOptions.status_handle = statbar.handle;
        roOptions.closestatbar = 0;
    end

    ctrlStruct = mpt_removeOverlaps(onestep, roOptions);
    ctrlStruct.probStruct=probStruct;
    ctrlStruct.probStruct.subopt_lev=0;
    ctrlStruct = rmfield(ctrlStruct,'details');
    ctrlStruct.details.runTime = cputime - startTime;
    ctrlStruct.details.TsetTime = TsetTime;
    ctrlStruct.details.PWATime=0;
    ctrlStruct.details.PWQTime=0;
    if Options.verbose > -1,
        fprintf('%d regions generated\n',length(ctrlStruct.Pn));
    end
    
    min_progress = max_progress;
    max_progress = 1;
    lyapOptions = Options;
    if statusbar,
        lyapOptions.statusbar = 1;
        lyapOptions.status_min = min_progress;
        lyapOptions.status_max = max_progress;
        lyapOptions.status_handle = statbar.handle;
        lyapOptions.closestatbar = 0;
    end
        

    if strcmpi(Options.lyapunov_type, 'none'),
        % don't compute any lyapunov function
        feasible = 0;
    elseif mpt_isnoise(sysStruct.noise)
        % system is subject to additive noise, can compute only common lyapunov
        % function
        disp('System is subject to additive noise, calculating Quadratic Lyapunov function...');
        try
            [lyapunovP, drho, feasible] = mpt_getQuadLyapFct(ctrlStruct, lyapOptions);
        catch
            feasible = 0;
        end
        if feasible
            ctrlStruct.details.lyapP = lyapunovP;
        end
        
    else
        % compute PWA or PWQ Lyapunov function to guarantee stability
        feasible_pwq = 0;
        feasible_pwa = 0;
        
        if strcmpi(Options.lyapunov_type, 'pwq') | strcmpi(Options.lyapunov_type, 'any'),
            if Options.verbose > -1,
                disp('Computing PWQ Lyapunov function...');
            end
            feasible_pwq=0;
            try
                startt = cputime;
                [lQ, lL, lC, feasible_pwq] = mpt_getPWQLyapFct(ctrlStruct, lyapOptions);
                ctrlStruct.details.PWQTime = cputime-startt;
            catch
                feasible_pwq=0;
            end
            if(feasible_pwq)
                ctrlStruct.details.PWQlyapQ=1;
                if Options.details,
                    ctrlStruct.details.lyapPWQ_Q = lQ;
                    ctrlStruct.details.lyapPWQ_L = lL;
                    ctrlStruct.details.lyapPWQ_C = lC;
                end
            else
                if Options.verbose > 0,
                    fprintf('PWQ Lyapunov function not found!\n');
                end
                ctrlStruct.details.PWQlyapQ=0;
            end
        end
        if (~feasible_pwq & strcmpi(Options.lyapunov_type, 'any')) | ...
                strcmpi(Options.lyapunov_type, 'pwa'),
            if Options.verbose > -1,
                fprintf('Computing PWA Lyapunov function...\n');
            end
            try
                startt = cputime;
                [lL, lC, feasible_pwa] = mpt_getPWALyapFct(ctrlStruct, lyapOptions);
                ctrlStruct.details.PWATime = cputime-startt;
            catch
                feasible_pwa = 0;
            end
            if feasible_pwa,
                ctrlStruct.details.PWAlyapL=1;
                if Options.details,
                    ctrlStruct.details.lyapPWA_L = lL;
                    ctrlStruct.details.lyapPWA_C = lC;
                end
            else
                if Options.verbose > 0,
                    fprintf('PWA Lyapunov function not found!\n');
                end                
                ctrlStruct.details.PWAlyapL=0;    
            end
        end
        feasible = feasible_pwq | feasible_pwa;
    end
    
    if Options.verbose > -1,
        if strcmpi(Options.lyapunov_type, 'none'),
            fprintf('\nNo Lyapunov function was computed since Options.lyapunov_type=''none''.\n');
            fprintf('Note that the closed-loop system may be unstable!\n\n');
        elseif feasible,
            fprintf('\nLyapunov function found, system is stable\n\n');
        else
            fprintf('\nLyapunov function NOT found, stability is not guaranteed!\n\n');
        end     
    end
    
    if closestatbar,
        mpt_statusbar;
    end
    ctrlStruct.probStruct=origProbStruct;
    return
end

if Options.verbose > -1,
    disp(sprintf('%d hulls to be explored',nTargets));
end

%% AL
if Options.PWA_warmend | Options.PWA_savemode
    % Save results thus far
    disp('Saving all hulls..')
    eval(['save ' Options.PWA_savefile 'hulls' num2str(horizon)])
    disp('...done')
end
% Should also add possibility to jump here and load an earlier result file....

if Options.verbose > -1,
    disp('Generating regions');
end

startMPtime = cputime;
nR = 0;
PP = [];
Pdyn = cell(1,nPWA);

if ~TsetMerged,
    Pnorig = Pn;
    Fiorig = Fi;
    Giorig = Gi;
end

targetMatrices = cell(1,nPWA);
for ii=1:nPWA,
    targetMatrices{ii} = {};
end
for ii=1:length(Step),
    for jj=1:length(Step{ii}),
        dyns = Step{ii}{jj}.dynamics;
        for kk=1:length(dyns),
            M = Step{ii}{jj}.Matrices{kk};
            if ~isempty(M),
                targetMatrices{dyns(kk)}{end+1} = M;
            end
        end
    end
end

if isfield(Options, 'projection'),
    Options = rmfield(Options, 'projection');
end

nRtot = 0;
pwalist = cell(1,nPWA);
min_progress = max_progress;
max_progress = 0.8;
for dyn=1:nPWA,
    
    prog_min = (max_progress - min_progress) * (dyn-1) / nPWA + min_progress;
    prog_max = (max_progress - min_progress) * dyn / nPWA + min_progress;
    
    pwalist{dyn} = {};
    PP = [];
    targets = targetMatrices{dyn};
    lenT = length(targets);
    for ii=1:lenT,
        
        if statusbar,
            if mod(ii, 3)==0,
                if isempty(mpt_statusbar(statbar.handle, (ii-1)/lenT, prog_min, prog_max)),
                    mpt_statusbar;
                    error('Break...');
                end     
            end
        end
        
        if ii==1 | ii==lenT | mod(ii,round(lenT/3))==0,
            if Options.verbose > -1,
                fprintf('dynamics %d, hull %d/%d        \r', dyn, ii, lenT);
            end
        end
        Matrices = targets{ii};
        targets{ii} = {}; % we don't need this element any more, erase it to save memory
        if isempty(Matrices),
            continue
        end
        if Options.iterative
            tmpProbStruct = probStruct;
            tmpProbStruct.N = 1;
            tmpProbStruct.Tset = Matrices.Tset;
            tmpProbStruct.Tconstraint = 2;
            locOpt = Options;
            locOpt.ispwa = 1;
            locOpt.pwa_index = dyn;
            locOpt.verbose = 0;
            locOpt.ssInfo = ssInfo;
            [ccs] = mpt_iterative(sysStruct, tmpProbStruct, locOpt);
            Pn = fliplr(ccs.Pn);
            Fi = fliplr(ccs.Fi);
            Gi = fliplr(ccs.Gi);
            nR = length(Fi);
            details.Ai=cell(1,nR);    % value function
            details.Bi=cell(1,nR);
            details.Ci=cell(1,nR);
            for region=1:nR;
                details.Ai{region} = zeros(nx);
                details.Bi{region} = zeros(1,nx);
                details.Ci{region} = Matrices.setdist;             % cost is simply the step-distance to the invariant set
            end
        else
            if probStruct.norm==2,
                mpqpOptions = Options;
                mpqpOptions.verbose = 0;
                [Pn,Fi,Gi,activeConstraints,Pfinal]=mpt_mpqp(Matrices,mpqpOptions);
                nR = length(Fi);
                details.Ai=cell(1,nR);    % value function
                details.Bi=cell(1,nR);
                details.Ci=cell(1,nR);
                for region=1:nR;
                    details.Ai{region} = zeros(nx);
                    details.Bi{region} = zeros(1,nx);
                    details.Ci{region} = Matrices.setdist;             % cost is simply the step-distance to the invariant set
                end
            else
                mplpOptions = Options;
                mplpOptions.verbose = 0;
                mplpOptions.nu = size(sysStruct.B{1},2);
                [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
                if ~isfulldim(Pn),   % mplp is numerically VERY sensitive! in case it did not find a solution, we recompute it, since it has to exist!
                    [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
                end
                nR = length(Fi);
                details.Ai=cell(1,nR);
                details.Bi=cell(1,nR);
                details.Ci=cell(1,nR);
                nx=dimension(Pn);
                for region=1:nR,
                    details.Ai{region}=zeros(nx);
                    details.Bi{region}=zeros(1,nx);
                    details.Ci{region}=Matrices.setdist;     % cost is simply the step-distance to the invariant set
                end
            end
        end
        if ~isempty(Fi) & isfulldim(Pn(1)),
            list = {};
            list.Pfinal = Pfinal;
            list.Pn = Pn;                       % transitional controller partition
            list.Fi = Fi;                       % control law
            list.Gi = Gi;                       % control law
            list.Ai = details.Ai;               % value function, quadratic term
            list.Bi = details.Bi;               % value function, linear term
            list.Ci = details.Ci;               % value function, constant term
            list.setdistance=Matrices.setdist;
            pwalist{dyn}{end+1} = list;
            nRtot = nRtot + length(list.Pn);
        end
    end
end
mpTime = cputime - startMPtime;

%%AL
if Options.PWA_warmend | Options.PWA_savemode
    % Save results thus far. 
    disp('Saving all regions..')
    eval(['save ' Options.PWA_savefile 'regions' num2str(horizon)])
    disp('...done')
end
% Should also add possibility to jump here and load an earlier result file....
if Options.verbose > -1,
    fprintf('\n%d regions generated\n',nRtot);
end

startRedTime = cputime;
if Options.verbose > -1,
    disp('Removing unused regions (backwards)...');
end
cs = {};
min_progress = 0.8;
max_progress = 1;
for dyn=1:nPWA,
    
    prog_min = (max_progress - min_progress) * (dyn-1) / nPWA + min_progress;
    prog_max = (max_progress - min_progress) * dyn / nPWA + min_progress;
    
    stack_pos = intInfo.dyns_links(dyn,2);
    sPn = [];
    sFi = {};
    sGi = {};
    sAi = {};
    sBi = {};
    sCi = {};
    list = pwalist{dyn};
    clear newlist
    lenL=length(list);
    for ii=lenL:-1:1,
        %%AL
        
        if statusbar,
            if mod(ii, 3)==0,
                if isempty(mpt_statusbar(statbar.handle, (lenL - ii)/lenL, prog_min, prog_max)),
                    mpt_statusbar;
                    error('Break...');
                end     
            end
        end
        
        if ii==1 | ii==lenL | mod(ii,round(lenL/3))==0,
            if Options.verbose>0,
                if Options.verbose > -1,
                    fprintf('dynamics %d, hull %d/%d           \r', dyn, ii, lenL);
                end
            end
        end
        
        Pn = list{ii}.Pn;
        nR = length(list{ii}.Fi);
        setdist = list{ii}.setdistance;
        for jj = 1:nR
            if setdist==1 & userTset,
                sPn = [sPn Pn(jj)];
                sFi{end+1} = list{ii}.Fi{jj};
                sGi{end+1} = list{ii}.Gi{jj};
                sAi{end+1} = list{ii}.Ai{jj};
                sBi{end+1} = list{ii}.Bi{jj};
                sCi{end+1} = list{ii}.Ci{jj};
                continue
            end
            R = regiondiff(Pn(jj),maxPfinalStep{setdist}{stack_pos},OptionsXU);
            if isfulldim(R),
                sPn = [sPn Pn(jj)];
                sFi{end+1} = list{ii}.Fi{jj};
                sGi{end+1} = list{ii}.Gi{jj};
                sAi{end+1} = list{ii}.Ai{jj};
                sBi{end+1} = list{ii}.Bi{jj};
                sCi{end+1} = list{ii}.Ci{jj};
            end
        end
    end
    newlist.Pn = sPn;
    newlist.Fi = sFi;
    newlist.Gi = sGi;
    newlist.Ai = sAi;
    newlist.Bi = sBi;
    newlist.Ci = sCi;
    cs{dyn} = newlist;
end
fprintf('\n');

if statusbar,
    if isempty(mpt_statusbar(statbar.handle, 1)),
        mpt_statusbar;
        error('Break...');
    end     
end

if ~userTset,
    % if the initial target sets were merged at the beginning, we put them back to particular list
    %for ii=1:length(dynamicsorig),
    for ii=1:length(Pnorig),
        actdyn = 0;
        for jj=1:nPWA,
            Pdyn = intInfo.Pdyn{jj};
            if isfulldim(Pnorig(ii) & Pdyn),
                actdyn = jj;
                break
            end
        end
        if actdyn==0,
            if statusbar,
                mpt_statusbar;
            end
            error('mpt_iterativePWA: no dynamics associated to target set! Contact mpt@control.ee.ethz.ch if this happens.');
        end
        cs{actdyn}.Pn=[cs{actdyn}.Pn Pnorig(ii)];
        cs{actdyn}.Fi{end+1}=Fiorig{ii};
        cs{actdyn}.Gi{end+1}=Giorig{ii};
        cs{actdyn}.Ai{end+1}=zeros(nx);
        cs{actdyn}.Bi{end+1}=zeros(1,nx);
        cs{actdyn}.Ci{end+1}=0;
    end
end

Pn = emptypoly;
nRtot = 0;
for dyn=1:nPWA,
    nRtot = nRtot + length(cs{dyn}.Fi);
end
Fi = cell(1,nRtot);
Gi = cell(1,nRtot);
Ai = cell(1,nRtot);
Bi = cell(1,nRtot);
Ci = cell(1,nRtot);
dynamics = zeros(1,nRtot);
count = 0;
for ii=0:length(Step)-1,
    for dyn=1:nPWA,
        for jj=length(cs{dyn}.Fi):-1:1,
            if cs{dyn}.Ci{jj}==ii,
                count = count+1;
                Pn = [Pn cs{dyn}.Pn(jj)];
                Fi{count} = cs{dyn}.Fi{jj};
                Gi{count} = cs{dyn}.Gi{jj};
                Ai{count} = cs{dyn}.Ai{jj};
                Bi{count} = cs{dyn}.Bi{jj};
                Ci{count} = cs{dyn}.Ci{jj};
                dynamics(count) = dyn;
            end
        end
    end
end
redTime = cputime - startRedTime;

ctrlStruct.sysStruct = sysStruct;
ctrlStruct.probStruct = origProbStruct;
ctrlStruct.Pn = Pn;
ctrlStruct.Fi = Fi;
ctrlStruct.Gi = Gi;
ctrlStruct.Ai = Ai;
ctrlStruct.Bi = Bi;
ctrlStruct.Ci = Ci;
ctrlStruct.Pfinal = emptypoly;
ctrlStruct.Pfinal = [ctrlStruct.Pfinal maxPfinal{:}];
ctrlStruct.dynamics = dynamics;
ctrlStruct.overlaps = 1;
ctrlStruct.details.runTime = cputime - startTime;
ctrlStruct.details.TsetTime = TsetTime;
ctrlStruct.details.mpTime = mpTime;
ctrlStruct.details.reductionTime = redTime;
ctrlStruct.details.N = length(Step)-1;
ctrlStruct.details.invTset = invTset;
if Options.details,
    ctrlStruct.details.maxHull = maxxHull;
    ctrlStruct.details.maxHullStep = maxPfinalStep;
end

if Options.ifa_special,
    ctrlStruct.details.mpStorage = mpStorage;
    ctrlStruct.details.pontryagin_complexity = pontryagin_complexity;
end

if closestatbar,
    mpt_statusbar;
end

return
% and we are done!


% ===========================================================================================
function [expanded,Pfinal] = sub_minSwitch(sysStruct, probStruct, Tset, BMatrices, nx,nu,source_dyn,cmOptions, pOptions)

Matrices = {};

notconverged = 1;
expanded=0;
emptypoly = polytope;
Pfinal = emptypoly;

while notconverged,
    [M,mfeas] = mpt_addTset(sysStruct, BMatrices, Tset,nx,nu,source_dyn);
    if ~mfeas,
        return
    end
    
    W = M.W;
    if isinf(-W),
        % problem is infeasible
        return
    end
    G = M.G;
    E = M.E;
    P=polytope([-E G],W,0,2);
    if isfulldim(P),
        Pfinal=projection(P,1:size(E,2),pOptions);
    else
        break
    end
    if isfulldim(Pfinal),
        if ~isminrep(Pfinal),
            Pfinal = reduce(Pfinal);
        end
        if Pfinal > Tset,
            expanded=1;
            Tset = Pfinal;
        else
            break
        end
    else
        break
    end
    notconverged=expanded;
end



% ===========================================================================================
function [Pun,keep] = sub_mpt_expandlist(Pf,Pu,Options),
% given a polytope Pf and a polyarray Pu, using a subset check removes all fields of Pu
% which are covered by Pf

Options.reduce=0;        % for fast subset check
Options.constructset=0;  % for fast subset check
Options.elementwise=1;   % go elementwise, i.e. Pu(1)<=Pf, Pu(2)<=Pf, ... , Pu(n)<=Pf
expands = 1;
if ~isfulldim(Pu(1)),
    % if Pu is empty, returns Pf
    Pun=Pf;
    keep=1;
    return
end

PuExt=[Pu Pf];
lenPuExt = length(PuExt);
keep = [];
if lenPuExt>1,
    keep=(~le(Pu,Pf,Options))'; % returns indices of polytopes in Pu which are not a subset of Pf
    keep=[keep 1];
else
    keep=1;
end

Pun=PuExt(find(keep~=0));    % the output consists of polytopes which are not a subset of Pf
return

% ===========================================================================================
function [Pun,keep] = sub_mpt_expandlist2(Pf,Pu,Options)
% given a polytope Pf and a polyarray Pu, using regiondiff removes all fields of Pu
% which are completely covered by Pf

Options.reduce=0;         % for fast subset check
Options.constructset=0;   % for fast subset check
Options.reduce_output=0;  % for fast subset check

if ~isfulldim(Pu(1)),
    % if Pu is empty, returns Pf
    Pun=Pf;
    return
end

Pu=[Pu Pf];
keep=ones(1,length(Pu));
emptypoly=polytope;
lenPu=length(Pu);

for ii=1:lenPu-1,
    Pun=emptypoly;
    Pun=Pu([1:ii-1,ii+1:lenPu]);         % all indices except of the current one
    P=regiondiff(Pu(ii),Pun,Options);    % if output of regiondiff is empty, Pu(ii) is NOT covered by Pu(I-ii)
    if ~isfulldim(P(1)),
        keep(ii)=0;
    end
end

Pun=Pu(find(keep~=0));
return

% ===========================================================================================
function [Pun,keep] = sub_mpt_reduce_expandlist2(Pu,Options)
% given a (possibly non-convex) union of polytopes Pu, using regiondiff removes
% all elements of Pu which are completely covered by the other regions
%
% i.e. assume we have 3 polytopes:
%      p1 = polytope([0 -1; -1 0; 0 1; 1 0],[1;0;1;1])
%      p2 = polytope([0 -1; -1 0; 0 1; 1 0],[1;1;1;0])
%      p3 = polytope([0 -1; -1 0; 0 1; 1 0],[1;1;1;1])
%
% then if Pu=[p1 p2 p3], this function removes polytopes p1 and p2, since they are completely
% covered by a larger polytope (p3 in this case)
keep = 1;
Options.reduce=0;         % for fast subset check
Options.constructset=0;   % for fast subset check
Options.reduce_output=0;  % for fast subset check



lenPu=length(Pu);
if lenPu==1,
    % nothing to reduce if just one element in Pu
    Pun=Pu;
    return
end

%sort polyarray according to chebychev radius
polyradius = [];
for i=1:length(Pu)
    [xc,rc]=chebyball(Pu(i));
    polyradius = [polyradius; rc];
end
[polyradius,index] = sort(polyradius);   %sort in ascending order

keep=ones(1,lenPu);     % initialize the vector, 1 means keep the region, 0 kick it out

for ii=1:lenPu,
    ind1=[1:ii-1,ii+1:lenPu];     % all indices except of the current one
    ind2=find(keep==0);           % indices of regions which were kicked out
    ind=setdiff(ind1,ind2);       % remove indices of regions, which were already kicked out
    if isempty(ind),              % if the set of indices is empty, that means that we kicked out all polytope except of the last one
        Pun=Pu(ii);               % so return it
        return
    end
    Pum=Pu(ind);                         % polyarray which could be possibly fully covering Pu(ii)
    
    keepthis = 1;
    for jj=1:length(Pum),
        index_le = le(Pu(ii),Pum(jj),struct('elementwise',1));
        if index_le,
            keepthis = 0;
            break
        end
    end
    %%index_le = le(Pu(ii),Pum,struct('elementwise',1));
    %%keepthis = all(index_le==0);
    if ~keepthis,
        keep(ii)=0;                      % if the solution is not empty, means that Pu(ii) is fully covered by the remaining polytopes, we kick it out
    end
end

Pun=Pu(find(keep==1));           % the result is then a subset of Pu where indices are not zero

return
