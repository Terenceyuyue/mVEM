function ctrlStruct=mpt_mixedMinTime(sysStruct,probStruct,Options)
%MPT_MIXEDMINTIME Computes minimum time controller for systems with discrete and continuous inputs
%
% ctrlStruct = mpt_mixedMinTimePWA(sysStruct,probStruct)
% ctrlStruct = mpt_mixedMinTimePWA(sysStruct,probStruct,Options)
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
% Inputs may be discrete or continuous, as defined in sysStruct.Uset
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
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrlStruct    - Controller structure
%
% see also MPT_CONTROL, MPT_ITERATIVE
%

% Copyright is with the following author(s):
%
% (C) 2004-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch

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
if ~isfield(Options,'bothLyapFct'),
    Options.bothLyapFct=0;
end
if ~isfield(Options,'oldTset')
    Options.oldTset=0;
end

if probStruct.subopt_lev==2,
    Options.finalOneStep = 1;
else
    Options.finalOneStep = 0;
end

if probStruct.tracking==1
    error('Tracking with deltaU formulation is not allowed! Please use "probStruct.tracking=2".');
end
if isfield(sysStruct, 'dumode'),
    error('Sorry, deltaU constraints are not supported!');
end
if isfield(probStruct,'Rdu')
    if ~all(probStruct.Rdu==0),
        error('mpt_optMixedCtrl: non-zero probStruct.Rdu is not supported in this version!');
    end
end

if ~iscell(sysStruct.A),
    % LTI system passed, convert it to PWA
    sysStruct = mpt_lti2pwa(sysStruct);
end
[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);

emptypoly=polytope;

% % originindynamics=[];
% % if isfield(probStruct,'xref'),
% %     origin = probStruct.xref;
% % else
% %     origin = zeros(nx,1);
% % end
% % ctr=0;
% % for ii=1:nPWA
% %     if(all(sysStruct.f{ii}==0)) | isfield(probStruct,'xref') | isfield(probStruct,'uref')
% %         %.... othewise the origin would not be an equilibrium point
% %         if all(sysStruct.guardU{ii}==0) & max(sysStruct.guardX{ii}*origin - sysStruct.guardC{ii}) <= Options.abs_tol,
% %             originindynamics = [originindynamics ii];
% %         elseif(any(sysStruct.guardU{ii}~=0))
% %             tempP=polytope(sysStruct.guardU{ii},-sysStruct.guardX{ii}*origin + sysStruct.guardC{ii});
% %             if(isfulldim(tempP))
% %                 originindynamics = [originindynamics ii];
% %             end
% %         end
% %     end
% % end
% % 
% % if(isempty(originindynamics))
% %     error('mpt_iterativePWA: No dynamic has the origin as an equilibrium point !! Aborting computation...');
% % end
% % if Options.verbose>=1,
% %     disp(['origin included in: ' num2str(originindynamics)]);
% % end
% % 
%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%%  COMPUTE TARGET SET
%%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(~isfield(probStruct,'Tset') | ~isfulldim(probStruct.Tset))
    error('Terminal set must be provided in probStruct.Tset');
    userTset=0; % true if user provided the terminal set, false otherwise
    % compute target set
    [Pn,Fi,Gi,dynamics,probStruct]=mpt_computePWATset(sysStruct,probStruct,originindynamics,Options);
else
    userTset=1; % true if user provided the terminal set, false otherwise
    Pn=probStruct.Tset;
    fprintf('\n')
    disp('Using User-Defined Target Set...')
    disp('WARNING: this may not guarantee stability')
    if ~isfield(probStruct, 'P_N'),
        error('Penalty on final state "probStruct.P_N" must be given!');
    end
    fprintf('\n')
    for ii=1:length(Pn)
        %%dynamics(ii) = originindynamics(1);
        Fi{ii}=zeros(1,nx);
        Gi{ii}=zeros(nu,1);
    end 
end

if any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax)),
    dUconstraints = 1;
else
    dUconstraints = 0;
end
if dUconstraints,
    error('mpt_mixedMinTime: Constraints on delta U not supported in this version!');
end
if isfield(probStruct,'Rdu')
    if ~all(probStruct.Rdu==0),
        error('mpt_mixedMinTime: non-zero probStruct.Rdu is not supported in this version!');
    end
end

invTset = Pn;

for ii=1:nPWA,
    if ~isempty(intInfo.Xintersect{ii}),
        if length(intInfo.Xintersect{ii})>1,
            %disp(sprintf('dynamics %d overlaps with dynamics %s in the X space',ii,mat2str(intInfo.Xintersect{ii})));
            disp(sprintf('dynamics %d overlaps with dynamics %s in the X space',ii,mat2str(setdiff(intInfo.Xintersect{ii},ii))));
        end
    end
end

if ~exist('dynamics','var'),
    dynamics = 1:nPWA;
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

originindynamics=[];
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
    CSstorage{ctr}.Ubool = {};
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
    CSstorage{dyn_stack}.Ubool{end+1} = zeros(nbool,1);
    % %     if ~userTset,
    % %         maxPfinal{dyn_stack} = [maxPfinal{dyn_stack} Pn(ii)];
    % %     end
end
Step{1} = CSstorage;
for ii=1:intInfo.stacks,
    maxPfinalStep{1}{ii} = CSstorage{ii}.Pfinals;
end

Mstorage = cell(1,nPWA);
for dyn=1:nPWA,
    Mstorage{dyn} = {};
end

orig_gU = sysStruct.guardU;
orig_gC = sysStruct.guardC;

% determine parts of state-space where dynamics are active
Pdynamics = intInfo.Pdyn;
[bndA, bndb] = double(sysStruct.Pbnd);
pureX_gX = {}; pureX_gU = {}; pureX_gC = {};
pureU_gX = {}; pureU_gU = {}; pureU_gC = {};
for ii=1:nPWA,
    gX = sysStruct.guardX{ii};
    gC = sysStruct.guardC{ii};
    gU = sysStruct.guardU{ii};
    zerorows = [];
    for jj=1:size(gX,1)
        if all(gX(jj,:)==0),
            zerorows = [zerorows; jj];
        end
    end
    nonzerorows = setdiff(1:size(gX,1),zerorows);
    gXnz = gX(nonzerorows,:);
    gCnz = gC(nonzerorows,:);
    gUnz = gU(nonzerorows,:);
    
    gXz = gX(zerorows,:);
    gUz = gU(zerorows,:);
    gCz = gC(zerorows,:);
    pureX_gX{ii} = gXnz;
    pureX_gU{ii} = gUnz;
    pureX_gC{ii} = gCnz;
    pureU_gX{ii} = gXz;
    pureU_gU{ii} = gUz;
    pureU_gC{ii} = gCz;
end

Udiscr = {};
for ii=1:nbool,
    Udiscr{ii} = sysStruct.Uset{ubool(ii)};
end
Ucombs = mpt_allcombs(Udiscr);

% indices of continuous inputs;
ucont = setdiff(1:nu, ubool);
ncont = length(ucont);
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
disp('Generating feasible set');
startTsetTime = cputime;

OptionsXU=Options;
OptionsXU.reduce=0;
OptionsXU.constructset=0;
OptionsXU.reduce_output=0;
startTime = cputime;
feasible = 0;

if any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax)),
    dUconstraints = 1;
else
    dUconstraints = 0;
end
dUconstraints = 0;
dUctr = 0;

maxxHull = emptypoly;
for horizon = 2:Options.maxiterations+1,
    fprintf('new horizon: %d        \r', horizon-1);
    expand=zeros(1,nPWA);
    maxhull{horizon} = emptypoly;
    
    for ii=1:intInfo.stacks,
        PFstorage{ii}.stack = {};
    end

    CSstorage = Step{horizon-1};
    newCSstorage = emptyCS;

    for possibleU = 1:size(Ucombs,1),
        U = Ucombs(possibleU,:)';
        DDyn=[];
        for ii=1:nPWA,
            % obtain index of dynamics associated to a given U
            if isempty(pureU_gU{ii}),
                DDyn=[DDyn ii];
            elseif all(pureU_gU{ii}(:,ubool)*U - pureU_gC{ii}<=Options.abs_tol),
                DDyn = [DDyn ii];
            end
        end
        if isfield(sysStruct,'Uexcl'),
            if all(ismember(U(:)',sysStruct.Uexcl,'rows')),
                % exclude this U
                disp(['Excluding input ' mat2str(U')]);
                continue
            end
        end
        if isempty(DDyn),
            disp(['No guardline associated to input ' mat2str(U')]);
            continue
        end
        
        if isaddnoise,
            maxSet=minus([maxPfinalStep{horizon-1}{:}],sysStruct.noise,Options);      %take minkowski difference of the target set
            
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
        
        
        for source_dyn = DDyn,
            fprintf('U=%s, dynamics=%d             \r', mat2str(U'),source_dyn);
            
            AllTsets = [];
            Pstorage = emptypoly;
            source_dyn_stack = intInfo.dyns_links(source_dyn,2);
            
            % construct matrices for 1-step poblem
            % we insert target set later
            tmpProbStruct = probStruct;
            tmpProbStruct.Tset = emptypoly;
            tmpProbStruct.Tconstraint=0;
            tmpProbStruct.N=1;
            tmpProbStruct.norm=2;
            tmpProbStruct.R = probStruct.R(ucont,ucont);
            if isfield(tmpProbStruct,'uref'),
                tmpProbStruct.uref = probStruct.uref(ucont);
            end
            
            tmpSysStruct = sysStruct;
            tmpSysStruct.B{source_dyn} = sysStruct.B{source_dyn}(:,ucont);
            tmpSysStruct.f{source_dyn} = sysStruct.f{source_dyn} + sysStruct.B{source_dyn}(:,ubool)*U;
            tmpSysStruct.D{source_dyn} = sysStruct.D{source_dyn}(:,ucont);
            tmpSysStruct.g{source_dyn} = sysStruct.g{source_dyn} + sysStruct.D{source_dyn}(:,ubool)*U;
            
            % eliminate equality constraints from guards on U
            tmpSysStruct.guardC{source_dyn} = pureX_gC{source_dyn};
            tmpSysStruct.guardU{source_dyn} = pureX_gU{source_dyn}(:,ucont);
            tmpSysStruct.guardX{source_dyn} = pureX_gX{source_dyn};
            tmpSysStruct.umax = sysStruct.umax(ucont);
            tmpSysStruct.umin = sysStruct.umin(ucont);
            tmpSysStruct.dumax = sysStruct.dumax(ucont);
            tmpSysStruct.dumin = sysStruct.dumin(ucont);
            
            localoptions=Options;
            localoptions.pwa_index=source_dyn;
            localoptions.verbose=0;
            
            [BMatrices]=mpt_constructMatrices(tmpSysStruct,tmpProbStruct,localoptions);
            if isinf(-BMatrices.W), 
                % if transition infeasible, continue with next target
                continue
            end
            
            if isaddnoise,
                PF_prev = maxSet;
            else
                PF_prev = emptypoly;
                Ubool_prev = {};
                
                % first extract targets regions from dynamics source_dyn
                
                source_dyn_ind = find(CSstorage{source_dyn_stack}.dynamics==source_dyn);
                PF = CSstorage{source_dyn_stack}.Pfinals(source_dyn_ind);
                Ubool_prev = {CSstorage{source_dyn_stack}.Ubool{source_dyn_ind}};
                PF_prev = PF;
                % now add the rest
                other_dyn_ind = setdiff(1:length(CSstorage{source_dyn_stack}.dynamics),source_dyn_ind);
                PF = CSstorage{source_dyn_stack}.Pfinals(other_dyn_ind);
                UB = {CSstorage{source_dyn_stack}.Ubool{other_dyn_ind}};
                PF_prev = [PF_prev PF];
                Ubool_prev = {Ubool_prev{:},UB{:}};
                
                for ii=setdiff(1:intInfo.stacks,source_dyn_stack)
                    PF_prev = [PF_prev CSstorage{ii}.Pfinals];
                    Ubool_prev = {Ubool_prev{:},CSstorage{ii}.Ubool{:}};
                end
                
                % PF_prev now contains all target sets to explore in this iteration
            end
            
            % %             if ~isfulldim(PF_prev)
            % %                 % no more target sets to explore
            % %                 expand(source_dyn)=0;
            % %                 %%continue
            % %             end
            
            for reg=1:length(PF_prev)
                Tset = PF_prev(reg);
                
                deltaU = 0;
                if horizon>2,
                    % control move valid in target region
                    targetU = Ubool_prev{reg};
                    deltaU = U - targetU;
                    if any(deltaU<sysStruct.dumin(ubool)) | any(deltaU>sysStruct.dumax(ubool)),
                        % new control action violates slew constraints, continue
                        % with next sequence of inputs
                        continue
                    end
                end
                
                [Matrices,mfeas] = mpt_addTset(tmpSysStruct, BMatrices, Tset,nx,ncont,source_dyn);
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
                    Pfinal = Pbnd;
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
                        cmOptions = Options;
                        cmOptions.pwa_index = source_dyn;
                        cmOptions.includeLQRset = 0;
                        cmOptions.verbose = 0;
                        cmOptions.noReduce = 0;
                        
                        [expanded,Pfinal] = sub_minSwitch(tmpSysStruct, tmpProbStruct, Tset, BMatrices, nx,ncont,source_dyn,cmOptions, pOptions);
                    end
                    if isfulldim(Pfinal),
                        if ~isminrep(Pfinal),
                            Pfinal = reduce(Pfinal);
                        end
                        
                        if probStruct.norm~=2,
                            % construct new matrices if 1/Inf norm is
                            % required (remember earlier we've set the norm
                            % to 2 to get faster projection computation
                            tmpProbStruct.Tset = Tset;
                            tmpProbStruct.Tconstraint = 2;
                            tmpProbStruct.norm=probStruct.norm;
                            cmOptions = Options;
                            cmOptions.pwa_index = source_dyn;
                            cmOptions.includeLQRset = 0;
                            cmOptions.verbose = 0;
                            cmOptions.noReduce = 0;
                            Matrices=mpt_constructMatrices(tmpSysStruct,tmpProbStruct,cmOptions); 
                        end
                        
                        Matrices.setdist = horizon - 1;
                        Matrices.Tset = Tset;
                        
                        newCSstorage{source_dyn_stack}.dynamics = [newCSstorage{source_dyn_stack}.dynamics source_dyn];
                        newCSstorage{source_dyn_stack}.setdist = [newCSstorage{source_dyn_stack}.setdist horizon-1];
                        newCSstorage{source_dyn_stack}.Matrices{end+1} = Matrices;
                        newCSstorage{source_dyn_stack}.Ubool{end+1} = U;
                        [newCSstorage{source_dyn_stack}.Pfinals,keep] = sub_mpt_expandlist(Pfinal, newCSstorage{source_dyn_stack}.Pfinals);
                        notkept = find(keep==0);
                        if ~isempty(notkept),
                            kept = setdiff(1:length(keep),notkept);
                            newCSstorage{source_dyn_stack}.dynamics(notkept) = [];
                            newCSstorage{source_dyn_stack}.setdist(notkept) = [];
                            newCSstorage{source_dyn_stack}.Matrices = {newCSstorage{source_dyn_stack}.Matrices{kept}};
                            newCSstorage{source_dyn_stack}.Ubool = {newCSstorage{source_dyn_stack}.Ubool{kept}};
                        end
                    end
                    expand(source_dyn) = 1;
                end % expands 
            end % region
        end %source_dyn
    end
    for ii=1:intInfo.stacks,
        [newCSstorage{ii}.Pfinals,keep] = reduceunion(newCSstorage{ii}.Pfinals,Options);
        %%[newCSstorage{ii}.Pfinals,keep] = sub_mpt_reduce_expandlist2(newCSstorage{ii}.Pfinals,Options);
        notkept = find(keep==0);
        if ~isempty(notkept),
            kept = setdiff(1:length(keep),notkept);
            newCSstorage{ii}.dynamics(notkept) = [];
            newCSstorage{ii}.setdist(notkept) = [];
            newCSstorage{ii}.Matrices = {newCSstorage{ii}.Matrices{kept}};
            newCSstorage{ii}.Ubool = {newCSstorage{ii}.Ubool{kept}};
        end
        nTargets = nTargets + length(newCSstorage{ii}.setdist);
        for jj=1:length(newCSstorage{ii}.Pfinals),
            maxPfinal{ii} = sub_mpt_expandlist(newCSstorage{ii}.Pfinals(jj),maxPfinal{ii});
        end
        if isaddnoise,
            maxPfinal{ii} = reduceunion(maxPfinal{ii},Options);
        end
        %%maxPfinal{ii} = sub_mpt_reduce_expandlist2(maxPfinal{ii},Options);
        
        maxPfinalStep{horizon}{ii} = maxPfinal{ii};
    end
    
    Step{horizon} = newCSstorage;
    
    if all(expand==0),
        break
    end
end % horizon

if ~feasible,
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

maxxHull = [maxPfinal{:}];

if Options.feasset,
    truct = maxxHull;
    return
end

disp(sprintf('%d hulls to be explored',nTargets));

%% AL
if Options.PWA_warmend | Options.PWA_savemode
    % Save results thus far
    disp('Saving all hulls..')
    eval(['save ' Options.PWA_savefile 'hulls' num2str(horizon)])
    disp('...done')
end
% Should also add possibility to jump here and load an earlier result file....

disp('Generating regions');

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

pwalist = cell(1,nPWA);
for dyn=1:nPWA,
    pwalist{dyn} = {};
    PP = [];
    targets = targetMatrices{dyn};
    lenT = length(targets);
    for ii=1:lenT,
        if ii==1 | ii==lenT | mod(ii,round(lenT/3))==0,
            fprintf('dynamics %d, hull %d/%d        \r', dyn, ii, lenT);
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
        for qq=1:nR,
            nR=nR+1;
            FFi = [];
            GGi = [];
            boolctr = 0;
            contctr = 0;
            for rr=1:nu,
                if any(rr==ubool),
                    % this input is discrete
                    boolctr = boolctr+1;
                    FFi = [FFi; zeros(1,nx)];
                    GGi = [GGi; U(boolctr)];
                else
                    % this input is continuous
                    contctr = contctr+1;
                    FFi = [FFi; Fi{qq}(contctr,:)];
                    GGi = [GGi; Gi{qq}(contctr)];
                end
            end
            Fi{qq} = FFi;
            Gi{qq} = GGi;
        end
        
        if length(Fi)>0 & isfulldim(Pn(1)),
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


startRedTime = cputime;
disp('Removing unused regions (backwards)...');
cs = {};
for dyn=1:nPWA,
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
        if ii==1 | ii==lenL | mod(ii,round(lenL/3))==0,
            if Options.verbose>0,
                %disp(sprintf('dynamics %d, hull %d/%d', dyn, ii, lenL));
                fprintf('dynamics %d, hull %d/%d         \r', dyn, ii, lenL);
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


if ~userTset,
    % if the initial target sets were merged at the beginning, we put them back to particular list
    for ii=1:length(dynamicsorig),
        cs{dynamicsorig(ii)}.Pn=[cs{dynamicsorig(ii)}.Pn Pnorig(ii)];
        cs{dynamicsorig(ii)}.Fi{end+1}=Fiorig{ii};
        cs{dynamicsorig(ii)}.Gi{end+1}=Giorig{ii};
        cs{dynamicsorig(ii)}.Ai{end+1}=zeros(nx);
        cs{dynamicsorig(ii)}.Bi{end+1}=zeros(1,nx);
        cs{dynamicsorig(ii)}.Ci{end+1}=0;
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
ctrlStruct.details.maxHull = maxxHull;
ctrlStruct.details.N = length(Step)-1;
ctrlStruct.details.invTset = invTset;


return
% and we are done!


% ===========================================================================================
function [expanded,Pfinal] = sub_minSwitch(sysStruct, probStruct, Tset, BMatrices, nx,nu,source_dyn,cmOptions, pOptions),

Matrices = {};

notconverged = 1;
expanded=0;
emptypoly = polytope;
Pfinal = emptypoly;
while notconverged,
    % %     tmpProbStruct = probStruct;
    % %     tmpProbStruct.N = 1;
    % %     tmpProbStruct.Tset = Tset;
    % %     tmpProbStruct.Tconstraint = 2;
    % %     M=mpt_constructMatrices(sysStruct,tmpProbStruct,cmOptions); 

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
function [Pun,keep] = sub_mpt_expandlist2(Pf,Pu,Options),
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
    ctr=0;
    Pun=Pu([1:ii-1,ii+1:lenPu]);         % all indices except of the current one
    P=regiondiff(Pu(ii),Pun,Options);    % if output of regiondiff is empty, Pu(ii) is NOT covered by Pu(I-ii)
    if length(P)==1 & ~isfulldim(P(1)),
        keep(ii)=0;
    end
end

Pun=Pu(find(keep~=0));
return

% ===========================================================================================
function [Pun,keep] = sub_mpt_reduce_expandlist2(Pu,Options),
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
