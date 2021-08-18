function ctrlStruct = mpt_optInfControlPWA(sysStruct, probStruct, Options)
%MPT_OPTINFCONTROLPWA Solves the infinite-time constrained optimal control problem for PWA systems
%
% ctrlStruct=mpt_optInfControlPWA(sysStruct,probStruct,Options),
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the solution to the constraned infinite time optimal control problem
% for a given PWA system 
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
% Options.verbose     - Level of verbosity (see help mpt_init for more details)
% Options.maxiterations - Maximum number of iterations (default is 1000)
% Options.coregen     - If set to 1 (default), core of the infinite time
%                       solution (i.e. the starting part) is generated
%                       automatically by the algorithm
% Options.onlycore    - Breaks and returns once the core has been found
%                       (Default: 0)
% Options.coreiters   - Maximum allowed number of iterations in core generation
% Options.core        -  Starting core can be provided here. It has to be a
%                        valid controller structure object!
% Options.details     - If set to 1, solution of each iteration is stored in the
%                       details fields of the resulting controller structure 
%                       (Default: 0)
% Options.lowmem      - defines memory saving mode
%                         0 - no memory saving - fast computation (default)
%                         1 - slight memory saving
%                         2 - heavy memory saving (slow computation)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrlStruct    - controller structure with the following fields:
%   Pn,Fi,Gi    - for region Pn(i).H*x <= Pn.K(i) computed input is U=Fi{i}*x+Gi{i}   
%   Ai,Bi,Ci    - cost associated to each region (x'Aix + Bi*x + Ci)
%                 Note that Ai and Bi are zero matrices, Ci contains the
%                 step distance to the origin
%   Pfinal      - The maximum control invariant set as a polytope or a polyarray
%   dynamics    - Dynamics active in region Pn(i)
%   details     - A structure with additional details about the solution
%     details.runTime - total run time of the algorithm
%     details.Horizon - a cell array of ctrlStruct's corresponding to each
%                       time step of the algorithm
%
%
% see also MPT_CONTROL, MPT_OPTCONTROLPWA, MPT_ITERATIVEPWA
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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

if nargin<3,
    Options = mptOptions;
end

if ~isfield(Options,'verbose'),
    Options.verbose = mptOptions.verbose;
end

if ~isfield(Options,'abs_tol'),
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options,'maxiterations'),
    Options.maxiterations = 1000;
end
if ~isfield(Options,'details'),
    Options.details = mptOptions.details;
end
if ~isfield(Options,'coregen'),
    Options.coregen=1;
end
if ~isfield(Options,'coreiters'),
    Options.coreiters=100;
end
if ~isfield(Options,'checkcore'),
    % check if core from previous iteration is contained in the newly generated
    % core on the subsequent iteration.
    Options.checkcore=0;
end
if ~isfield(Options,'onlycore')
    Options.onlycore=0;
end
if ~isfield(Options,'oldstuff')
    % this option is for testing purposes only, please do not change this value
    % (or the algorithm runs slower)
    Options.oldstuff=0;
end
if ~isfield(Options,'oldTset')
    % this option is for testing purposes only, please do not change this value
    % (or the algorithm runs slower)
    Options.oldTset=0;
end

if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end

if ~isinf(probStruct.N)
    error('mpt_optInfControlPWA: Horizon must be Inf!');
end

if probStruct.subopt_lev > 0,
    error('mpt_optInfControlPWA: Level of sub-optimality must be 0 !');
end

if probStruct.norm==2,
    error('mpt_optInfControlPWA: Sorry, only linear performance index supported by this function.');
end

origSysStruct = sysStruct;
origProbStruct = probStruct;
if ~iscell(sysStruct.A),
    % LTI system passed, convert it to PWA
    sysStruct = mpt_lti2pwa(sysStruct);
end
[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);
if ubool>0,
    error('mpt_optInfControlPWA: discrete inputs not allowed in this function...');
end

starttime = cputime;

Options.noNoiseOnTset=1;

Step = {};
cs.sysStruct = sysStruct;
cs.probStruct = probStruct;
cs.Pfinal = polytope;
cs.Pn = polytope;
cs.Fi = {};
cs.Gi = {};
cs.Ai = {};
cs.Bi = {};
cs.Ci = {};
cs.dynamics = [];
cs.details.runTime = 0;

cs.overlaps = 0;
emptyCS = cs;
emptypoly = polytope;

CSstorage = {};
for ctr=1:intInfo.stacks,
    CSstorage{ctr}.stack = {};
    CSstorage{ctr}.dynamics = intInfo.dyns_stack{ctr};
end

if ~isfield(Options, 'statusbar')
    Options.statusbar = 0;
end
if ~isfield(Options, 'status_min')
    Options.status_min = 0;
end
if ~isfield(Options, 'status_max')
    Options.status_max = 1;
end
if ~isfield(Options, 'closestatbar'),
    Options.closestatbar = 1;
end
statusbar = Options.statusbar;
closestatbar = Options.closestatbar;
Options.closestatbar = 0;
Options.statusbar = 0;
if statusbar,
    Options.verbose = -1;
end

if statusbar,
    if ~isfield(Options, 'status_handle')
        Options.status_handle = mpt_statusbar('Computing...');
        closestatbar = 1;
    end
end

for ii=1:nPWA,
    % use the Pbnd polytope as initial target set
    inittarget = sysStruct.Pbnd; %& polytope(sysStruct.guardX{ii},sysStruct.guardC{ii});
    cs.Pn = [cs.Pn inittarget];
    cs.Fi{end+1} = zeros(nu,nx);
    cs.Gi{end+1} = Inf;
    cs.Ai{end+1} = {};
    cs.Bi{end+1} = zeros(1,nx);
    cs.Ci{end+1} = 0;
    cs.dynamics = [cs.dynamics ii];
end
Step{1} = cs;

if isfield(probStruct,'xref'),
    origin = probStruct.xref;
else
    origin = zeros(nx,1);
end

converged = 0;
firsttime = 0;
abs_tol = Options.abs_tol;

mplpOptions = Options;
mplpOptions.nu = nu;
mplpOptions.verbose = 0;
roOptions = Options;
inOptions = Options;
inOptions.abs_tol = 0;

if isfield(Options,'core'),
    Options.coregen=0;
end
if Options.onlycore,
    Options.coregen=1;
end
if Options.coregen,
    % generate "core" of the infinite time solution
    if Options.verbose>-1,
        fprintf('\nGenerating the "core"\n\n');
    end
    corestarttime = cputime;
    
    min_progress = 0;
    max_progress = 0.4;
    
    for nn=2:Options.coreiters+1,
        
        prog_min = min_progress + (max_progress - min_progress) * mod(nn-2, 10) / 10;
        prog_max = min_progress + (max_progress - min_progress) * mod(nn-1, 10) / 10;
        
        dyn_min = prog_min;
        dyn_max = prog_min + (prog_max - prog_min) * 0.4;
        ovl_min = dyn_max;
        ovl_max = prog_max;
        
        Step{nn}.cs = {};
        for ii=1:intInfo.stacks,
            CSstorage{ii}.stack = {};
        end
        csctr = 0;
        stepstarttime = cputime;
        if Options.verbose>-1,
            fprintf('---Iteration %d---\n',(nn-1));
        end
        for dyn = 1:nPWA,
            progress = 0;
            p_min = dyn_min + (dyn_max - dyn_min) * (dyn-1) / nPWA;
            p_max = dyn_min + (dyn_max - dyn_min) * dyn / nPWA;
            if statusbar,
                if isempty(mpt_statusbar(Options.status_handle, (dyn-1) / nPWA, dyn_min, dyn_max)),
                    mpt_statusbar;
                    error('Break...');
                end
            end
            
            [isin, inwhich] = isinside(intInfo.Pdyn{dyn}, origin, inOptions);
            if ~isin,
                % dynamics does not contain the origin in it's interior, we can
                % skip it
                continue
            end
            
            CSstack_pos = intInfo.dyns_links(dyn,2);
            if Options.verbose>-1,
                fprintf('exploring dynamics %d     \r', dyn);
            end

            if ~Options.oldTset,
                % form matrices of the 1-step problem
                tmpProbStruct = probStruct;
                tmpProbStruct.Tconstraint=0;
                tmpProbStruct.Tset = emptypoly;
                tmpProbStruct.N=1;
                if nn>2 | 1,
                    if isfield(probStruct,'Qy'),
                        tmpProbStruct.P_N=zeros(ny);
                    else
                        tmpProbStruct.P_N=zeros(nx);
                    end
                end
                localoptions=Options;
                localoptions.verbose=0;
                localoptions.pwa_index=dyn;
                [BMatrices]=mpt_constructMatrices(sysStruct,tmpProbStruct,localoptions);
                if isinf(-BMatrices.W), 
                    % if transition infeasible, continue with next target
                    continue
                end
            end
            
            for reg=1:length(Step{nn-1}.Bi),
                mplpstarttime = cputime;
                
                if mod(reg, 3)==0,
                    progress = (reg-1) / length(Step{nn-1}.Bi);
                    if statusbar,
                        if isempty(mpt_statusbar(Options.status_handle, progress, p_min, p_max)),
                            mpt_statusbar;
                            error('Break...');
                        end
                    end
                end
                
                % targed dynamics
                targetdyn = Step{nn-1}.dynamics(reg);
                
                % target region
                targetregion = Step{nn-1}.Pn(reg);
            
                if Options.oldTset,
                    % form matrices of the 1-step problem
                    tmpProbStruct = probStruct;
                    tmpProbStruct.Tset = targetregion;
                    tmpProbStruct.Tconstraint=2;
                    tmpProbStruct.N=1;
                    if nn>2 | 1,
                        if isfield(probStruct,'Qy'),
                            tmpProbStruct.P_N=zeros(ny);
                        else
                            tmpProbStruct.P_N=zeros(nx);
                        end
                    end
                    localoptions=Options;
                    localoptions.verbose=0;
                    localoptions.pwa_index=dyn;
                    [Matrices]=mpt_constructMatrices(sysStruct,tmpProbStruct,localoptions);
                    
                    if isinf(-Matrices.W), 
                        % if transition infeasible, continue with next target
                        continue
                    end
                else
                    [Matrices,mfeas] = mpt_addTset(sysStruct, BMatrices, targetregion,nx,nu,dyn);
                    if ~mfeas,
                        continue
                    end
                end
                
                % introduce cost to go of the target region
                %S%Matrices.H = Matrices.H + Step{nn-1}.Bi{jj}*[sysStruct.B{ii} zeros(nx,size(Matrices.H,2)-nu)];
                Matrices.H(:,1:nu) = Matrices.H(:,1:nu) + Step{nn-1}.Bi{reg} * sysStruct.B{dyn};
                
                try
                    % solve the 1-step mpLP
                    [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
                catch
                    if Options.verbose>1,
                        disp('Infeasible transition');
                    end
                    continue
                end
                nRinthis = length(Fi);
                if nRinthis==0,
                    if Options.verbose > 1,
                        disp('no regions');
                    end
                    continue
                end
                
                cs = emptyCS;
                cs.Pfinal = Pfinal;
                cs.Pn = Pn;
                cs.Fi = Fi;
                cs.Gi = Gi;
                cs.Ai = cell(1,nRinthis);
                cs.Bi = cell(1,nRinthis);
                cs.Ci = cell(1,nRinthis);
                cs.dynamics = dyn*ones(1,nRinthis);
                
                for qq=1:nRinthis,
                    cs.Ai{qq} = zeros(nx);
                    % add cost-to-go
                    cs.Bi{qq} = details.Bi{qq}(:)' + Step{nn-1}.Bi{reg}*sysStruct.A{dyn};
                    cs.Ci{qq} = details.Ci{qq} + Step{nn-1}.Ci{reg} + Step{nn-1}.Bi{reg}*sysStruct.f{dyn}; 
                end
                cs.overlaps = 0;
                mplpendtime = cputime;
                cs.details.runTime = mplpendtime - mplpstarttime;
                if Options.oldstuff,
                    csctr = csctr + 1;
                    Step{nn}.cs{csctr} = cs;
                else
                    CSstorage{CSstack_pos}.stack{end+1} = cs;
                end
            end  % go through all regions
        end % go through all dynamics
        
        if statusbar,
            if isempty(mpt_statusbar(Options.status_handle, 1, dyn_min, dyn_max)),
                mpt_statusbar;
                error('Break...');
            end
        end
        
        clear nonovlCS
        if Options.oldstuff,
            
            if statusbar,
                roOptions.statusbar = 1;
                roOptions.status_min = ovl_min;
                roOptions.status_max = ovl_max;
                roOptions.status_handle = Options.status_handle;
                roOptions.closestatbar = 0;
            end
            
            Step{nn}.overlaps = 1;
            Step{nn} = mpt_removeOverlaps(Step{nn}.cs, roOptions);
        else
            ctr_ctr = 0;
            for ctr=1:intInfo.stacks,
                ctr_ctr = ctr_ctr + 1;
                if statusbar,
                    roOptions.statusbar = 1;
                    roOptions.status_min = ovl_min + (ctr_ctr-1)*(ovl_max - ovl_min) / length(intInfo.stacks);
                    roOptions.status_max = ovl_min + (ctr_ctr) * (ovl_max - ovl_min) / length(intInfo.stacks);
                    roOptions.status_handle = Options.status_handle;
                    roOptions.closestatbar = 0;
                end
                if ~isempty(CSstorage{ctr}.stack),
                    Cstack = {};
                    for kk=1:length(CSstorage{ctr}.stack),
                        if mpt_isValidCS(CSstorage{ctr}.stack{kk},struct('nowarnings',1)),
                            Cstack{end+1} = CSstorage{ctr}.stack{kk};
                        end
                    end
                    nonovl = mpt_removeOverlaps(Cstack, roOptions);
                    if ~exist('nonovlCS','var'),
                        nonovlCS = nonovl;
                    else
                        nonovlCS = mpt_mergeCS({nonovlCS,nonovl});
                    end
                end
            end
            if ~exist('nonovlCS','var'),
                if statusbar,
                    mpt_statusbar;
                end
                error('mpt_optInfControlPWA: no regions found in this iteration (numerical problems)');
            end
            Step{nn} = nonovlCS;
        end
        
        stependtime=cputime;
        totalregions = length(Step{nn}.Pn);
        if Options.verbose>-1,
            fprintf('Regions: %d\n',totalregions);
        end
        coreregions = 0;

        [in, inwhich] = isinside(Step{nn}.Pn, origin, inOptions);
        if ~in,
            if statusbar
                mpt_statusbar;
            end
            error('no region contains the origin!');
        end

        Step_nn = Step{nn};
        Step_nn_1 = Step{nn-1};

        cs = emptyCS;
        cs.Pn = Step_nn.Pn(inwhich);
        cs.Pfinal = cs.Pn;
        cs.dynamics = Step_nn.dynamics(inwhich);
        cs.overlaps = 0;
        cs.Fi = {Step_nn.Fi{inwhich}};
        cs.Gi = {Step_nn.Gi{inwhich}};
        cs.Ai = {Step_nn.Ai{inwhich}};
        cs.Bi = {Step_nn.Bi{inwhich}};
        cs.Ci = {Step_nn.Ci{inwhich}};
        cs.details = Step_nn.details;
        Step{nn} = cs;

        Step{nn}.details.runTime = stependtime-stepstarttime;
        Step{nn}.details.incore = zeros(size(Step{nn}.dynamics));

        Step_nn = Step{nn};
        Step_nn_1 = Step{nn-1};

        coreregions = 0;
        Rexplored = [];
        for ii=1:length(Step_nn.Fi)
            for jj=1:length(Step_nn_1.Fi),
                if ~isempty(Rexplored),
                    isthere = find(Rexplored==jj);
                    if ~isempty(isthere),
                        continue
                    end
                end
                % check if two regions are the same in terms of shape and value
                % function. If so, this region is a part of the "core"
                if abs(Step_nn.Ci{ii} - Step_nn_1.Ci{jj}) <= abs_tol,
                   if all(all(abs(Step_nn.Bi{ii} - Step_nn_1.Bi{jj}) <= abs_tol)),
                      if Step_nn.Pn(ii) == Step_nn_1.Pn(jj),
                         Step{nn}.details.incore(ii) = 1;
                         coreregions = coreregions + 1;
                         Rexplored = [Rexplored jj];
                         break
                      end
                   end
                end
            end
        end
        
        if coreregions > 0,
            % First check if every region that touches the origin is in the core
            % If it is then the core computation has finished
            if any(Step{nn}.details.incore~=1)
                if Options.verbose>-1,
                    fprintf('\nCore has at least %d regions but we did not yet fully converge\n\n',coreregions);
                end
                coreregions = 0;
            end
        end

        if coreregions > 0,
            Core = emptyCS;
            incore = find(Step{nn}.details.incore==1);
            Core.Pn = Step_nn.Pn(incore);
            Core.Pfinal = Core.Pn;
            Core.Fi = {Step_nn.Fi{incore}};
            Core.Gi = {Step_nn.Gi{incore}};
            Core.Ai = {Step_nn.Ai{incore}};
            Core.Bi = {Step_nn.Bi{incore}};
            Core.Ci = {Step_nn.Ci{incore}};
            Core.dynamics = Step_nn.dynamics(incore);
            Core.details.runTime = cputime - corestarttime;
            Core.details.iterations = nn;
            Core.overlaps = 0;
            converged = 1;
            if Options.verbose>-1,
                fprintf('\nCore consists of %d regions\n',length(Core.Pn));
            end
            
            if Options.onlycore
                ctrlStruct=Core;
                return
            end
            
            % call the same function again with the core
            
            if statusbar,
                Options.statusbar = 1;
                Options.closestatbar = 1;
                Options.status_min = max_progress;
                Options.status_max = 1;
            end
            Options.core = Core;
            clear Step_nn Step_nn_1
            ctrlStruct = mpt_optInfControlPWA(sysStruct, probStruct, Options);
            return
        end
        clear Step_nn Step_nn_1
    end
end

if isfield(Options,'core'),
    % core provided, test if it is a valid controller structure
    Core = Options.core;
    if isa(Core, 'mptctrl'),
        if ~isexplicit(Core),
            error('mpt_optInfControlPWA: Options.core must be an explicit controller!');
        end
        Core = struct(Core);
    end
    try
        msg = evalc('status = mpt_isValidCS(Core);');
    catch
        error('mpt_optInfControlPWA: Options.core is not a valid controller structure!');
        msg = 0;
    end
    if status~=1,
        if statusbar,
            mpt_statusbar;
        end
        error('mpt_optInfControlPWA: Options.core is not a valid controller structure!');
    end
    
    % use the core for exploration
    Step = {};
    Step{1} = Core;
    Options.coregen = 1;
    converged = 1;
end

nocore = 0;

if isfield(Options, 'status_min')
    min_progress = Options.status_min;
else
    min_progress = 0;
end
if isfield(Options, 'status_max')
    max_progress = Options.status_max;
else
    max_progress = 1;
end

if Options.coregen==0 | converged==0,
    % no core generation requested, or core generation failed, or core not provided
    % by options; we use the outside-in approach, i.e. the first target set is
    % the whole state-space of interest
    
    Step={};
    if Options.coregen==1 & converged==0,
        if Options.verbose>-1,
            fprintf('\nCore computation failed to converge, using alternative approach to find it (may take longer)...\n');
        end
        Options.coregen = 0;
    end
    
    Core = emptyCS;
    
    Core.Pn = [Core.Pn sysStruct.Pbnd];
    Core.Fi{1} = zeros(nu,nx);
    Core.Gi{1} = zeros(nu,nx);
    Core.Ai{1} = zeros(nx);
    Core.Bi{1} = zeros(1,nx);
    Core.Ci{1} = 0;
    Core.dynamics = [Core.dynamics 1];

    Step = {};
    Step{1} = Core;
    Ring = Core;
    Ring.Ci{1} = 0;
    Core.Ci{1} = 1e6;
    nocore = 1;
end

origCore = Core;
laststep = Step{end};
laststep.Pfinal = laststep.Pn;
Step = {};
Step{1} = laststep;
Step{1}.details.Core = Core;
if nocore
    Step{1}.details.Ring = Ring;
else
    Step{1}.details.Ring = Core;
end
converged = 0;
finalhorizon = 0;

CSstorage = {};
for ctr=1:intInfo.stacks,
    CSstorage{ctr}.stack = {};
    CSstorage{ctr}.dynamics = intInfo.dyns_stack{ctr};
end

if Options.verbose>-1,
    fprintf('\nExploring the state-space\n\n');
end
for nn = 2:Options.maxiterations,
    if Options.coregen,
        for ii=1:intInfo.stacks,
            CSstorage{ii}.stack = {};
            
            ovl_dyns = [];
            for jj=intInfo.dyns_stack{ii},
                fdyns = [ovl_dyns find(Step{nn-1}.dynamics==jj)];
            end
            % regions from previous step which belong to one of the dynamics
            % which overlap in this segment
            
            CSstorage{ii}.stack{1}.Pn = Step{nn-1}.Pn(fdyns);
            CSstorage{ii}.stack{1}.Pfinal = CSstorage{ii}.stack{1}.Pn;
            CSstorage{ii}.stack{1}.Fi = {Step{nn-1}.Fi{fdyns}};
            CSstorage{ii}.stack{1}.Gi = {Step{nn-1}.Gi{fdyns}};
            CSstorage{ii}.stack{1}.Ai = {Step{nn-1}.Ai{fdyns}};
            CSstorage{ii}.stack{1}.Bi = {Step{nn-1}.Bi{fdyns}};
            CSstorage{ii}.stack{1}.Ci = {Step{nn-1}.Ci{fdyns}};
            CSstorage{ii}.stack{1}.dynamics = Step{nn-1}.dynamics(fdyns);
            CSstorage{ii}.stack{1}.details = Step{nn-1}.details;
            CSstorage{ii}.stack{1}.overlaps = Step{nn-1}.overlaps;
            CSstorage{ii}.stack{1}.sysStruct = Step{nn-1}.sysStruct;
            CSstorage{ii}.stack{1}.probStruct = Step{nn-1}.probStruct;
            if ~isfulldim(CSstorage{ii}.stack{1}.Pn),
                CSstorage{ii}.stack = {};
            end
        end
        if Options.oldstuff,
            Step{nn}.cs{1} = Step{nn-1};
            csctr = 1;
        end
    else
        if Options.oldstuff,
            Step{nn}.cs = {};
            csctr = 0;
        end
        for ii=1:intInfo.stacks,
            CSstorage{ii}.stack = {};
        end
    end
    
    prog_min = min_progress + (max_progress - min_progress) * mod(nn-2, 10) / 10;
    prog_max = min_progress + (max_progress - min_progress) * mod(nn-1, 10) / 10;
    
    dyn_min = prog_min;
    dyn_max = prog_min + (prog_max - prog_min) * 0.4;
    ovl_min = dyn_max;
    ovl_max = prog_max;
    
    if Options.verbose>-1,
        fprintf('---Iteration %d---\n',nn-1);
    end
    stepstarttime=cputime;
    for dyn = 1:nPWA,
        
        progress = 0;
        p_min = dyn_min + (dyn_max - dyn_min) * (dyn-1) / nPWA;
        p_max = dyn_min + (dyn_max - dyn_min) * dyn / nPWA;
        if statusbar,
            if isempty(mpt_statusbar(Options.status_handle, (dyn-1) / nPWA, dyn_min, dyn_max)),
                mpt_statusbar;
                error('Break...');
            end
        end
        
        CSstack_pos = intInfo.dyns_links(dyn,2);
        
        if Options.verbose>-1,
            fprintf('exploring dynamics %d        \r', dyn);
        end
        
        if ~Options.oldTset,
            tmpProbStruct = probStruct;
            tmpProbStruct.Tconstraint=0;
            tmpProbStruct.Tset = emptypoly;
            tmpProbStruct.N=1;
            if nocore==1 & nn==2,
                % keep original P_N at the first iteration
            else
                % final state penalty must be zero as we add cost to go later
                tmpProbStruct.P_N = zeros(nx);
            end
            localoptions=Options;
            localoptions.verbose=0;
            localoptions.pwa_index=dyn;
            [BMatrices]=mpt_constructMatrices(sysStruct,tmpProbStruct,localoptions);
            if isinf(-BMatrices.W), 
                % if transition infeasible, continue with next target
                continue
            end
        end
        
        for reg=1:length(Step{nn-1}.details.Ring.Fi),
            
            if mod(reg, 3)==0,
                progress = (reg-1) / length(Step{nn-1}.details.Ring.Fi);
                if statusbar,
                    if isempty(mpt_statusbar(Options.status_handle, progress, p_min, p_max)),
                        mpt_statusbar;
                        error('Break...');
                    end
                end
            end
            
            mplpstarttime = cputime;
            
            % targed dynamics
            targetdyn = Step{nn-1}.details.Ring.dynamics(reg);
            
            % target region
            targetregion = Step{nn-1}.details.Ring.Pn(reg);
            
            if Options.oldTset,
                % form matrices of the 1-step problem
                tmpProbStruct = probStruct;
                tmpProbStruct.Tset = targetregion;
                tmpProbStruct.Tconstraint=2;
                tmpProbStruct.N=1;
                
                if nocore==1 & nn==2,
                    % keep original P_N at the first iteration
                else
                    % final state penalty must be zero as we add cost to go later
                    tmpProbStruct.P_N = zeros(nx);
                end
                localoptions=Options;
                localoptions.verbose=0;
                localoptions.pwa_index=dyn;
                [Matrices]=mpt_constructMatrices(sysStruct,tmpProbStruct,localoptions);
                if isinf(-Matrices.W), 
                    % if transition infeasible, continue with next target
                    continue
                end
            else
                [Matrices,mfeas] = mpt_addTset(sysStruct, BMatrices, targetregion,nx,nu,dyn);
                if ~mfeas,
                    continue
                end
            end
            % introduce cost to go of the target region
            %S%Matrices.H = Matrices.H + Step{nn-1}.details.Ring.Bi{reg}*[sysStruct.B{dyn} zeros(nx,size(Matrices.H,2)-nu)];
            Matrices.H(:,1:nu) = Matrices.H(:,1:nu) + Step{nn-1}.details.Ring.Bi{reg} * sysStruct.B{dyn};
            
            try
                % solve the 1-step mpLP
                [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
            catch
                if Options.verbose>1,
                    disp('Infeasible transition');
                end
                continue
            end
            nRinthis = length(Fi);
            if nRinthis==0,
                if Options.verbose > 0,
                    if Options.verbose>1,
                        disp('no regions');
                    end
                end
                continue
            end
            
            cs = emptyCS;
            
            cs.Pfinal = Pfinal;
            cs.Pn = Pn;
            cs.Fi = Fi;
            cs.Gi = Gi;
            cs.Ai = cell(1,nRinthis);
            cs.Bi = cell(1,nRinthis);
            cs.Ci = cell(1,nRinthis);
            cs.dynamics = dyn*ones(1,nRinthis);
            for qq=1:nRinthis,
                cs.Ai{qq} = zeros(nx);
                % add cost-to-go
                cs.Bi{qq} = details.Bi{qq}(:)' + Step{nn-1}.details.Ring.Bi{reg}*sysStruct.A{dyn};
                cs.Ci{qq} = details.Ci{qq} + Step{nn-1}.details.Ring.Ci{reg} + Step{nn-1}.details.Ring.Bi{reg}*sysStruct.f{dyn}; 
            end
            cs.overlaps = 0;
            mplpendtime = cputime;
            cs.details.runTime = mplpendtime - mplpstarttime;

            if Options.oldstuff,
                csctr = csctr + 1;
                Step{nn}.cs{csctr} = cs;
            else
                CSstorage{CSstack_pos}.stack{end+1} = cs;
            end
        end % go through all regions
    end % go through all dynamics

    roOptions = Options;

    if statusbar,
        if isempty(mpt_statusbar(Options.status_handle, 1, dyn_min, dyn_max)),
            mpt_statusbar;
            error('Break...');
        end
    end
    
    clear nonovlCS
    ctr_ctr = 0;
    for ctr=1:intInfo.stacks,
        ctr_ctr = ctr_ctr + 1;
        if statusbar,
            roOptions.statusbar = 1;
            roOptions.status_min = ovl_min + (ctr_ctr-1)*(ovl_max - ovl_min) / length(intInfo.stacks);
            roOptions.status_max = ovl_min + (ctr_ctr) * (ovl_max - ovl_min) / length(intInfo.stacks);
            roOptions.status_handle = Options.status_handle;
            roOptions.closestatbar = 0;
        end
        
        if ~isempty(CSstorage{ctr}.stack),
            Cstack = {};
            for kk=1:length(CSstorage{ctr}.stack),
                if mpt_isValidCS(CSstorage{ctr}.stack{kk},struct('nowarnings',1)),
                    Cstack{end+1} = CSstorage{ctr}.stack{kk};
                end
            end
            nonovl = mpt_removeOverlaps(Cstack, roOptions);
            if ~exist('nonovlCS','var'),
                nonovlCS = nonovl;
            else
                nonovlCS = mpt_mergeCS({nonovlCS,nonovl});
            end
        end
    end
    Step{nn} = nonovlCS;

    Step{nn}.overlaps = 0;
    
    Step{nn}.details.runTime = cputime-stepstarttime;
    totalregions = length(Step{nn}.Fi);
    if Options.verbose>-1,
        fprintf('Regions: %d\n',totalregions);
    end
    
    Step{nn}.details.incore = zeros(size(Step{nn}.dynamics));
    coreregions = 0;
    Step_nn = Step{nn};
    Step_nn1 = Step{nn-1};
    Rexplored = [];
    for ii=1:length(Step_nn.Fi),
        Step_nn_Bi_ii = Step_nn.Bi{ii};
        Step_nn_Ci_ii = Step_nn.Ci{ii};
        Step_nn_Fi_ii = Step_nn.Fi{ii};
        Step_nn_Gi_ii = Step_nn.Gi{ii};
        for jj=1:length(Step_nn1.Fi),
            if ~isempty(Rexplored),
                isthere = find(Rexplored==jj);
                if ~isempty(isthere),
                    continue
                end
            end
            % check if two regions are the same in terms of shape and value
            % function. If so, this region is a part
            % of the infinite time solution
            if abs(Step_nn_Ci_ii - Step_nn1.Ci{jj}) <= abs_tol,
                if all(all(abs(Step_nn_Bi_ii - Step_nn1.Bi{jj}) <= abs_tol)),
                    %if all(all(abs(Step_nn_Fi_ii - Step_nn1.Fi{jj}) <= abs_tol)),
                        %if all(abs(Step_nn_Gi_ii - Step_nn1.Gi{jj}) <= abs_tol),
                            if Step_nn.Pn(ii) == Step_nn1.Pn(jj),
                                % regions satisfying the above conditions belog to the "core"
                                Step{nn}.details.incore(ii) = 1;
                                Rexplored = [Rexplored jj];
                                break
                                %continue
                            end
                            %end
                            % end
                end
            end
        end
    end 
    clear Rexplored
    
    incore = find(Step{nn}.details.incore==1);
    inring = find(Step{nn}.details.incore==0);

    Core = emptyCS;

    Core.Pn = Step{nn}.Pn(incore);
    Core.Pfinal = Core.Pn;
    Core.dynamics = Step{nn}.dynamics(incore);
    Core.Fi = {Step{nn}.Fi{incore}};
    Core.Gi = {Step{nn}.Gi{incore}};
    Core.Ai = {Step{nn}.Ai{incore}};
    Core.Bi = {Step{nn}.Bi{incore}};
    Core.Ci = {Step{nn}.Ci{incore}};
    Core.overlaps = 0;
    
    Ring = emptyCS;
    Ring.Pn = Step{nn}.Pn(inring);
    Ring.Pfinal = Ring.Pn;
    Ring.dynamics = Step{nn}.dynamics(inring);
    Ring.Fi = {Step{nn}.Fi{inring}};
    Ring.Gi = {Step{nn}.Gi{inring}};
    Ring.Ai = {Step{nn}.Ai{inring}};
    Ring.Bi = {Step{nn}.Bi{inring}};
    Ring.Ci = {Step{nn}.Ci{inring}};
    Ring.overlaps = 0;
    
    if nocore==1 & sum(incore)>0,
        Options.core = Core;
        Options.coregen = 0;
        if Options.verbose>-1,
            disp('Core found, restarting algorithm...');
        end
        ctrlStruct = mpt_optInfControlPWA(sysStruct, probStruct, Options);
        return
    end
        
    Step{nn}.details.Core = Core;
    Step{nn}.details.Ring = Ring;
    
    if Options.verbose>-1,
        fprintf('Core: %d regions\n', length(Core.Fi));
        fprintf('Ring: %d regions\n', length(Ring.Fi));
    end
    if isempty(Ring.Fi),
        converged = 1;
        finalhorizon = nn;
        break
    end
end

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, 1, 0, 1)),
        mpt_statusbar;
        error('Break...');
    end
end

if ~converged,
    finalhorizon = nn;
    disp('WARNING: maximum number of iterations reached before convergence! Increase Options.maxiterations and try again');
else
    if Options.verbose>-1,
        fprintf('\n\nAlgorithm finished after %d iterations\n',finalhorizon-1);
    end
end
if nocore==1,
    fprintf('\n\nWARNING: Core not found, result is not infinite time solution!\n\n');
end


laststep = Step{end};

ctrlStruct = laststep;
ctrlStruct.sysStruct = origSysStruct;
ctrlStruct.probStruct = origProbStruct;
ctrlStruct.overlaps = 0;

endtime = cputime;
ctrlStruct.details.runTime = endtime-starttime;
ctrlStruct.details.iterations=finalhorizon-1;
ctrlStruct.details.invCore = origCore;
if Options.details,
    for ii=2:finalhorizon
        ctrlStruct.details.Horizon{ii-1} = Step{ii};
    end
else
    ctrlStruct.details = rmfield(ctrlStruct.details,'Core');
    ctrlStruct.details = rmfield(ctrlStruct.details,'Ring');
end

if closestatbar,
    mpt_statusbar;
end
