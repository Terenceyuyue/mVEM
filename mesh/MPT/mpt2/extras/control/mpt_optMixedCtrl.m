function ctrlStruct = mpt_optMixedCtrl(sysStruct, probStruct, Options)
%MPT_OPTMIXEDCTRL Computes optimal controller for systems discrete and continuous inputs
%
% ctrlStruct=mpt_optMixedCtrl(sysStruct,probStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the solution of a constraned finite time optimal control problem for
% a given PWA system 
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
% sysStruct         - System structure in the sysStruct format
% probStruct        - Problem structure in the probStruct format
%
% Options.verbose   - Level of verbosity (see help mpt_init for more details)
% Options.details   - If set to 1, solution of each iteration is stored in the
%                       details fields of the resulting controller structure 
%                       (0 by default)
% Options.lowmem    - defines memory saving mode
%                       0 - no memory saving - fast computation (default)
%                       1 - slight memory saving
%                       2 - heavy memory saving (slow computation)
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
% see also MPT_CONTROL, MPT_OPTBOOLCTRL, MPT_OPTCONTROLPWA
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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
    Options = [];
end

if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'details'),
    Options.details = mptOptions.details;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options,'oldstuff')
    % this option is for testing purposes only, please do not change this value
    % (or the algorithm runs slower)
    Options.oldstuff=0;
end
if ~isfield(Options, 'guierrors'),
    Options.guierrors = 0;
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

if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end

Options.noNoiseOnTset=1;

if probStruct.norm==2,
    error('mpt_optMixedCtrl: Sorry, only linear performance index supported by this function.');
end
if isinf(probStruct.N),
    error('mpt_optMixedCtrl: Prediction horizon must be finite!');
end
if ~isfield(sysStruct,'Uset'),
    error('mpt_optMixedCtrl: sysStruct.Uset must be given!');
end

starttime = cputime;

origSysStruct = sysStruct;
origProbStruct = probStruct;
if ~iscell(sysStruct.A),
    % convert LTI system to PWA format
    sysStruct = mpt_lti2pwa(sysStruct);
end

[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);

if nbool==0,
    error('mpt_optMixedCtrl: At least one input must be discrete!');
end

Step = {};
cs.sysStruct = sysStruct;
cs.probStruct = probStruct;
cs.Pfinal = polytope;
cs.Pn = [];
cs.Fi = {};
cs.Gi = {};
cs.Ai = {};
cs.Bi = {};
cs.Ci = {};
cs.dynamics = [];
emptyCS = cs;
emptypoly = polytope;

Pn = sysStruct.Pbnd;
ii=1;
if probStruct.Tconstraint==2 & isfulldim(probStruct.Tset),
    inittarget = probStruct.Tset;
else
    % use the Pbnd polytope as initial target set
    inittarget = Pn;
end
cs.Pn = [cs.Pn inittarget];
for ii=1:length(inittarget),
    cs.Fi{end+1} = {};
    cs.Gi{end+1} = {};
    cs.Ai{end+1} = {};
    cs.Bi{end+1} = zeros(1,nx);
    cs.Ci{end+1} = 0;
    cs.dynamics = [cs.dynamics ii];
end

Step{1}=cs;

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

if statusbar,
    statbar.handle = mpt_statusbar('Computing...');
    statbar.progress = 0;
    Options.verbose = -1;
end

mplpOptions = Options;
mplpOptions.verbose = 0;
mplpOptions.nu = ncont;

CSstorage = {};
for ctr=1:intInfo.stacks,
    CSstorage{ctr}.stack = {};
    CSstorage{ctr}.dynamics = intInfo.dyns_stack{ctr};
end

if any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax)),
    dUconstraints = 1;
else
    dUconstraints = 0;
end
if dUconstraints,
    error('mpt_optMixedCtrl: Constraints on delta U not supported in this version!');
end
if probStruct.tracking==1
    if Options.guierrors,
        error('Tracking with delta U formulation for boolean inputs not supported. Set "tracking without delta U formulation" option in "Advanced options".');
    else
        error('Tracking with delta U formulation for boolean inputs not supported. Please use "probStruct.tracking=2".');
    end
end
if isfield(sysStruct, 'dumode'),
    error('Sorry, deltaU constraints are not supported for systems with boolean inputs!');
end
if isfield(probStruct,'Rdu')
    if ~all(probStruct.Rdu==0),
        error('mpt_optMixedCtrl: non-zero probStruct.Rdu is not supported in this version!');
    end
end

dUconstraints = 0;
dUctr = 0;

for horizon = 2:probStruct.N+1,
    Step{horizon}.cs = {};
    if dUconstraints,
        if dUctr==0,
            for ii=1:intInfo.stacks,
                CSstorage{ii}.stack = {};
            end
        end
    else
        for ii=1:intInfo.stacks,
            CSstorage{ii}.stack = {};
        end
    end
    
    csctr = 0;
    Pfinals = emptypoly;
    if Options.verbose>-1,
        fprintf('--- Step %d (Horizon %d) ---\n',horizon-1,probStruct.N+2-horizon);
    end
    stepstarttime=cputime;
    nR = 0;
    nHulls = 0;
    
    min_progress = (horizon - 2)/probStruct.N;
    max_progress = (horizon - 1)/probStruct.N;
    
    prog_min = min_progress;
    prog_max = min_progress + 0.1*max_progress;
    p_cnt = 0;
    progress = 0;

    for possibleU = 1:size(Ucombs,1),
        
        p_cnt = p_cnt + 1;
        progress = (p_cnt - 1) / size(Ucombs, 1);
        if statusbar,
            if isempty(mpt_statusbar(statbar.handle, progress, prog_min, prog_max)),
                mpt_statusbar;
                error('Break...');
            end     
        end
        
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
                if Options.verbose>-1,
                    disp(['Excluding input ' mat2str(U')]);
                end
                continue
            end
        end
        if isempty(DDyn) & Options.verbose>-1,
            disp(['No guardline associated to input ' mat2str(U')]);
        end
        for dyn = DDyn,
            CSstack_pos = intInfo.dyns_links(dyn,2);
            if Options.verbose>-1,
                fprintf('U=%s, dynamics=%d             \r', mat2str(U'),dyn);
            end
            
            % construct matrices for 1-step poblem
            % we insert target set later
            tmpProbStruct = probStruct;
            tmpProbStruct.Tset = emptypoly;
            tmpProbStruct.Tconstraint=0;
            tmpProbStruct.N=1;
            
            if statusbar,
                if isempty(mpt_statusbar(statbar.handle, progress, prog_min, prog_max)),
                    mpt_statusbar;
                    error('Break...');
                end     
            end

            if horizon>2,
                if isfield(probStruct,'Qy'),
                    tmpProbStruct.P_N=zeros(1,ny);
                else
                    tmpProbStruct.P_N=zeros(1,nx);
                end
            end
            tmpProbStruct.R = probStruct.R(ucont,ucont);
            if isfield(tmpProbStruct,'uref'),
                tmpProbStruct.uref = probStruct.uref(ucont);
            end
            
            tmpSysStruct = sysStruct;
            tmpSysStruct.B{dyn} = sysStruct.B{dyn}(:,ucont);
            tmpSysStruct.f{dyn} = sysStruct.f{dyn} + sysStruct.B{dyn}(:,ubool)*U;
            tmpSysStruct.D{dyn} = sysStruct.D{dyn}(:,ucont);
            if isfield(sysStruct, 'Dy'),
                tmpSysStruct.Dy{dyn} = sysStruct.Dy{dyn}(:,ucont);
            end
            tmpSysStruct.g{dyn} = sysStruct.g{dyn} + sysStruct.D{dyn}(:,ubool)*U;
            
            % eliminate equality constraints from guards on U
            tmpSysStruct.guardC{dyn} = pureX_gC{dyn};
            tmpSysStruct.guardU{dyn} = pureX_gU{dyn}(:,ucont);
            tmpSysStruct.guardX{dyn} = pureX_gX{dyn};
            tmpSysStruct.umax = sysStruct.umax(ucont);
            tmpSysStruct.umin = sysStruct.umin(ucont);
            tmpSysStruct.dumax = sysStruct.dumax(ucont);
            tmpSysStruct.dumin = sysStruct.dumin(ucont);
            
            localoptions=Options;
            localoptions.pwa_index=dyn;
            localoptions.verbose=0;
            
            [BMatrices]=mpt_constructMatrices(tmpSysStruct,tmpProbStruct,localoptions);
            if isinf(-BMatrices.W), 
                % if transition infeasible, continue with next target
                continue
            end
                
            for reg=1:length(Step{horizon-1}.Bi),
                mplpstarttime = cputime;
                % target region
                targetregion = Step{horizon-1}.Pn(reg);
                
                if statusbar,
                    if isempty(mpt_statusbar(statbar.handle, progress, prog_min, prog_max)),
                        mpt_statusbar;
                        error('Break...');
                    end     
                end
                
                deltaU = 0;
                if horizon>2,
                    % control move valid in target region
                    targetU = Step{horizon-1}.Gi{reg}(ubool);
                    deltaU = U - targetU;
                    if any(deltaU<sysStruct.dumin(ubool)) | any(deltaU>sysStruct.dumax(ubool)),
                        % new control action violates slew constraints, continue
                        % with next sequence of inputs
                        continue
                    end
                end
                
                [Matrices,mfeas] = mpt_addTset(tmpSysStruct, BMatrices, targetregion,nx,ncont,dyn);
                if ~mfeas,
                    continue
                end
                
                % include cost-to-go
                %Matrices.H(:,1:nu) = Matrices.H(:,1:nu) + Step{horizon-1}.Bi{reg} * sysStruct.B{dyn};
                Matrices.H(:,1:ncont) = Matrices.H(:,1:ncont) + Step{horizon-1}.Bi{reg} * tmpSysStruct.B{dyn};
                
                [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
                
                if statusbar,
                    if isempty(mpt_statusbar(statbar.handle, progress, prog_min, prog_max)),
                        mpt_statusbar;
                        error('Break...');
                    end     
                end
                
                nRinthis = length(Fi);
                if nRinthis<1,
                    if Options.verbose>-1,
                        disp('No regions');
                    end
                    continue
                end
                    
                nHulls = nHulls+1;
                cs = emptyCS;
                cs.Pfinal = Pfinal;
                cs.Pn = Pn;
                cs.Fi = cell(1,nRinthis);
                cs.Gi = cell(1,nRinthis);
                cs.Ai = cell(1,nRinthis);
                cs.Bi = cell(1,nRinthis);
                cs.Ci = cell(1,nRinthis);
                cs.dynamics = dyn*ones(1,nRinthis);
                for qq=1:nRinthis,
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
                    cs.Fi{qq} = FFi;
                    cs.Gi{qq} = GGi;
                    cs.Ai{qq} = zeros(nx);
                    % add cost-to-go
                    cs.Bi{qq} = details.Bi{qq} + Step{horizon-1}.Bi{reg}*sysStruct.A{dyn};
                    CCi = norm(probStruct.R(ubool,ubool)*U,probStruct.norm);
                    if isfield(probStruct,'Rdu'),
                        CCi = CCi + norm(probStruct.Rdu(ubool,ubool)*deltaU,probStruct.norm);
                    end
                    cs.Ci{qq} =  CCi + details.Ci{qq} + Step{horizon-1}.Ci{reg} + Step{horizon-1}.Bi{reg}*sysStruct.f{dyn}; 
                end
                cs.overlaps = 0;
                cs.details.runTime = cputime - mplpstarttime;
                csctr = csctr + 1;
                if Options.oldstuff,
                    Step{horizon}.cs{csctr} = cs;
                else
                    CSstorage{CSstack_pos}.stack{end+1} = cs;
                end
                
                clear cs
            end
        end
    end
    if Options.verbose>-1,
        fprintf('%d regions in %d feasible transitions\n',nR,nHulls);
    end
    if nR==0,
        break
    end
    roOptions = Options;

    if statusbar,
        if isempty(mpt_statusbar(statbar.handle, progress, prog_min, prog_max)),
            mpt_statusbar;
            error('Break...');
        end     
    end

    if dUconstraints,
        dUctr = dUctr;
        if dUctr==2,
            % remove overlaps
            for ctr=1:intInfo.stacks,
                if statusbar,
                    prog_min = prog_max;
                    prog_max = prog_min + (max_progress-prog_min) * (ctr/(length(intInfo.stacks)+1));
                    roOptions.statusbar = 1;
                    roOptions.status_min = prog_min;
                    roOptions.status_max = prog_max;
                    roOptions.status_handle = statbar.handle;
                    roOptions.closestatbar = 0;
                end
                
                if ~isempty(CSstorage{ctr}.stack),
                    nonovl = mpt_removeOverlaps(CSstorage{ctr}.stack, roOptions);
                    if ctr==1,
                        nonovlCS = nonovl;
                    else
                        nonovlCS = mpt_mergeCS({nonovlCS,nonovl});
                    end
                end
            end
            Step{horizon} = nonovlCS;
            dUctr = 0;
        else
            % just add current partitions to controllers computed in previous
            % step
            for ctr=1:intInfo.stacks,
                if ~isempty(CSstorage{ctr}.stack),
                    if ctr==1,
                        ovlCS = CSstorage{ctr}.stack;
                    else
                        ovlCS = mpt_mergeCS({ovlCS,CSstorage{ctr}.stack});
                    end
                    %%dUstorage{ctr}.stack = mpt_mergeCS({dUstorage{ctr}.stack,CSstorage{ctr}.stack});
                end
            end
            Step{horizon} = ovlCS;
            dUctr = dUctr+1;
        end
    else
        clear nonovlCS
        for ctr=1:intInfo.stacks,
            
            if statusbar,
                prog_min = prog_max;
                prog_max = prog_min + (max_progress-prog_min) * (ctr/(length(intInfo.stacks)+1));
                roOptions.statusbar = 1;
                roOptions.status_min = prog_min;
                roOptions.status_max = prog_max;
                roOptions.status_handle = statbar.handle;
                roOptions.closestatbar = 0;
            end
            
           if ~isempty(CSstorage{ctr}.stack),
                Cstack = {};
                for kk=1:length(CSstorage{ctr}.stack),
                    if mpt_isValidCS(CSstorage{ctr}.stack{kk},struct('nowarnings',1)),
                        Cstack{end+1} = CSstorage{ctr}.stack{kk};
                    end
                end
                %%nonovl = mpt_removeOverlaps(CSstorage{ctr}.stack, roOptions);
                nonovl = mpt_removeOverlaps(Cstack, roOptions);
                if ~exist('nonovlCS','var'),
                    nonovlCS = nonovl;
                else
                    nonovlCS = mpt_mergeCS({nonovlCS,nonovl});
                end
            end
        end
        Step{horizon} = nonovlCS;
    end
    
    if statusbar,
        if isempty(mpt_statusbar(statbar.handle, 1, min_progress, max_progress)),
            mpt_statusbar;
            error('Break...');
        end     
    end
    
    stependtime=cputime;
    Step{horizon}.details.runTime = stependtime-stepstarttime;
    totalregions = length(Step{horizon}.Pn);
    if Options.verbose>-1,
        disp(['Regions: ' num2str(totalregions)]);
    end
    laststep = Step{horizon};    
end

if nR==0,
    laststep = Step{end-1};
else
    laststep = Step{end};
end

ctrlStruct = laststep;
ctrlStruct.sysStruct = origSysStruct;
ctrlStruct.probStruct = origProbStruct;
ctrlStruct.overlaps = 0;

endtime = cputime;
ctrlStruct.details.runTime = endtime-starttime;
if Options.details,
    for ii=2:probStruct.N+1
        ctrlStruct.details.Horizon{ii-1} = Step{ii};
    end
end

if closestatbar,
    mpt_statusbar;
end