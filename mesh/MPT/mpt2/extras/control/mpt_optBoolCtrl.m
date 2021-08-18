function ctrlStruct = mpt_optBoolCtrl(sysStruct, probStruct, Options)
%MPT_OPTBOOLCTRL Computes optimal controller for systems with discrete inputs
%
% ctrlStruct = mpt_optBoolCtrl(sysStruct,probStruct)
% ctrlStruct = mpt_optBoolCtrl(sysStruct,probStruct,Options)
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
% All inputs are assumed to be discrete as defined in sysStruct.Uset
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
%                       0 - no memory saving - fast computation
%                       1 - slight memory saving (default)
%                       2 - heavy memory saving (slow computation)
%                     If set to 1 or 2, intermediate data are stored to disk
%                     before entering the critical part of the algorithm.
%                     Data is stored in a file in your temporary directory and
%                     deleted if no error occures. Ensure that you have write
%                     permission to your temp directory (use command tempdir)
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
%   details     - A structure with additional details about the solution
%     details.runTime - total run time of the algorithm
%     details.Horizon - a cell array of ctrlStruct's corresponding to each
%                       time step of the algorithm
%
%
% see also MPT_CONTROL, MPT_OPTCONROLPWA
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
if ~isfield(Options,'rel_tol'),
    Options.rel_tol = mptOptions.rel_tol;
end
if ~isfield(Options, 'guierrors'),
    Options.guierrors = 0;
end
if ~isfield(Options,'oldstuff')
    % this option is for testing purposes only, please do not change this value
    % (or the algorithm runs slower)
    Options.oldstuff=0;
end
if ~isfield(Options,'y0boundshack')
    Options.y0boundshack=0;
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
origSysStruct = sysStruct;

Options.noNoiseOnTset=1;

if probStruct.norm==2,
    error('mpt_optBoolCtrl: Sorry, only linear performance index supported by this function.');
end
if isinf(probStruct.N),
    error('mpt_optBoolCtrl: Prediction horizon must be finite!');
end
if ~isfield(sysStruct,'Uset'),
    error('mpt_optBoolCtrl: sysStruct.Uset must be given!');
end
if isfield(probStruct,'P_N'),
    if any(any(probStruct.P_N~=0)),
        error('mpt_optBoolCtrl: probStruct.P_N must be a zero matrix!');
    end
end

starttime = cputime;

origSysStruct = sysStruct;
origProbStruct = probStruct;
if ~iscell(sysStruct.A),
    % convert LTI system to PWA format
    sysStruct = mpt_lti2pwa(sysStruct);
end

[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);

if nbool~=nu,
    error(['mpt_optBoolCtrl: All inputs must be integer! sysStruct.Uset has to be a cell with ' num2str(nu) ' elements.']);
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
            
if iscell(sysStruct.Uset),
    if length(sysStruct.Uset)~=nu,
        error(['mpt_optBoolCtrl: All inputs must be integer! sysStruct.Uset has to be a cell with ' num2str(nu) ' elements.']);
    end
elseif nu~=1,
    error('mpt_optBoolCtrl: All inputs must be integer! Please set sysStruct.Uset');
end

Ucombs = mpt_allcombs(sysStruct.Uset);

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
    error('mpt_optBoolCtrl: Constraints on delta U not supported in this version!');
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
        error('mpt_optBoolCtrl: non-zero probStruct.Rdu is not supported in this version!');
    end
end
dUctr = 0;

if statusbar,
    statbar.handle = mpt_statusbar('Computing...');
    statbar.progress = 0;
    Options.verbose = -1;
end

mplpOptions = Options;
mplpOptions.verbose = 0;
mplpOptions.nu = nu;
dmOptions = Options;
dmOptions.noReduce = 1;

if dUconstraints,
    finalHorizon = probStruct.N+2;
else
    finalHorizon = probStruct.N+1;
end

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
            elseif all(pureU_gU{ii}*U - pureU_gC{ii}<=Options.abs_tol),
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
        end
        for dyn = DDyn,
            if Options.verbose>-1,
                fprintf('U=%s, dynamics=%d             \r', mat2str(U'),dyn);
            end
            CSstack_pos = intInfo.dyns_links(dyn,2);

            if statusbar,
                if isempty(mpt_statusbar(statbar.handle, progress, prog_min, prog_max)),
                    mpt_statusbar;
                    error('Break...');
                end     
            end

            A = sysStruct.A{dyn};
            B = sysStruct.B{dyn};
            f = sysStruct.f{dyn};
            
            % include constraints on 'y'
            C = sysStruct.C{dyn};
            D = sysStruct.D{dyn};
            g = sysStruct.g{dyn};
            
            tmpSysStruct = sysStruct;
            tmpSysStruct.f{dyn} = f + B*U;
            tmpSysStruct.B{dyn} = zeros(size(B));
            
            
            % eliminate equality constraints from guards on U
            tmpSysStruct.guardC{dyn} = pureX_gC{dyn};
            tmpSysStruct.guardU{dyn} = pureX_gU{dyn};
            tmpSysStruct.guardX{dyn} = pureX_gX{dyn};
            tmpSysStruct.guardU{dyn} = [tmpSysStruct.guardU{dyn}; eye(nu); -eye(nu)];
            tmpSysStruct.guardC{dyn} = [tmpSysStruct.guardC{dyn}; U+Options.rel_tol*2; -(U-Options.rel_tol*2)];
            tmpSysStruct.guardX{dyn} = [tmpSysStruct.guardX{dyn}; zeros(2*nu,nx)];
            
            PPfinal = Pdynamics{dyn};
            if Options.y0boundshack & probStruct.y0bounds,
                [Hf,Kf]=double(PPfinal);
                Hn = [Hf; C; -C];
                Kn = [Kf; sysStruct.ymax-D*U-g; -sysStruct.ymin+D*U+g];
                PPfinal = polytope(Hn,Kn);
            end
            
            for reg=1:length(Step{horizon-1}.Bi),
                
                if statusbar,
                    if isempty(mpt_statusbar(statbar.handle, progress, prog_min, prog_max)),
                        mpt_statusbar;
                        error('Break...');
                    end     
                end

                mplpstarttime = cputime;
                % target region
                targetregion = Step{horizon-1}.Pn(reg);

                deltaU = 0;
                if horizon>2,
                    % control move valid in target region
                    targetU = Step{horizon-1}.Gi{reg};
                    deltaU = U - targetU;
                    if any(deltaU<sysStruct.dumin) | any(deltaU>sysStruct.dumax),
                        % new control action violates slew constraints, continue
                        % with next sequence of inputs
                        continue
                    end
                end
                
                if Options.y0boundshack & ~probStruct.y0bounds,
                    [Hf,Kf]=double(targetregion);
                    Hn = [Hf; C; -C];
                    Kn = [Kf; sysStruct.ymax-D*U-g; -sysStruct.ymin+D*U+g];
                    targetregion = polytope(Hn,Kn);
                end
                
                % compute all states which can be driven from everywhere where
                % dynamics 'dyn' is active to 'targetregion' in one time step
                % using control move 'U'

                [Pret,kr,feasible] = domain(targetregion,A,B*U+f,PPfinal,1,dmOptions);
                % if this set is not empty, we know there exists set of states which
                % can be driven into 'targetregion' using fixed U
  
                if feasible,
                    Pret = reduce(Pret);
                    tmpProbStruct = probStruct;
                    tmpProbStruct.Tset = targetregion;
                    tmpProbStruct.Tconstraint=2;
                    tmpProbStruct.N=1;
                    if horizon>2,
                        if isfield(probStruct,'Qy'),
                            tmpProbStruct.P_N=zeros(1,ny);
                        else
                            tmpProbStruct.P_N=zeros(1,nx);
                        end
                    end
                    tmpProbStruct.R=zeros(nu);
                    tmpProbStruct.y0bounds=0;
                    
                    tmpSysStruct.Pbnd = Pret;
                    
                    localoptions=Options;
                    localoptions.pwa_index=dyn;
                    localoptions.verbose=0;
                    
                    [Matrices]=mpt_constructMatrices(tmpSysStruct,tmpProbStruct,localoptions);
                    if isinf(-Matrices.W), 
                        % if transition infeasible, continue with next target
                        continue
                    end
   
                    % include cost-to-go
                    %Matrices.H = Matrices.H + Step{horizon-1}.Bi{jj}*[sysStruct.B{dyn} zeros(nx,size(Matrices.H,2)-nu)];
                    Matrices.H(:,1:nu) = Matrices.H(:,1:nu) + Step{horizon-1}.Bi{reg} * sysStruct.B{dyn};
                    
                    [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
                    
                    if statusbar,
                        if isempty(mpt_statusbar(statbar.handle, progress, prog_min, prog_max)),
                            mpt_statusbar;
                            error('Break...');
                        end     
                    end
                    
                    nRinthis = length(Pn);
                    if nRinthis<1,
                        disp('No regions');
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

                        cs.Fi{qq} = zeros(nu,nx);
                        cs.Gi{qq} = U;
                        cs.Ai{qq} = zeros(nx);

                        % add cost-to-go
                        cs.Bi{qq} = details.Bi{qq} + Step{horizon-1}.Bi{reg}*sysStruct.A{dyn};
                        CCi = norm(probStruct.R*U,probStruct.norm);
                        if isfield(probStruct,'Rdu'),
                            CCi = CCi + norm(probStruct.Rdu*deltaU,probStruct.norm);
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
    end
    if Options.verbose > -1
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
        if statusbar,
            prog_min = prog_max;
            prog_max = max_progress;
            roOptions.statusbar = 1;
            roOptions.status_min = prog_min;
            roOptions.status_max = prog_max;
            roOptions.status_handle = statbar.handle;
            roOptions.closestatbar = 0;
        end
        
        dUctr = dUctr;
        if dUctr==1 | horizon==2,
            % remove overlaps
            for ctr=1:intInfo.stacks,
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
                        ovlCS = mpt_mergeCS(CSstorage{ctr}.stack);
                    else
                        ovlCS = mpt_mergeCS({ovlCS,CSstorage{ctr}.stack{:}});
                    end
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
try
    if Options.details,
        for ii=2:probStruct.N+1
            ctrlStruct.details.Horizon{ii-1} = Step{ii};
        end
    end
end

if closestatbar,
    mpt_statusbar;
end