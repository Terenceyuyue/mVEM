function ctrlStruct = mpt_optControlPWA(sysStruct, probStruct, Options)
%MPT_OPTCONTROLPWA Solves the CFTOC problem for a given PWA system
%
% ctrlStruct=mpt_optControlPWA(sysStruct,probStruct,Options),
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
% Options.checkinf  - If set to 1, regions which belong to infinite time
%                       solution will be identified. Indices of such regions
%                       will then be returned in ctrlStruct.details.incore
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
% see also MPT_CONTROL, MPT_OPTINFCONTROLPWA, MPT_ITERATIVEPWA
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
    Options = [];
end

if ~isfield(Options,'verbose'),
    Options.verbose = mptOptions.verbose;
end
if ~isfield(Options,'details'),
    Options.details = mptOptions.details;
end
if ~isfield(Options,'checkinf')
    Options.checkinf=0;
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options,'oldTset')
    % this option is for testing purposes only, please do not change this value
    % (or the algorithm runs slower)
    Options.oldTset=0;
end
if ~isfield(Options,'maxTime')
    Options.maxTime = Inf;
end


if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end

if isinf(probStruct.N)
    error('mpt_optControlPWA: Horizon must be finite!');
end

if probStruct.subopt_lev > 0,
    error('mpt_optControlPWA: Level of sub-optimality must be 0 !');
end

if probStruct.norm==2,
    error('mpt_optControlPWA: Sorry, only linear performance index supported by this function.');
end

origSysStruct = sysStruct;
origProbStruct = probStruct;
if ~iscell(sysStruct.A),
    % LTI system passed, convert it to PWA
    sysStruct = mpt_lti2pwa(sysStruct);
end
[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);
if ubool>0,
    error('mpt_optControlPWA: discrete inputs not allowed, use mpt_optBoolCtrl or mpt_optMixedCtrl instead...');
end

starttime = cputime;

Options.noNoiseOnTset=1;

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

if probStruct.Tconstraint==2 & isfulldim(probStruct.Tset),
    inittarget = probStruct.Tset;
else
    % use the Pbnd polytope as initial target set
    inittarget = sysStruct.Pbnd;
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

CSstorage = {};
for ctr=1:intInfo.stacks,
    CSstorage{ctr}.stack = {};
    CSstorage{ctr}.dynamics = intInfo.dyns_stack{ctr};
end

mplpOptions = Options;
mplpOptions.nu = nu;
mplpOptions.verbose = 0;
roOptions = Options;

timebreak = 0;

for horizon = 2:probStruct.N+1,
    Step{horizon} = {};
    csctr = 0;
    Pfinals = emptypoly;

    fprintf('--- Step %d (Horizon %d) ---\n',horizon-1,probStruct.N+2-horizon);
    %disp(['---Step ' num2str(horizon-1) '---']);
    stepstarttime=cputime;
    nR = 0;
    nHulls = 0;
    
    for ii=1:intInfo.stacks,
        CSstorage{ii}.stack = {};
    end
    
    for dyn = 1:nPWA,
        CSstack_pos = intInfo.dyns_links(dyn,2);
        fprintf('exploring dynamics %d         \r', dyn);
        
        tmpProbStruct = probStruct;
        if iscell(probStruct.Q),
            tmpProbStruct.Q = origProbStruct.Q{probStruct.N+2-horizon};
        end
        if iscell(probStruct.R),
            tmpProbStruct.R = origProbStruct.R{probStruct.N+2-horizon};
        end
        tmpProbStruct.Tconstraint=0;
        tmpProbStruct.Tset = emptypoly;
        tmpProbStruct.N=1;
        if horizon>2,
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
            disp('Problem infeasible');
            continue
        end
        
        for reg=1:length(Step{horizon-1}.Bi),
            mplpstarttime = cputime;
            
            % targed dynamics
            targetdyn = Step{horizon-1}.dynamics(reg);
            
            % target region
            targetregion = Step{horizon-1}.Pn(reg);
        
            [Matrices,mfeas] = mpt_addTset(sysStruct, BMatrices, targetregion,nx,nu,dyn);
            if ~mfeas,
                continue
            end
                        
            % introduce cost to go of the target region
            %S%Matrices.H = Matrices.H + Step{horizon-1}.Bi{jj}*[sysStruct.B{ii} zeros(nx,size(Matrices.H,2)-nu)];
            Matrices.H(:,1:nu) = Matrices.H(:,1:nu) + Step{horizon-1}.Bi{reg} * sysStruct.B{dyn};
            
            try
                % solve the 1-step mpLP
                [Pn,Fi,Gi,activeConstraints,Pfinal,details]=mpt_mplp(Matrices,mplpOptions);
            catch
                disp('Infeasible transition');
                continue
            end
            if isempty(Fi),
                %if length(Fi)==0 | ~isfulldim(Pn),
                if Options.verbose > 0,
                    disp('no regions');
                end
                continue
            end
            nRinthis = length(Fi);
            nR = nR+nRinthis;
            nHulls = nHulls+1;
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
                cs.Bi{qq} = details.Bi{qq}(:)' + Step{horizon-1}.Bi{reg}*sysStruct.A{dyn};
                cs.Ci{qq} = details.Ci{qq} + Step{horizon-1}.Ci{reg} + Step{horizon-1}.Bi{reg}*sysStruct.f{dyn}; 
            end
            cs.overlaps = 0;
            cs.details.runTime = cputime - mplpstarttime;
            CSstorage{CSstack_pos}.stack{end+1} = cs;
        end % go through all regions
    end % go through all dynamics
    disp(sprintf('%d regions in %d feasible transitions',nR,nHulls));
    
    clear nonovlCS
    for ctr=1:intInfo.stacks,
        if ~isempty(CSstorage{ctr}.stack),
            if nPWA>1,
                Cstack = {};
                for kk=1:length(CSstorage{ctr}.stack),
                    if mpt_isValidCS(CSstorage{ctr}.stack{kk},struct('nowarnings',1)),
                        Cstack{end+1} = CSstorage{ctr}.stack{kk};
                    end
                end
                nonovl = mpt_removeOverlaps(Cstack, roOptions);
            else
                nonovl = CSstorage{ctr}.stack{1};
            end
            if ~exist('nonovlCS','var'),
                nonovlCS = nonovl;
            else
                nonovlCS = mpt_mergeCS({nonovlCS,nonovl});
            end
        end
    end
    Step{horizon} = nonovlCS;

    stependtime=cputime;
    Step{horizon}.details.runTime = stependtime-stepstarttime;
    Step{horizon}.overlaps = 0;
    totalregions = length(Step{horizon}.Pn);
    disp(['Regions: ' num2str(totalregions)]);
    laststep = Step{horizon};    
    
    if Options.checkinf & horizon>1,
        nn = horizon;
        Step{nn}.details.incore = zeros(size(Step{nn}.dynamics));
        coreregions = 0;
        Step_nn = Step{nn};
        Step_nn1 = Step{nn-1};
        for ii=1:length(Step_nn.Fi),
            Step_nn_Bi_ii = Step_nn.Bi{ii};
            Step_nn_Ci_ii = Step_nn.Ci{ii};
            Step_nn_Fi_ii = Step_nn.Fi{ii};
            Step_nn_Gi_ii = Step_nn.Gi{ii};
            for jj=1:length(Step_nn1.Fi),
                % check if two regions are the same in terms of shape, value
                % function and associated control law. If so, this region is a part
                % of the infinite time solution
                if abs(Step_nn_Ci_ii - Step_nn1.Ci{jj}) <= Options.abs_tol,
                    if all(all(abs(Step_nn_Bi_ii - Step_nn1.Bi{jj}) <= Options.abs_tol)),
                        if all(all(abs(Step_nn_Fi_ii - Step_nn1.Fi{jj}) <= Options.abs_tol)),
                            if all(abs(Step_nn_Gi_ii - Step_nn1.Gi{jj}) <= Options.abs_tol),
                                if Step_nn.Pn(ii) == Step_nn1.Pn(jj),
                                    % regions satisfying the above conditions belog to the "core"
                                    Step{nn}.details.incore(ii) = 1;
                                    continue
                                end
                            end
                        end
                    end
                end
            end
        end 
        fprintf('Regions in infinite-time solution: %d\n',sum(Step{nn}.details.incore));
        Step{nn}.details.incore = find(Step{nn}.details.incore==1);
    end
    if cputime - starttime > Options.maxTime,
        timebreak = 1;
        break
    end
end

laststep = Step{end};

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
ctrlStruct.details.timeBreak = timebreak;
