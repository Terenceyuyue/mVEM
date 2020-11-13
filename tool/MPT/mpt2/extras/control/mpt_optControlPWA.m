function ctrlStruct = mpt_optControlPWA(sysStruct, probStruct, Options)
% MPT_OPTCONTROLPWA Solves the CFTOC problem for a given PWA system
%
% ctrlStruct = mpt_optControlPWA(sysStruct,probStruct)
% ctrlStruct = mpt_optControlPWA(sysStruct,probStruct,Options)
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
%                     details fields of the resulting controller structure
%                     (0 by default)
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
% ctrlStruct            - Controller structure with following fields:
%   Pn,Fi,Gi            - for region Pn(i).H*x <= Pn.K(i) computed input is U=Fi{i}*x+Gi{i}
%   Ai,Bi,Ci            - cost associated to each region (x'Aix + Bi*x + Ci)
%                         Note that Ai and Bi are zero matrices, Ci contains the
%                         step distance to the origin
%   Pfinal              - The maximum control invariant set as a polytope or a polyarray
%   dynamics            - Dynamics active in region Pn(i)
%   details             - A structure with additional details about the solution
%   details.runTime     - total run time of the algorithm
%   details.mplpRunTime - time spent for solving mpLP
%   details.roRuntime   - time spent for detecting and removing overlaps
%   details.Horizon     - a cell array of ctrlStruct's corresponding to each
%                         time step of the algorithm (only returned if
%                         Options.details > 0)
%
% see also MPT_CONTROL, MPT_OPTINFCONTROLPWA, MPT_ITERATIVEPWA

% Copyright is with the following author(s):
%
% (C) 2004 Miroslav Baric, Automatic Control Laboratory, ETH Zurich,
%          baric@control.ee.ethz.ch
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
%
global mptOptions;

true = 1;
false = 0;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<2,
    error('mpt_optControlPWA: Wrong number of input arguments!');
end

if nargin < 3,
    Options = mptOptions;
end

if ~isfield(Options, 'mplpver')
    Options.mplpver = Inf;  % choose the fastest
end
if ~isfield(Options, 'qpsolver'),
    Options.qpsolver = mptOptions.qpsolver;
end

if Options.mplpver < 4 | Options.qpsolver < 0,
    fprintf('WARNING: older version of mpLP solver specified, switching to different algorithm...\n');
    ctrlStruct = mpt_optControlPWAold(sysStruct, probStruct, Options);
    return
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


if ~isfield(Options,'details'),
    Options.details = mptOptions.details;
end
if ~isfield(Options,'verbose'),
    Options.verbose = mptOptions.verbose;
end
if ~isfield(Options,'recedingHorizon'),
    %Options.recedingHorizon = false;
    Options.recedingHorizon = (Options.details == 0);
elseif ( Options.recedingHorizon ),
    disp ('Calculating receding horizon control policy.');
end

if statusbar,
    statbar.handle = mpt_statusbar('Computing...');
    statbar.progress = 0;
    Options.verbose = -1;
end

if ~isfield(sysStruct,'verified')
    sysStruct = mpt_verifySysStruct(sysStruct);
end
if ~isfield(probStruct,'verified')
    probStruct = mpt_verifyProbStruct(probStruct);
end

if isfield(sysStruct, 'noise'),
    if mpt_isnoise(sysStruct.noise),
        error('No optimal solution for systems with additive noise available.');
    end
end

if isinf(probStruct.N)
    error('mpt_optControlPWA: Horizon must be finite!');
end

if probStruct.subopt_lev > 0,
    error('mpt_optControlPWA: Level of sub-optimality must be 0 !');
end
if ~iscell(sysStruct.A),
    % LTI system passed, convert it to PWA
    sysStruct = mpt_lti2pwa(sysStruct);
end
%
if probStruct.norm==2,
    error('mpt_optControlPWA: Sorry, only linear performance index supported.');
end

if isfield(probStruct,'P_N'),
    tmpProbStruct.P_N = probStruct.P_N;
else
    disp('No terminal cost defined. Using the default.');
end

%
% If no terminal region is given ...
%
if isfield(probStruct,'Tset')
    probStruct.Tconstraint = 2; % user defined terminal set
else
    probStruct.Tset = polytope;
end
%
% polytopic uncertainty
%
if isfield(sysStruct,'Aunc') | isfield(sysStruct,'Bunc'),
    error('mpt_optControlPWA: Parametric uncertainty not supported for PWA systems!');
end
%
% additive uncertainty
%
Options.noNoiseOnTset = 1;
if isfield(sysStruct,'noise') == 1,
    if ~mpt_isnoise(sysStruct.noise),
        disp('Noise polytope empty.');
        isNoisy = false;
    else
        if isa(sysStruct.noise, 'polytope'),
            noisedim = dimension(sysStruct.noise);
        else
            % NOTE! remember that the noise in V-representation has vertices
            % stored column-wise!!!
            noisedim = size(sysStruct.noise, 1);
        end
        if noisedim ~= size(sysStruct.A{1},1),
            error('Noise polytope of the wrong dimension!');
        else
            isNoisy = true;
            if isa(sysStruct.noise, 'polytope'),
                Vnoise = extreme(sysStruct.noise);
            else
                % NOTE! remember that the noise in V-representation has vertices
                % stored column-wise!!! therefore we need to convert it to a
                % row-wise format, since that is what extreme() returns!
                Vnoise = sysStruct.noise';
            end
            Options.noNoiseOnTset = 0;
        end
    end
else
    isNoisy = false;
end

[nx,nu,ny,nPWA] = mpt_sysStructInfo(sysStruct);
horizon = probStruct.N;

% initialize locally used "option" structures
%
localOptions     = Options;
roOptions        = Options;
mplpOptions      = Options;
mplpOptions.nu   = nu;
mplpOptions.verbose = -1;
if ( Options.recedingHorizon ),
    mplpOptions.skipDualDegenerate = true;
end
%
tmpProbStruct   = probStruct;
tmpProbStruct.N = 1;           % horizon=1 (step of dyn. programming)

% create basic data structures
%
cs.sysStruct  = sysStruct;
cs.probStruct = probStruct;
cs.Pfinal     = [];
cs.details    = struct('mplpRunTime',0,'roRunTime',0,'runTime',0);
cs.dynamics   = [];
cs.overlaps   = 0;
cs.Pn         = [];
cs.Ai         = {};
cs.Bi         = {};
cs.Ci         = {};
cs.Fi         = {};
cs.Gi         = {};
%
Step = cell(1,horizon+1);
for ii = 1:horizon,
    Step{ii}.ctrlStruct{1}       = cs;
    Step{ii}.mergedCtrlStruct    = cs;
    Step{ii}.nPartitions    = 0;
    Step{ii}.validPartitions    = [];
end
%
empty_polytope = polytope;
%
% constraint information - contains adjacency information and
% target set constraints
%
adjInfo        = struct ('adjacencyList', [], ...
    'tSetList',      []);
%
constraintInfo = struct ('idxValueFConstr' , [], ...
    'idxTSetConstr'   , [], ...
    'idxSysConstr',     [], ...
    'adjacencyInfo', adjInfo);
%
%
%  Last step: cost to go is zero, target region(s) defined by
%  Tset - Tset should be array of convex polytope objects
%
nTargetRegions = max([1,length(probStruct.Tset)]);
if isfield(probStruct, 'P_N')
    if iscell(probStruct.P_N),
        if length(probStruct.P_N)~=nTargetRegions,
            error('"probStruct.P_N" must have as many elements as "probStruct.Tset"!');
        end
    end
end

Step{horizon+1}.ctrlStruct = cell(1,nTargetRegions);
auxCstruct = cs;
isCostOnTset = isfield(probStruct,'Bi_N');
for ii = 1:nTargetRegions,
    auxCstruct.Pn = probStruct.Tset(ii);
    auxCstruct.Pfinal = probStruct.Tset(ii);
    if ( isCostOnTset ),
        auxCstruct.Bi = {probStruct.Bi_N{ii}};
        auxCstruct.Ci = {probStruct.Ci_N{ii}};
    else
        auxCstruct.Bi = {zeros(1,nx)};
        auxCstruct.Ci = {0};
    end
    Step{horizon+1}.ctrlStruct{ii} = auxCstruct;
end
mergedCstruct    = cs;
mergedCstruct.Pn = probStruct.Tset;
Step{horizon+1}.nPartitions          = nTargetRegions;
Step{horizon+1}.mergedCtrlStruct     = cs; 
Step{horizon+1}.validPartitions      = [1:nTargetRegions];
Step{horizon+1}.adjacencyInfo        = adjInfo;
Pfinal=polytope;
Pn=polytope;
% set this flag if partitions from diferent dynamics can overlap -
% this is used to speed-up intersect and compare procedure
%
if ~isfield(sysStruct,'guardU'),
    isGuardU = false;
else
    isGuardU = ~logical(isempty(sysStruct.guardU)) & ...
        (max(max(cat(1,sysStruct.guardU{:}))) ~= 0);
end

%
% this is for debugging purposes - defines the steps of the DP we
% want to save
saveSteps = [];

% select mplp solver - this is for testing only
%
%mplpSolver = @mpt_mplp_ver4;

%
% first, we assume min-max
%
isMinMax = true;

startTime = cputime;
mplp_runtime = 0;
tsetcon={};
for k = horizon:-1:1,
    %
    merged = {};
    nRegions = length(Step{k+1}.mergedCtrlStruct.Pn);

    roTime = 0;   % time spent in Intersect-And-Compare procedure
    stepMplpRunTime = 0;
    stepStartTime = cputime;
    nConvexPartitions = Step{k+1}.nPartitions;
    nTargets          = nConvexPartitions;
    doMinMax          = true;
    if Options.verbose > -1
        fprintf(1,['*** Step %d, targets partitions: %d, regions in the' ...
                ' previous step %d\n'],k, nConvexPartitions, nRegions);
    end
    %
    Step{k}.ctrlStruct = num2cell(repmat(cs,1,nTargets));
    old_dyn_idx = 1;
    part_idx = 0;      % number of SS partitions in step k
    idxMerged = [];
    idxKickOut = []; % pertitions to kick out after removing overlaps
    skipRemoveOverlaps = ( Options.recedingHorizon & (k > 1) );
    
    min_progress = (horizon - k)/(horizon);
    max_progress = (horizon - k + 1)/(horizon);
    for iDyn = 1:nPWA,
        if Options.verbose > -1,
            disp(['dynamic ' num2str(iDyn)]);
        end
        
        if statusbar,
            dynprogress = (iDyn - 1)/nPWA*0.25;
            prog_min = min_progress + (iDyn - 1) * (max_progress - min_progress)/(nPWA);
            prog_max = min_progress + iDyn * (max_progress - min_progress)/(nPWA);
            if isempty(mpt_statusbar(statbar.handle, dynprogress, prog_min, prog_max)),
                mpt_statusbar;
                error('Break...');
            end     
            prog_cnt = 0;
            prog_total = length(Step{k+1}.validPartitions);
        end
        
        for iTarget = Step{k+1}.validPartitions,
        
            if statusbar,
                targetprogress = prog_cnt/prog_total;
                if ~skipRemoveOverlaps,
                    targetprogress = 0.25*targetprogress;
                end
                prog_cnt = prog_cnt + 1;
                if statusbar,
                    if isempty(mpt_statusbar(statbar.handle, targetprogress, prog_min, prog_max)),
                        mpt_statusbar;
                        error('Break...');
                    end     
                end
            end
            
            if doMinMax,
                %                 tmpProbStruct.Tset =polytope;
                tmpProbStruct.Tset =Step{k+1}.ctrlStruct{iTarget}.Pfinal;
            else
                tmpProbStruc.Tset  = Step{k+1}.mergedCtrlStruct.Pn(iTarget);
            end
            localOptions.pwa_index     = iDyn;
            localOptions.noNoiseOnTset = ~isNoisy;
            mplpStartTime = cputime;
            %
            % transform the original problem into form suitable for mpLP
            %
            if k ~= horizon,
                if isfield(probStruct,'Qy'),
                    tmpProbStruct.P_N = zeros(1,ny);
                else
                    tmpProbStruct.P_N = zeros(1,nx);
                end
            else
                if isfield(probStruct, 'P_N'),
                    if iscell(probStruct.P_N),
                        tmpProbStruct.P_N = probStruct.P_N{iTarget};
                    end
                end
            end

            localOptions.noConstraintReduction = true;
            localOptions.includeLQRset = 0;
            [Matrices] =  mpt_constructMatrices(sysStruct, tmpProbStruct, ...
                localOptions);
            if isinf(-Matrices.W), continue, end;
            %
            % MIN-MAX
            %
            if doMinMax,
                %
                if ( k < horizon ),
                    constraintInfo.idxSysConstr = [1:size(Matrices.G,1)];
                end
                %
                %------------------------------------------------------------------
                %
                % additional constraints: mu >= B_i * x_{k+1} + C_i
                %
                if ( (k < horizon) | isCostOnTset | isNoisy ),
                Bi = cat(1,Step{k+1}.ctrlStruct{iTarget}.Bi{:});
                Ci = cat(1,Step{k+1}.ctrlStruct{iTarget}.Ci{:});
                rowsBiCi  = size(Bi,1);
                tmpZeroG1 = zeros(size(Matrices.G,1),1);
                tmpZeroG2 = zeros(rowsBiCi,size(Matrices.G,2)-nu);
                tmpOnesG  = ones(rowsBiCi,1);
                %
                % add noise offset to "cost-to-go"
                %
                if isNoisy,
                    costOffset = getNoiseOffset([],[],Vnoise,Bi);
                    Ci = Ci + costOffset;
                end
                Ci = Ci + Bi*sysStruct.f{iDyn};
                %
                % modify Matrices
                %
                    % add new variable representing "cost to go"
                    Matrices.H = [Matrices.H 1]; 
                Matrices.G = [Matrices.G tmpZeroG1; ...
                    Bi*sysStruct.B{iDyn}, tmpZeroG2, ...
                    -tmpOnesG];
                Matrices.W = [Matrices.W; ...
                    -Ci];
                Matrices.E = [Matrices.E; ...
                    -Bi*sysStruct.A{iDyn}];
                end
                %
                if ( k < horizon),
                    constraintInfo.idxValueFConstr = ...
                        [length(constraintInfo.idxSysConstr)+1:size(Matrices.G,1)];
                    %
                    % add target set constraints - we'll take care for
                    % the noisy case later
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % This part of code is supposed to include a reduced
                    % set of target set constraints based on adjacency
                    % information from the previous step. However, it's
                    % not working properly due to ... eh, hm,
                    % ... "numerical problems".
                    % Until this is fixed this code should remain
                    % commented.
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % $$$                     nSysVFConstr = size(Matrices.G,1);
                    % $$$                     [Htset,Ktset] = double(Step{k+1}.ctrlStruct{iTarget}.Pfinal);
                    % $$$                     tmpZeroG = zeros(size(Htset,1),size(Matrices.G,2)-nu);
                    % $$$                     Matrices.G = [Matrices.G; Htset * sysStruct.B{iDyn}, tmpZeroG];
                    % $$$                     Matrices.W = [Matrices.W; Ktset - Htset * sysStruct.f{iDyn}];
                    % $$$                     Matrices.E = [Matrices.E; -Htset * sysStruct.A{iDyn}];
                    % $$$                     constraintInfo.idxTSetConstr = [nSysVFConstr+1:size(Matrices.G,1)];
                    constraintInfo.idxTSetConstr = [];
                end
                %
                % add adjacency information from the previous step
                %
                if ( k < horizon ),
                    constraintInfo.adjacencyInfo = Step{k+1}.ctrlStruct{iTarget}.adjacencyInfo;
                    Matrices.constraintInfo = constraintInfo;
                end
                %------------------------------------------------------------------
                %
            else % not MIN-MAX
                %
                % add cost to go of the target region
                %
                Matrices.H(:,1:nu) = Matrices.H(:,1:nu) + ...
                    Step{k+1}.mergedCtrlStruct.Bi{iTarget} * ...
                    sysStruct.B{iDyn};
            end
            %
            if ~isempty(saveSteps) & ~isempty(find(saveSteps==k)),
                fname = sprintf('test_di%d.mat',k);
                fprintf(1,'Saving data to file %s\n',fname);
                save(fname,'Matrices','mplpOptions');
            end
            %
            if ( Options.recedingHorizon & k==1 )
                mplpOptions.skipDualDegenerate = false;
            end
            %
            % solve mplp
            %
            %            try,
            
            if statusbar,
                if isempty(mpt_statusbar(statbar.handle, targetprogress, prog_min, prog_max)),
                    mpt_statusbar;
                    error('Break...');
                end     
            end
            
            [Pn,Fi,Gi,activeConstraints,Pfinal,details] = mpt_mplp(Matrices,mplpOptions);
            
            if statusbar,
                if isempty(mpt_statusbar(statbar.handle, targetprogress, prog_min, prog_max)),
                    mpt_statusbar;
                    error('Break...');
                end     
            end
            
            %            catch,
            %		rethrow(lasterror);
            %	    end
            %
            mplpEndTime = cputime;
            mplpTime  = mplpEndTime - mplpStartTime;
            mplp_runtime = mplp_runtime + mplpTime;
            stepMplpRunTime = stepMplpRunTime + mplpTime;
            if ( ~details.feasible ),
                if Options.verbose > 1
                    fprintf('Infeasible transition\n');
                end
                continue;
            end

            %           figure, plot(Pn);
            % ----------------------------------------------------
            %
            % check if newly generated partition can be discarded,
            % or if it discards some old one
            %
            if ( ~isfulldim(Pfinal) ),
                if Options.verbose > 1,
                    disp(['Discarding lower-dimensional' ...
                        ' partition.']);
                end
                continue;
            end
            if ( part_idx >= old_dyn_idx ),
                [Hfin,kfin] = double(Pfinal);
                newJStarH = [cat(1,details.Bi{:}),-ones(length(Fi),1);...
                    Hfin, zeros(size(Hfin,1),1)];
                newJStarK = [-cat(1,details.Ci{:});kfin];
                newJStar  = polytope(newJStarH, newJStarK);
                if ( dimension(newJStar) < (nx+1) ),
                    %
                    % it can happen
                    %
                    continue;
                end
                skipThisOne = false;
                isStored    = false;
                discarded_parts = [];
                validRange = intersect([old_dyn_idx:part_idx], ...
                    Step{k}.validPartitions);
                
                for valid_idx = validRange,
                    
                    if statusbar,
                        if isempty(mpt_statusbar(statbar.handle, targetprogress, prog_min, prog_max)),
                            mpt_statusbar;
                            error('Break...');
                        end     
                    end
                    
                    
                    if Pfinal >= Step{k}.ctrlStruct{valid_idx}.Pfinal,
                        [Hfin,kfin] = double(Step{k}.ctrlStruct{valid_idx}.Pfinal);
                        Bi = cat(1,Step{k}.ctrlStruct{valid_idx}.Bi{:});
                        Ci = cat(1,Step{k}.ctrlStruct{valid_idx}.Ci{:});
                        oldJStarH = [Bi -ones(size(Bi,1),1);...
                            Hfin, zeros(size(Hfin,1),1)];
                        oldJStarK = [-Ci;kfin];
                        someOldJStar = polytope(oldJStarH,oldJStarK,0,1);
                        if newJStar > someOldJStar, % superset
                            if Options.verbose > 1,
                                disp(''); disp(['Discarding old' ...
                                    ' partition.']);
                            end
                            if ~isStored,
                                Step{k}.ctrlStruct{valid_idx}.Bi     = details.Bi;
                                Step{k}.ctrlStruct{valid_idx}.Ci     = details.Ci;
                                Step{k}.ctrlStruct{valid_idx}.Fi     = Fi;
                                Step{k}.ctrlStruct{valid_idx}.Gi     = Gi;
                                Step{k}.ctrlStruct{valid_idx}.Pn     = Pn;
                                Step{k}.ctrlStruct{valid_idx}.Pfinal = Pfinal;
                                Step{k}.ctrlStruct{valid_idx}.Ai = ...
                                    cell(1, length(Pn));
                                [Step{k}.ctrlStruct{valid_idx}.Ai{:}] = ...
                                    deal(zeros(nx));
                                Step{k}.ctrlStruct{valid_idx}.dynamics = ...
                                    repmat(iDyn,1,length(Pn));
                                Step{k}.ctrlStruct{valid_idx}.adjacencyInfo=details.adjacencyInfo;

                                isStored = true;
                            else
                                discarded_parts(end+1) = valid_idx;
                            end
                        end
                    elseif Pfinal <= Step{k}.ctrlStruct{valid_idx}.Pfinal,
                        Bi = cat(1,Step{k}.ctrlStruct{valid_idx}.Bi{:});
                        Ci = cat(1,Step{k}.ctrlStruct{valid_idx}.Ci{:});
                        someOldJStar = polytope([Bi -ones(size(Bi,1),1)],-Ci,0,1);
                        %
                        if newJStar < someOldJStar, % subset
                            skipThisOne = true;
                            break;
                        end
                    end
                end
                %
                if skipThisOne,
                    if Options.verbose > 1,
                        disp('Discarding the fresh partition.');
                    end
                    continue;
                end
                %
                % update indices of valid partitions
                %
                Step{k}.validPartitions = setdiff(Step{k}.validPartitions, ...
                    discarded_parts);
                if isStored,
                    %
                    % we've already stored new partition
                    %
                    continue;
                end
            end
            %
            part_idx = part_idx + 1; % increase partition index

            %
            % store the partition: it would be VERY, VERY nice if
            % we could eliminate some of partitions here when using MIN-MAX
            %
            Step{k}.ctrlStruct{part_idx}    = cs;
            Step{k}.ctrlStruct{part_idx}.Pn = Pn;
            Step{k}.ctrlStruct{part_idx}.Fi = Fi;
            Step{k}.ctrlStruct{part_idx}.Gi = Gi;
            Step{k}.ctrlStruct{part_idx}.Ai = cell(1,length(Pn));
            [Step{k}.ctrlStruct{part_idx}.Ai{:}] = deal(zeros(nx));
            Step{k}.ctrlStruct{part_idx}.Bi = details.Bi;
            Step{k}.ctrlStruct{part_idx}.Ci = details.Ci;
            Step{k}.ctrlStruct{part_idx}.Pfinal = Pfinal;
            Step{k}.ctrlStruct{part_idx}.dynamics = repmat(iDyn,1,length(Pn));
            Step{k}.ctrlStruct{part_idx}.overlaps = 1;
            Step{k}.ctrlStruct{part_idx}.sysStruct = sysStruct;
            Step{k}.ctrlStruct{part_idx}.probStruct = probStruct;
            Step{k}.ctrlStruct{part_idx}.details.mplpRunTime = mplp_runtime;
            Step{k}.ctrlStruct{part_idx}.adjacencyInfo=details.adjacencyInfo;
            % store partition index in the array of 'valid indices'
            %
            Step{k}.validPartitions(end+1) = part_idx;
            
            if statusbar,
                targetprogress = prog_cnt/prog_total;
                prog_cnt = prog_cnt + 1;
                if statusbar,
                    if isempty(mpt_statusbar(statbar.handle, targetprogress, prog_min, prog_max)),
                        mpt_statusbar;
                        error('Break...');
                    end     
                end
            end

        end
        
        if statusbar,
            if skipRemoveOverlaps,
                targetprogress = 1;
            end
            if isempty(mpt_statusbar(statbar.handle, targetprogress, prog_min, prog_max)),
                mpt_statusbar;
                error('Break...');
            end     
        end
        
        % number of partitions doesn't necesserily corespond to
        % length of partition array (Step{}.ctrlStruct) since some
        % partitions could be empty.
        %
        % check if partitions from different dynamics can overlap -
        % if not, merge partitions for each dynamic separately
        %
        %disp(' ');
        if ( skipRemoveOverlaps ),
            continue;
        end
        if ~isGuardU & nPWA>1,
            roStartTime = cputime;
            mergeRange = intersect([old_dyn_idx:part_idx],...
                Step{k}.validPartitions);
            if isempty(mergeRange),
                % no regions for this dynamic
                continue;
            end
            if statusbar,
                roOptions.statusbar = 1;
                roOptions.status_min = prog_min + 0.25*(prog_max-prog_min);
                roOptions.status_max = prog_max;
                roOptions.status_handle = statbar.handle;
                roOptions.closestatbar = 0;
            end
            merged{iDyn} = mpt_removeOverlaps(Step{k}.ctrlStruct(mergeRange), ...
                roOptions);
            
            % chek which partitions have been kept in the solution and
            % remove the rest
            %
            if ( length(merged{iDyn}.details.keptParts) < length(mergeRange) ),
                %
                % some partitions don't enter the solution
                %
                idxKept      = merged{iDyn}.details.keptParts;
                nKickedOut   = length(mergeRange) - length(mergeRange);
                kickOutParts = setdiff(mergeRange,mergeRange(idxKept));
                idxKickOut   = [idxKickOut kickOutParts];
            end
            idxMerged(end+1) = iDyn;
            old_dyn_idx = part_idx + 1;
            roEndTime = cputime;
            roTime = roTime + roEndTime - roStartTime;
            Step{k}.mergedCtrlStruct.details.roRunTime = ...
                Step{k}.mergedCtrlStruct.details.roRunTime + roTime;
        end
    end
    if ( length(idxKickOut) ),
        idxValid = Step{k}.validPartitions;
        idxValid = setdiff(idxValid,idxKickOut);
        Step{k}.validPartitions = idxValid;
    end
    Step{k}.nPartitions = length(Step{k}.validPartitions);

    %---------------------------------
    % INTERSECT & COMPARE procedure
    %---------------------------------
    %
    % If we are looking for RH solution, we can skip removing
    % overlaps for all steps of the DP, but the last one.
    %
    
    if statusbar,
        targetprogress = 1;
        if statusbar,
            if isempty(mpt_statusbar(statbar.handle, targetprogress, min_progress, max_progress)),
                mpt_statusbar;
                error('Break...');
            end     
        end
    end
    
    if ( ~skipRemoveOverlaps ),
        if statusbar,
            roOptions.statusbar = 1;
            roOptions.status_min = min_progress;
            roOptions.status_max = max_progress;
            roOptions.status_handle = statbar.handle;
            roOptions.closestatbar = 0;
        end
        roStartTime = cputime;
        if ~isGuardU & (nPWA > 1),
            % partitions from different dynamics cannot overlap - just
            % glue them together
            %
            if ~isempty(idxMerged),
                Step{k}.mergedCtrlStruct = cs;
            end
            for merged_idx = idxMerged,
                Step{k}.mergedCtrlStruct.Pn = [Step{k}.mergedCtrlStruct.Pn, ...
                    merged{merged_idx}.Pn];
                Step{k}.mergedCtrlStruct.Pfinal = [Step{k}.mergedCtrlStruct.Pfinal, ...
                    merged{merged_idx}.Pfinal];
                Step{k}.mergedCtrlStruct.Ai = {Step{k}.mergedCtrlStruct.Ai{:}, ...
                    merged{merged_idx}.Ai{:}};
                Step{k}.mergedCtrlStruct.Bi = {Step{k}.mergedCtrlStruct.Bi{:}, ...
                    merged{merged_idx}.Bi{:}};
                Step{k}.mergedCtrlStruct.Ci = {Step{k}.mergedCtrlStruct.Ci{:}, ...
                    merged{merged_idx}.Ci{:}}; ...
                    Step{k}.mergedCtrlStruct.Fi = {Step{k}.mergedCtrlStruct.Fi{:}, ...
                    merged{merged_idx}.Fi{:}};
                Step{k}.mergedCtrlStruct.Gi = {Step{k}.mergedCtrlStruct.Gi{:}, ...
                    merged{merged_idx}.Gi{:}};
                Step{k}.mergedCtrlStruct.dynamics = [Step{k}.mergedCtrlStruct.dynamics, ...
                    merged{merged_idx}.dynamics];
            end
            clear merged;
        elseif nPWA > 1,
            % do the thorough intersect and compare procedure
            %
            Step{k}.mergedCtrlStruct = ...
                mpt_removeOverlaps(Step{k}.ctrlStruct(Step{k}.validPartitions),...
                roOptions);
        elseif nPWA == 1,
            if isempty(Step{k}.ctrlStruct),
                error('Problem is infeasible!');
            end
            Step{k}.mergedCtrlStruct = Step{k}.ctrlStruct{1};
        end
        roEndTime = cputime;
        Step{k}.mergedCtrlStruct.details.roRunTime   = roEndTime - ...
            roStartTime;
        roTime = roTime + Step{k}.mergedCtrlStruct.details.roRunTime;
    end
    stepEndTime = cputime; % just for the clarity of the code
    Step{k}.mergedCtrlStruct.details.mplpRunTime = stepMplpRunTime;
    Step{k}.mergedCtrlStruct.details.runTime     = stepEndTime - stepStartTime;
    
    if Options.statusbar,
        progress = (horizon - k + 1)/horizon;
        if isempty(statusbar(progress, statb))
            error('Break');
        end
    end
end

if statusbar,
    if isempty(mpt_statusbar(statbar.handle, 1)),
        mpt_statusbar;
        error('Break...');
    end     
end

endTime = cputime;

% fill up ctrlStruct
%
ctrlStruct = cs;
ctrlStruct = Step{1}.mergedCtrlStruct;
ctrlStruct.details.runTime = endTime - startTime;
ctrlStruct.details.roRunTime = roTime;
for k = 1:horizon,
    ctrlStruct.details.Horizon{k} = ...
        Step{horizon-k+1}.mergedCtrlStruct;
    if ~isfield(ctrlStruct.details.Horizon{k}, 'Pfinal'),
        ctrlStruct.details.Horizon{k}.Pfinal = polytope;
    elseif ~isa(ctrlStruct.details.Horizon{k}.Pfinal, 'polytope'),
        ctrlStruct.details.Horizon{k}.Pfinal = polytope;
    end
    for pfin_idx = 1:Step{horizon-k+1}.nPartitions,
        ctrlStruct.details.Horizon{k}.Pfinal = ...
            [ctrlStruct.details.Horizon{k}.Pfinal, ...
            Step{horizon-k+1}.ctrlStruct{pfin_idx}.Pfinal];
    end
end
ctrlStruct.Pfinal = ctrlStruct.details.Horizon{horizon}.Pfinal;

if Options.details<1,
    % remove open loop solution if details==0
    details = ctrlStruct.details;
    details = rmfield(details,'Horizon');
    ctrlStruct.details = details;
end

if nPWA==1,
    % indicate that there are no overlaps if we have just one dynamics
    ctrlStruct.overlaps = 0;
end

if closestatbar,
    mpt_statusbar;
end


%
% getUncertaintyOffset
%
function [offset] = getNoiseOffset(noiseH,noiseK,Wvert,Bi)
%
% calculates the uncertainty offset for affine "cost to go"
%   J{k+1} = Bi * x{k+1} + Ci
% The noise is bounded by polytope: noiseH * x <= noiseK
% Wvert is the matrix of the noise polytope vertices.
%
if ~isempty(Wvert),
    offset = max(Bi * Wvert',[],2);
else
    offset = zeros(size(Bi,1),1);
    %
    for ii = 1:length(offset),
        [xopt,fval,lambda,exitflag,how]= ...
            mpt_solveLP(-Bi(ii,:),noiseH,noiseK,[],[],[]);
        if ~strcmp(how,'ok'),
            error(['Something''s wrong! Cannot calculate uncertainty '...
                'offset for "cost to go"']);
        else
            offset(ii) = -fval;
        end
    end
end
