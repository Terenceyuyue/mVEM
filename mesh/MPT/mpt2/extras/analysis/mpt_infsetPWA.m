function [Pn,dynamics,invCtrl]=mpt_infsetPWA(Pn,A,f,Wnoise,Options)
%MPT_INFSETPWA Computes (robust) positive invariant subset for PWA systems
%
% [Pn,dynamics] = mpt_infsetPWA(ctrl)
% [Pn,dynamics] = mpt_infsetPWA(ctrl,Options)
% [Pn,dynamics] = mpt_infsetPWA(Pn,A,f,Wnoise,Options)
% [Pn,dynamics,invCtrlStruct] = mpt_infsetPWA(ctrlStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the maximal (robust) positive invariant set of a PWA System
% x(k+1)=A{i}x(k)+f{i}+w    for x(k) \in Pn(i), w \in Wnoise
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn       - Polytope array defining the area where the PWA system is defined
% A        - Cell array containing dynamic matrices A_i
% f        - Cell array containing dynamic matrices f_i
% Wnoise   - Polytope which bounds additive uncertainty w \in Wnoise (can be empty)
% ctrl     - Explicit controller (MPTCTRL object)
% Options.verbose   - Level of verbosity 0,1 or 2
% Options.nohull    - If set to 1, do not compute convex unions
% Options.maxIter   - maximum number of iterations. Set is not invariant if
%                     iteration is aborted prior to convergence  (default is 200)
% Options.useTmap   - If set to true (default is false), transition map will be
%                     computed to rule out certain transitions
% Options.sphratio  - Gives factor which governs maximum number of separating
%                     hyperplanes computed in transition maps. Number of
%                     separating  hyperplnaes computed at each step is given by
%                     length(Pn)*length(targetPn) / Options.ratio
%                     Default value is 20.
%                     Set this option to 0 if you don't want to impose any limit
%                     on number of separating hyperplanes.
% Options.mergefinal - If set to true (default), tries to simplify
%                      final result by merging regions
% Options.simplify_target - If set to true, tries to merge regions between
%                           iteration steps, after Minkowski calculations, 
%                           and after detecting a non-convex union of more than
%                           two polytopes in iterations. Default value is set to
%                           TRUE if system is subject to noise, and to FALSE
%                           otherwise.
% Options.simplify_method - valid only if .simplify_target set; 'greedy'
%                           (default) - uses greedy merging, other string -
%                           uses optimal merging
%
% NOTE: Length of Pn, A, f must be identical
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Pn        - Polytope array defining the (robust) positive invariant set
% dynamics  - Integer array defining the active dynamics A_i,f_i for each
%             polytope Pn(i)
% invCtrl   - optional; returns a controller which contains control
%                 laws associated to the positive invariant set
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% "Computation of Invariant Sets for Piecewise Affine Discrete Time Systems 
%  subject to Bounded Disturbances"
% S. Rakovic, P. Grieder, M. Kvasnica, D. Q. Mayne and M. Morari, 2003, submitted
% check http://control.ee.ethz.ch for latest info
%
% see also MPT_INFSET
%

% Copyright is with the following author(s):
%
% (C) 2005 Mario Vasak, FER, Zagreb
%          mario.vasak@fer.hr 
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2005 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
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

error(nargchk(1,5,nargin));

if nargout>2 & nargin>=3,
    error('Please provide MPTCTRL object as an input!');
end

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<4         % | isempty(Wnoise)) M.V. later if we find out that the noise is hidden in ctrlStruct, we override this
    Wnoise=polytope;
end

nargs = nargin;
if nargin==2
    if isa(Pn,'polytope') & nargin<3 %M.V. to fix the issue with giving a PWL model (i.e. only Pn and A)
        f=cell(length(A),1);
        Options={};
    elseif isstruct(A), % if we give ctrlStruct and Options, again everything is handeled
        Options = A;
        nargs = 1;
    end
    nargs = 1;
end

if (nargin>2 & nargin<5) | nargin==1
    Options={};
end

if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options,'nounion')
    % if regions should be merged or not
    Options.nounion=0;
end
if ~isfield(Options, 'useTmap'),
    % if set to 1, transition map will be computed to rule out certain
    % transitions
    Options.useTmap = 0;
end
if ~isfield(Options, 'maxIter'),
    % Set is not invariant if iteration is aborted prior to convergence
    Options.maxIter = 200;
end
if ~isfield(Options, 'mergefinal'),
    % If true, tries to merge final result
    mergefinal = 1;
else
    mergefinal = Options.mergefinal;
end
if ~isfield(Options, 'sphratio'),
    sphratio = 20;
else
    sphratio = Options.sphratio;
end
if sphratio == 0,
    % to avoid "division by zero" warnings
    sphratio = 1e-9;
end
if ~isfield(Options,'simplify_method')
    Options.simplify_method = 'greedy';
end

maxIter = Options.maxIter;

if nargs==1,
    ctrl = Pn;
    if isa(ctrl,'polytope')
        error('You should provide additional fields defining autonomous dynamics of a PWA system!')
    end
    if isa(ctrl, 'mptctrl')
        if ~isexplicit(ctrl),
            error('This function supports only explicit controllers!');
        end
        ctrlStruct = struct(ctrl);
        input_is_mptctrl = 1;
    else
        ctrlStruct = ctrl;
        input_is_mptctrl = 0;
    end

    if ~mpt_isValidCS(ctrlStruct)
        error('mpt_infsetPWA: First argument has to be a valid controller structure! See mpt_control for details.');
    end
    sysStruct = ctrlStruct.sysStruct;
    if mpt_isnoise(sysStruct.noise)
        Wnoise = sysStruct.noise;
    end
    Pn = ctrlStruct.Pn;
    Fi = ctrlStruct.Fi;
    Gi = ctrlStruct.Gi;
    [nx, nu] = mpt_sysStructInfo(ctrlStruct.sysStruct);
    if (iscell(sysStruct.A))
        pwasystem=1;
        disp('In PWA case, each dynamics has to be associated with exactly one region!');
        disp('Linking dynamics to regions...');
        Acell={};
        Fcell={};
        for ii=1:length(Pn)
            [x,R] = chebyball(Pn(ii));            % compute center of the chebyshev's ball
            for jj=1:length(sysStruct.A)          % go through all dynamics description
                if max(sysStruct.guardX{jj}*x + ...c
                        sysStruct.guardU{jj}*(Fi{ii}(1:nu, :)*x+Gi{ii}(1:nu))-sysStruct.guardC{jj}) ...
                        <Options.abs_tol,    % check which dynamics is active in the region
                    Acell{ii}=sysStruct.A{jj}+sysStruct.B{jj}*ctrlStruct.Fi{ii}(1:nu, :);
                    Fcell{ii}=sysStruct.B{jj}*ctrlStruct.Gi{ii}(1:nu)+sysStruct.f{jj};
                end
            end
        end
        A = Acell;
        f = Fcell;
    else
        error('system is not PWA!');
    end
else
    if(nargin<3 | isempty(f))
        nx=size(A{1},1);
        for i=1:length(A)
            f{i}=zeros(nx,1);
        end     
    end
end

nx=size(A{1},1);
origin = zeros(nx,1);
emptypoly=polytope;


%%%NOW COMPUTE INVARIANT SUBSET
iter=1;
dynamics=1:length(Pn);
dynamicsX=dynamics; %M.V. according to Rakovic, Grieder et al., see Algorithm 4.1 and relation (9)
PnX=Pn;             %M.V. according to Rakovic, Grieder et al.
notConverged=1;

mergeOpt.verbose=0;
mergeOpt.greedy = strcmpi(Options.simplify_method, 'greedy'); 

if ~isfield(Options,'simplify_target') 
    Options.simplify_target = mpt_isnoise(Wnoise);
end

if Options.simplify_target,
    P_merg=merge(Pn,mergeOpt);
else
    P_merg=Pn;
end

if mpt_isnoise(Wnoise),
    targetPn=P_merg-Wnoise;
    if Options.simplify_target
        targetPn=merge(targetPn,mergeOpt);
    end
else
    targetPn=P_merg;
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

progress = 0;
if statusbar,
    if ~isfield(Options, 'status_handle')
        Options.status_handle = mpt_statusbar('Computing invariant set...');
        closestatbar = 1;
    end
end

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, progress, Options.status_min, Options.status_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

empty_dynamics=[]; %M.V.

while(notConverged>0 & iter<maxIter)
    
    if statusbar,
        prog_min = mod((iter-1), 10)/10;
        prog_max = mod(iter, 10)/10;
        if isempty(mpt_statusbar(Options.status_handle, prog_min, prog_min, prog_max)),
            mpt_statusbar;
            error('Break...');
        end     
    end
    Pn_past=Pn;             %M.V. according to Rakovic, Grieder et al. (need this for convergence detection)
    dynamics_past=dynamics; %M.V. according to Rakovic, Grieder et al. (need this for convergence detection)
    
                     
    dynamics=setdiff(dynamicsX,empty_dynamics);     %M.V. according to Rakovic, Grieder et al.
    Pn=PnX(dynamics);                               %M.V. according to Rakovic, Grieder et al., we always do the transitions from X
                                                    %we remove dynamics
                                                    %that do not yield any
                                                    %transition
    
    if Options.verbose>-1,
        disp(['Iteration Number ' num2str(iter)])
    end
    
    transP=polytope;    %initialize to empty polytope
    tdyn=[];
    notConverged=0;

    lenPn = length(Pn);
    isfulldimPn = zeros(1, lenPn);
    for ii = 1:lenPn,
        if isfulldim(Pn(ii)),
            isfulldimPn(ii) = 1;
        end
    end
    lenTargetPn = length(targetPn);
    isfulldimTarget = zeros(1, lenTargetPn);
    for ii = 1:lenTargetPn,
        if isfulldim(targetPn(ii)),
            isfulldimTarget(ii) = 1;
        end
    end

    if Options.useTmap,
        % compute transition map
        
        % prepare Acell and Fcell
        Acell = cell(1, lenPn);
        Fcell = cell(1, lenPn);
        for ii = 1:lenPn,
            Acell{ii} = A{dynamics(ii)};
            Fcell{ii} = f{dynamics(ii)};
        end
        % compute the transition map
        if Options.verbose > -1,
            fprintf('Computing transition map...\n');
        end
        tmapOptions.targetPn = targetPn;
        tmapOptions.maxsph = ceil(lenPn*lenTargetPn/sphratio);
        tmap = mpt_transmap(Pn, Acell, Fcell, tmapOptions);
        
        if Options.verbose > 1,
            fprintf('Transition map discarded %.2f%% of possible transitions.\n', 100*(1 - nnz(tmap)/numel(tmap)));
        end
    end
    
    for i=1:lenPn
        notConverged_before=notConverged;   %M.V.
        if statusbar,
            if isempty(mpt_statusbar(Options.status_handle, (i-1)/length(Pn), prog_min, prog_max)),
                mpt_statusbar;
                error('Break...');
            end     
        end
        if ~isfulldimPn(i),
            continue
        end
        trans=0;        %initialize transition counter
        tP=emptypoly;   %initialize transition polyarray
        convCtr=0;      %this extra counter is needed for error checks
        for j=1:lenTargetPn
            if ~isfulldimTarget(j),
                continue
            end
            if Options.useTmap,
                possible_transition = (tmap(i, j) == 1);
            else
                possible_transition = 1;
            end
            if (possible_transition), 
                % possible transition exists
                [Px,dummy,feasible]=domain(targetPn(j),...
                    A{dynamics(i)},f{dynamics(i)},Pn(i), 1, Options); %compute set of states Pn(i)->targetPn(j)
                if(feasible) 
                    %transition exists
                    %if Px~=Pn(i), %M.V. we don't use this any more, since
                    %Pn is the polytope of initial X
                    convCtr=convCtr+1;
                    notConverged=notConverged+1;
                    %end
                    trans=trans+1;          %one more transition found
                    tP=[tP Px];             %store transition polytope
                end
            end
        end
        
        if statusbar,
            if isempty(mpt_statusbar(Options.status_handle, (i-1)/length(Pn), prog_min, prog_max)),
                mpt_statusbar;
                error('Break...');
            end     
        end
        
        %try to merge regions
        if(trans>0)

            if Options.verbose>1,
                fprintf('iteration=%d region=%d transitions=%d\n', iter, i, trans);
            end
            
            P_to_compare_candidates=Pn_past(find(dynamics_past==dynamics(i)));
%             P_to_compare=polytope;
%             for kk=1:length(P_to_compare_candidates)
%                 if dointersect(P_to_compare_candidates(kk),tP)
%                     P_to_compare=P_to_compare_candidates(kk);
%                     break;
%                 end
%             end
            
            if ~isfulldim(P_to_compare_candidates) %debug check - this shouldn't happen
                disp('Strange: \Omega_{k+1}\varsubsetneq\Omega_k');
                disp('Discarding this transition; assuming numerical difficulties.');
                notConverged=notConverged_before;
                continue;
            end
            
            if Options.nounion==1,
                how=0;
            else
                [Pu,how]=union(tP,Options);
                if how,
                    % union is convex
                    %PuPn_equal = (Pu == Pn(i));
                    PuPn_equal = (length(P_to_compare_candidates)==1) & (P_to_compare_candidates<=Pu);
                    %M.V. according to Rakovic, Grieder et al., this implies that P_to_compare_candidates==Pu,
                    %since the other direction P_to_compare_candidates>=Pu
                    %we know holds for this iterative way of constructing
                    %invariant set (see Rakovic,Grieder, page 3:
                    %\Omega_{k+1}\subseteq\Omega_k, and thus in each
                    %dynamics i we have
                    %\Omega_{k+1}\cup\Q_i^*\subseteq\Omega_{k}\cup\Q_i^*)
                end
            end
            if(how)
                %union is convex (or simplified using greedy merging) => overwrite old set
                %if trans>1 & PuPn_equal,
                if PuPn_equal   %allways if PuPn_equal do the subtarction, trans==1 is not a special case any more
                    notConverged=notConverged-convCtr;%union is identical to original set; reduce transition counter again
                end
                transP = [transP Pu];
                for kk=1:length(Pu),
                    tdyn(end+1)=dynamics(i);
                end
            else
                %union is not convex
                if Options.simplify_target & trans>=3      %M.V. issue b)
                    tP=merge(tP,mergeOpt);
                end
                
                if P_to_compare_candidates<=tP
                    %M.V. according to Rakovic, Grieder et al., this implies that P_to_compare_candidates==Pu,
                    %since the direction P_to_compare_candidates>=Pu
                    %we know holds for this iterative way of constructing
                    %invariant set (see Rakovic,Grieder, page 3: \Omega_{k+1}\subseteq\Omega_k,
                    %and thus in each dynamics i we have
                    %\Omega_{k+1}\cup\Q_i^*\subseteq\Omega_{k}\cup\Q_i^*)
                    notConverged=notConverged-convCtr;
                    if length(P_to_compare_candidates)<length(tP)
                        tP=P_to_compare_candidates;
                    end
                end

                transP = [transP tP];  
                for kk=1:length(tP)
                    tdyn(end+1)=dynamics(i);             %dynamics of new polytope are same as original polytope at t-1
                end
            end
        else
            empty_dynamics(end+1)=dynamics(i);
%            P_to_compare_candidates=Pn_past(find(dynamics_past==dynamicsX(i)));
%            if isfulldim(P_to_compare_candidates)
            notConverged=notConverged+1;        %M.V. if there are no transitions from some dynamics, and in the previous step there were some transitions, then convergence did not occur
 %           end
        end
    end
    
    Pn=transP;      %write results for next iteration (M.V. for convergence checking)
    dynamics=tdyn;

    
    if ~notConverged | iter==maxIter  %M.V. so as not to calculate Minkdiff when convergence is detected
        break;
    end

    if ~isfulldim(Pn)
        break;     %M.V. if there are no feasible transitions, invariant set is empty!
    end
    
    if Options.simplify_target %M.V. issue b)
        P_merg=merge(Pn,mergeOpt);
    else
        P_merg=Pn;
    end

	if mpt_isnoise(Wnoise),
        targetPn=P_merg-Wnoise;
        if Options.simplify_target
            targetPn=merge(targetPn,mergeOpt);
        end
    else
        targetPn=P_merg;
    end
    
    iter=iter+1;
end


if(iter>=maxIter)
    if closestatbar,
        mpt_statusbar;
    end
    error('mpt_infsetPWA: Computation aborted because maximum number of iterations was reached');
end


if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, 1)),
        mpt_statusbar;
        error('Break...');
    end     
end

%%TRY TO MERGE THE FINAL RESULT
if mergefinal,
    for dyn=1:max(dynamics)
        ctr=0;
        Pt=emptypoly;
        
        dynSet=find(dynamics==dyn); 
        Pt=Pn(dynSet);   %set of polytopes with dynamic dyn
        
        %try to merge polytopes
        [Pu,how]=union(Pt, Options); 
        
        if(how==1)
            %union is convex
            Pn = [Pn Pu];               %add new set 
            for kk=1:length(Pu)
                dynamics(end+1)=dyn;    %add new set
            end
            
            Pn(dynSet)=[];          %remove old sets 
            dynamics(dynSet)=[];    %remove old sets
        end
    end
end

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, (i-1)/length(Pn), prog_min, prog_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

dynamics = dynamics(:)';

if nargout>2,
    invCtrlStruct = ctrlStruct;
    invCtrlStruct.Pn = Pn;
    invCtrlStruct.Pfinal = Pn;
    invCtrlStruct.Fi = {ctrlStruct.Fi{dynamics}};
    invCtrlStruct.Gi = {ctrlStruct.Gi{dynamics}};
    invCtrlStruct.Ai = {ctrlStruct.Ai{dynamics}};
    invCtrlStruct.Bi = {ctrlStruct.Bi{dynamics}};
    invCtrlStruct.Ci = {ctrlStruct.Ci{dynamics}};
    invCtrlStruct.dynamics = ctrlStruct.dynamics(dynamics);
    invCtrlStruct.details.isinvariant = 1;
    invCtrl = mptctrl(invCtrlStruct);
end

if closestatbar,
    mpt_statusbar;
end

return
