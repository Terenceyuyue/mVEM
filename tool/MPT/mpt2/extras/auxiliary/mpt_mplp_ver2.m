function [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver2(Matrices,Options)
%MPT_MPLP Explicitly solves the given linear program (LP)
%
% [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp(Matrices,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Multiparametric linear programming
%
% Solves the problem
%   V(x) = min H U
%           U
%   s.t.   G U <= W + E x
%          bndA*x <= bndb
%
% As a solution we get 'nR' regions
%   Pn(i)={x : H x <= K}
% 
% with the optimal control law
%   U = Fi{i} x + Gi{i}
%
% and the corresponding cost function expression
%   V(x) = Bi{i} x + Ci{i}
%  
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Matrices - a struct with all the parameters which are needed.
%            See description above for explanation.
%   Matrices.H=H;
%   Matrices.G=G;   
%   Matrices.W=W;
%   Matrices.E=E;
%   Matrices.bndA=bndA;   Limits on exploration space, i.e. bndA*x<=bndb
%   Matrices.bndb=bndb;
%
% Options.verbose      - level of verbosity
% Options.lpsolver     - which LP solver to use (help mpt_solveLP)
% Options.max_iter     - maximum number of iterations of the algorithm
% Options.step_size    - length of step over a facet
% Options.f_perturb    - Perturbation of the optimization direction
% Options.nu           - How many elements to extract from the optimizer (to
%                        deal with slacks) 
% Options.debug_level  
%           Due to numerical problems tiny regions are sometimes difficult to
%           calculate, i.e. are not identified at all. This may create "gaps"
%           in the computed control law. For the exploration, these will be
%           jumped over and the exploration in the state space will continue.
%           "debug_level" can have three values:
%       
%           0: No debug done
%           1: A tolerance is given to find gap in the region partition,
%              small empty regions inside the region partition will be discarded.
%              Note that this is generally not a problem, since the feedback law 
%              is continuous and can therefore be interpolated easily.
%              Correction to the calculation of the outer hull.
%           2: Zero tolerance to find gap in the region partition, empty regions
%              if they exist, will be detected, i.e. the user will be notified.
%              Correction to the calculation of the outer hull.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Pn,Fi,Gi           - for region Pn(i).H*x <= Pn(i).K computed input is
%                      U=Fi{i}*x+Gi{i}   
% activeConstraints  - Cell Array which stores the active constraints 
%                      of the optimizer in each region.
% Phard              - Defines the feasible state space partition (i.e. union of
%                      all regions) as Phard.H*x<=Phard.K
% details            - a structure with the following fields:
%     nRegions  number of regions
%     Pn        polyhedral partition
%     Fi        control law
%     Gi        control law
%     BC        connection list
%     Bi        value function
%     Ci        value function
%     nHard     number of hard constraints
%     Phard     polytope given by hard constraints
%     nb        number of constraints for each region
%     LISTa     list of active constraints
%
% see also MPT_CONSTRUCTMATRICES, MPT_MPQP, MPT_OPTCONTROL, MPT_OPTCONTROLPWA

% Copyright is with the following author(s):
%
% (C) 2004 Miroslav Baric, Automatic Control Laboratory, ETH Zurich,
%     baric@control.ee.ethz.ch    
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%     baotic@control.ee.ethz.ch
% (C) 2002 Francesco Borrelli, Automatic Control Laboratory, ETH Zurich

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

error(nargchk(1,2,nargin));

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end
if nargin<2,
    Options=[];
end

if ~isfield(Options,'debug_level'),
    Options.debug_level=mptOptions.debug_level;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'rel_tol'),
    Options.rel_tol=mptOptions.rel_tol;
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
% %++++++++++++++++++++>>>>>>>>>>> MODIFIED/CORRECTED/ADDED BY MATO BAOTIC 
if ~isfield(Options,'qpsolver'),
    Options.qpsolver=mptOptions.qpsolver;
end
% %++++++++++++++++++++>>>>>>>>>>> MODIFIED/CORRECTED/ADDED BY MATO BAOTIC 

if ~isfield(Options,'max_iter'),
    Options.max_iter=20;       % maximum number of iterations of the algorithm
end
if ~isfield(Options,'f_perturb'),
    Options.f_perturb=0.5;     % Perturbation of the optimization direction up to 50% 
end
if ~isfield(Options,'step_size'),
    Options.step_size=mptOptions.step_size;
end
if ~isfield(Options,'ispwa'),
    Options.ispwa=0;
end
if ~isfield(Options,'nu'),
    % how many elements to extract from the optimizer (to deal with slacks)
    Options.nu=size(Matrices.G,2);
end

% new options - force "smooth" optimizer, i.e.
% the optimizer without discontinuities between crticial regions
%
if ~isfield(Options,'force_smooth')
    Options.force_smooth = 1;
end

lpsolver      = Options.lpsolver;
EMPTY_ROW_TOL = Options.abs_tol;         % Tolerance for declaring
                                         % that row is empty
CONSTR_TOL    = Options.rel_tol;     % Tolerance for declaring that
                                     % constr. is redundant
                                     
ALPHA         = max(Options.step_size,1e-7);
DEBUG         = (Options.debug_level==2); % DEBUG is true if
                                          % debug_level==2
ALPHAmax      = Options.step_size;        % maximum step                                               
ALPHAit       = 50;                       % number of iterations
ALPHAinc      = ALPHAmax/(ALPHA*ALPHAit);        

L1   = Matrices.H;
G    = Matrices.G;
S    = Matrices.E;
W    = Matrices.W;
bndA = Matrices.bndA;
bndb = Matrices.bndb;

nx = size(S,2);
nu = size(G,2);
nC = size(G,1);

emptypoly=polytope;

nRegions         = 0;       % number of regions
nHard            = 0;       % number of hard constraints
hardA            = [];      % hard constraints
hardb            = [];      % hard constraints
no_of_constr     = [];      % number of constraints for each region
list_active      = {};      % list of active constraints

Pn = emptypoly;
Hn = {};      % normalized polyhedron description
Kn = {};      % normalized polyhedron description
Fi = {};      % control law
Gi = {};      % control law
BC = {};      % connection list
Bi = {};      % value function
Ci = {};      % value function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   %%
%%      FIND THE STARTING POINT      %%
%% (IF ORIGINAL PROBLEM IS FEASIBLE) %%
%%                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if isempty(bndA),
    crA = [G -S];
    crb = W;
else
    crA = [G -S; zeros(size(bndA,1),size(G,2)+size(S,2)-size(bndA,2)) bndA];
    crb = [W; bndb];
end
%
% pick up the initial point as the chebyshev center of the (z,x) polytope
%
ZXpoly = polytope(crA,crb);
[xCheby,rCheby] = chebyball(ZXpoly);
if rCheby <= Options.abs_tol,
%
% no basic solution, i.e. (z,x) polytope doesn't exist
%
    error(['ERROR in MPLP: initial problem is infeasible or unbounded!' ...
           ' Check your constraint and problem matrices.']);
end

% set the center of coordinate system to the projection of
% the Chebyshev center to the parameter space
%
if ~isfield(Options,'center'),
    center = xCheby(nu+1:nu+nx);
end

%================================================%
% Solve optimization problem for point xFeasible %
%================================================%
%
% %++++++++++++++++++++>>>>>>>>>>> MODIFIED/CORRECTED/ADDED BY MATO BAOTIC 
xInit = xCheby(nu+1:nu+nx);
% xInit = [-7.6 3.3]';
%
% search for an initial point which lies inside the full
% dimensional critical region
%
initPointFound = 0;

xPert=zeros(nx,1);
for i = 1:ALPHAit,
    xFeasible = xInit + xPert;
    [actSetFound,idx_active,lambda,zopt,errstr] = findas(L1,G,W,S,[],xFeasible, ...
                                                      Options);
    if actSetFound == 1,
        cr = getCR(G, S, W, bndA, bndb, idx_active, xFeasible, ...
                        zopt, lambda, Options);
        if cr.type ~= 0, % ful dimensional CR
            initPointFound = 1;
            break;
        end
            
        % otherwise, try to perturb the initial point within the
        % borders of the Chebyshev ball
        %
        xPert = 2*rand(nx,1)-1;
        xPert = rand(1) * rCheby * xPert / sqrt(nx);
    end
end
% %++++++++++++++++++++>>>>>>>>>>> MODIFIED/CORRECTED/ADDED BY MATO BAOTIC 


if initPointFound == 0,
    Ri={};
    nRegions=1;
    nHard=nC;
    hardA=bndA;
    hardb=bndb;
    disp('MPLP: No feasible starting point!');
    
    % This is if you whant to store the initial region where either the
    % problem is infeasible or there are only flat CR
    %
    Ri.list_active{nRegions} = list_active;
    Ri.no_of_constr = no_of_constr;
    %Ri.Pn = polytope;
    Ri.Pn = mptOptions.emptypoly;
    Ri.Fi{nRegions}=repmat(inf,nu,nx);
    Ri.Gi{nRegions}=repmat(inf,nu,1);
    Ri.Bi{nRegions}=repmat(inf,nx,1);
    Ri.Ci{nRegions}=Inf;
    Ri.BC{nRegions}=BC;
    Ri.nRegions=0;
    Ri.nHard=nHard;
    Ri.Phard = polytope(hardA, hardb);
    activeConstraints=Ri.list_active;
    Phard=Ri.Phard;
    for i=1:Ri.nRegions
        Ri.Bi{i}=L1(:)'*Ri.Fi{i};
        Ri.Ci{i}=L1(:)'*Ri.Gi{i};
    end
    details=Ri;
    Pn=Ri.Pn;
    Fi=Ri.Fi;
    Gi=Ri.Gi;
    return;
end

nRegions = nRegions + 1; % new region

Pquadrant    = {};
Pquadrant    = cell(1,2^nx); % number of quadrant is 2^nx
q = sub_whichquadrant(cr.P, center);
for qqq=1:length(q),
    Pquadrant{q(qqq)}(end+1) = nRegions;
end

%
% store the initial region and continue exploration of the
% parameter space
%
list_active{end+1} = idx_active;
Fi{end+1}    = cr.Fi;
Gi{end+1}    = cr.Gi;
Pn(nRegions) = cr.P;
no_of_constr(nRegions) = nconstr(Pn(nRegions));
BC{nRegions} = zeros(no_of_constr(nRegions),1);

%%%%%%%%%%%%%%%%%%%%%%%%
%%                    %%
%% DO THE EXPLORATION %%
%%                    %%
%%%%%%%%%%%%%%%%%%%%%%%%
%
region = 1;
while region <= nRegions,
    if Options.verbose>1,
        disp(sprintf('Region: %d/%d, %d borders, %d hard', region,nRegions,no_of_constr(region),nHard));
    elseif Options.verbose==1,
        if mod(region,20)==0,
            disp(sprintf('Region: %d/%d', region,nRegions));
        end
    end

    [borderDirections, Knr] = double( Pn(region) );
    %
    % borderDirections are normalized
    %
    for border = 1:no_of_constr(region),
        %
        % Create a point close to the border
        %
        [xBorder, RBorder] = facetcircle(Pn(region), border, ...
                                         Options);
        xBeyond = xBorder + borderDirections(border,:)' * ALPHA;
        
        %
        % check if the point is outside the feasible region
        %
        if ( any((bndA * xBeyond - bndb) > 0) ),
            % 
            % we're outside the feasible region of parameters
            %
            continue;
        end
        
        %
        % Check if the point is inside of any existing polyhedra
        %
        isinOpt.abs_tol   = 0;
        isinOpt.fastbreak = 0;
        vi = (xBeyond - center) > 0;
        q = 1 + 2.^(0:nx-1) * vi;
        neighbor = [];
        [isin,neighbor] = isinside(Pn(Pquadrant{q}), xBeyond, isinOpt);
        nMatch = length(neighbor);
        if nMatch > 0,
            neighbor = neighbor(end);
            BC{region}(border) = neighbor; 
            if nMatch > 1,
                if ~Options.ispwa & Options.verbose>0,
                    disp('Polyhedra should not overlap');
                end
            end
            %
        elseif nHard==0 | all(hardA*xBeyond <= hardb),
            %==============================================%
            % Solve optimization problem for point xBeyond %
            %==============================================%
            bordvect = borderDirections(border,:)';
            for niter = 1:ALPHAit,
                [found,idx_active,lambda,zopt,errstr] = findas (L1,G,W,S,[], ...
                                                    xBeyond,Options);
                if ( found == 1 ),
                    cr = getCR (G, S, W, bndA, bndb, idx_active, ...
                                   xBeyond, zopt, lambda, Options);
                    if ( cr.type == 1 | cr.type == -1 ) 
                        %
                        % full-dimensional critical region has been
                        % found
                        %
                        break;
                    end
                else
                    break;
                end
                % perturb the initial point a little bit
                %
                rvect = randn(size(bordvect)); 
                dotpr = dot(rvect,bordvect);
                while ( (norm(rvect) < 10*Options.abs_tol) | ...
                        (abs(dotpr)  < 10*Options.abs_tol) ),
                    rvect = randn(size(bordvect)); 
                    dotpr = dot(rvect,bordvect);
                end
                rvect = rvect / norm(rvect);
                pertdir = rvect - dot(rvect,bordvect) * bordvect;
                pertdir = 0.5 * rand(1) * RBorder * pertdir / norm(pertdir);
                xBeyond = xBorder + ALPHA * (pertdir + bordvect);
            end
            %
            if found == -1 %| isequal(idx_active,list_active{region}),
                % 
                % infeasible LP subproblem
                %
                nHard = nHard + 1;
                hardA(nHard,:) = borderDirections(border,:);
                hardb(nHard,:) = Knr(border,:);
                %
            else
                if ( cr.type == 0 )
                    %
                    % flat region
                    %
                    %disp ('Skipping a flat region.');
                    continue;
                elseif ( cr.type == -1 )
                    %
                    % dual degenerate
                    %
                    if Options.verbose>1,
                        disp ('Storing dual degenerate subregion.');
                    end
                end
                %
                % store a new region
                %
                nRegions = nRegions + 1; 
                BC{region}(border) = nRegions;
                list_active{end+1} = idx_active;
                Pn(end+1) = cr.P;
                Fi{end+1} = cr.Fi;
                Gi{end+1} = cr.Gi;
                q = sub_whichquadrant(cr.P,center);
                for qqq = 1:length(q),
                    Pquadrant{q(qqq)}(end+1) = nRegions;
                end
                no_of_constr(end+1) = nconstr(cr.P);
                BC{end+1} = zeros(no_of_constr(end),1);
            end
        end % nMatch==0
    end % border
    region = region + 1;
end % END EXPLORATION

Ri.list_active = list_active;
Ri.no_of_constr=no_of_constr;
Ri.Pn = Pn;
Ri.Fi = Fi;
Ri.Gi = Gi;
Ri.Bi = Bi;
Ri.Ci = Ci;
Ri.BC = BC;
Ri.nRegions = nRegions;
if nHard == 0,
    Ri.nHard = length(bndb);
    if Ri.nHard ~= 0,
        Ri.Phard = polytope(bndA,bndb);
    else
        Ri.Phard = emptypoly;
    end
else
    hardA = [hardA; bndA];
    hardb = [hardb; bndb];
    Ri.Phard = polytope(hardA, hardb);
    Ri.nHard = nconstr(Ri.Phard);
end
Ri.activeConstraints = Ri.list_active;
activeConstraints    = Ri.list_active;
Phard                = Ri.Phard;
for i = 1:Ri.nRegions,
    Ri.Bi{i} = L1(:)' * Ri.Fi{i};
    Ri.Ci{i} = L1(:)' * Ri.Gi{i};
    Ri.Fi{i} = Ri.Fi{i}(1:Options.nu,:);
    Ri.Gi{i} = Ri.Gi{i}(1:Options.nu,:);
end
Fi = Ri.Fi;
Gi = Ri.Gi;
details = Ri;
if Options.verbose>0,
    disp(sprintf('mpt_mplp: %d regions', nRegions));
end

return;


%=============================================================================
%
%function [out,ii,indexmaxdiff]=findas(f,A,b,F,indeq,tq, ...
%                                     indexmaxdiff,Options)
function [out,idx_active,lambda,xoptLP,errstr] = findas(f,A,b,F,indeq,tq,Options)
%
%=============================================================================
%
% Find set of active constraints.
%
% out=1 => tq belong to a non-flat region and idx_active is list of active
% constraints defining a base
% out=-1 => problem infeasible
% f - cost
% A z <= b+F tq
%-----------------------------------------------------------------------------
    
    errstr = [];
    idx_active = [];         % Index of active constraints
    out = 1;                 % out = 1 if the algorithms finds a
                             % basis for for theta=tq.    
  
    indineq = setdiff((1:size(A,1)),indeq);
    nu = length(f);
    Aeq=A(indeq,:);
    Aineq=A(indineq,:);
    bigb2 = b + F*tq; 

    if ~isempty(indeq),
        bigb2eq=bigb2(indeq);
        bigb2ineq=bigb2(indineq);
    else
        bigb2ineq=bigb2;
        bigb2eq=[];
    end

    if ( ~Options.force_smooth ),
        %----------
        % SOLVE LP
        %----------
        [xoptLP,fval,lambda,exitflag,how]=mpt_solveLPi(f(:)',Aineq,bigb2ineq, ...
                                                        Aeq,bigb2eq,[], ...
                                                      Options.lpsolver);
        if ~strcmp(how,'ok'), 
            errorstr = ['MPLP error: cannot find active constraints due' ...
                        ' to the infeasible PRIMAL subproblem.'];
            out = -1;
            return;
        end
        % the problem is feasible so find the the set of active
        % constraints:
        %
        idx_active = find( abs(A*xoptLP-bigb2) < Options.abs_tol ); % Indices of active
                                                                    % constraints
    else
        %
        % formulate and solve dual problem
        %
        % $$$         dualAeq   = A';
        % $$$         dualBeq   = -f';
        % $$$         dualAineq = -eye(size(Aineq,1));
        % $$$         dualBineq = zeros(size(dualAineq,1),1);
        % $$$         dualf     = bigb2;
        % $$$         [lambda,dualfval,dlambda,dualexitflag,dualhow] = ...
        % $$$             mpt_solveLPi(dualf,dualAineq,dualBineq,dualAeq,dualBeq, ...
        % $$$                   [],Options.lpsolver);
        % $$$         if ~strcmp(dualhow,'ok'),
        % $$$             %
        % $$$             % dual problem is infeasible or unbounded for some
        % $$$             % reason
        % $$$             %
        % $$$             errstr = ['MPLP error: cannot find active constraints due to ' ...
        % $$$                       'the infeasible primal subproblem.'];
        % $$$             out=-1;
        % $$$             return;
        % $$$         end
        %
        % solve primal and get indices of active constraints
        %
        [xoptLP,fval,lambda,exitflag,how]=mpt_solveLPi(f(:)',Aineq,bigb2ineq, ...
                        Aeq,bigb2eq,[], Options.lpsolver);
        if ~strcmp(how,'ok'), 
            %
            % dual problem is infeasible or unbounded for some
            % reason 
            %
            errstr = ['MPLP error: cannot find active constraints due to ' ...      
                      'the infeasible primal subproblem.'];
            out = -1;
            return;
        end
        idx_active = find( abs(Aineq*xoptLP-bigb2ineq) < Options.abs_tol );  
    end
    
    return
%
%============= END findas ===========================
    
    
%============= getCR ================================
%
% get H-description of a critical region
function cr  = getCR (G, S, W, bndA, bndb, idx_active, x0, z0, lambda, Options)
%
%====================================================
%    
%
% Compute matrices the solution X and the polyhedra crA crb in this way:
% apply Gauss to the following system 
% Gtilde   -St     |z  =  Wt  
%                  |x  
% to obtain
%      
%  Iz       0    L     |z    =  |Wt 
%  0        Ix   N     |x      
%
    cr = struct('P',polytope,'Fi',[],'Gi',[],'type',0);
    nu = size(G,2); nx = size(S,2);
    
    %
    % remove all active constraints that don't satisfy strict
    % complementarity condition - this way we detect dual
    % degeneracy
    %
    compl_lambdas = find(abs(lambda(idx_active)) > Options.abs_tol);
%
    idx_inactive = 1:size(G,1); 
    idx_inactive(idx_active) = [];
    Ga  = G(idx_active,:);
    Sa  = S(idx_active,:);
    Wa  = W(idx_active,:);
    Gna = G(idx_inactive,:);
    Wna = W(idx_inactive);
    Sna = S(idx_inactive,:);
    %
    if ( length(idx_active) < nu ),
        cr.type = -1;  % dual degeneracy
        %
        %
        % just for the test we're doing the projection of the XZ
        % polytope to X subspace
        %
% $$$ 	rankGa = rank(Ga,Options.abs_tol);
% $$$ 	[Qt,Rt] = qr(Ga');
% $$$ 	P1      = Qt(:,1:size(Rt,2)) / (Rt(1:size(Rt,2),1:size(Rt,2))');
% $$$ 	P2      = null(Ga);
% $$$ 	HaKa    = Gna * P1 * ((Ga*P1) \ ([Sa Wa]));
% $$$ 	Hx = [HaKa(:,1:nx) - Sna;bndA];
% $$$ 	Hz = [Gna*P2;zeros(size(bndA,1),size(P2,2))];
% $$$ 	H = [Hx Hz];
% $$$ 	K = [Wna - HaKa(:,nx+1);bndb];
% $$$ 	XZPoly = polytope(H,K);
% $$$ 	if ~isfulldim(XZPoly),
% $$$ 	    cr.type = 0;
% $$$ 	    disp ('===> The region is flat!');
% $$$ 	    return;
% $$$ 	end
%        cr.P = projection(XZPoly,[1:nx]);
%       return;
        %
    else
        [Q,R] = qr([Ga, -Sa]);
        U = R(1:nu,1:nu); 
        P = R(1:nu,nu+1:nu+nx);
        if ( length(idx_active) > nu ),
            D = R(nu+1:size(R,1),nu+1:size(R,2));    
            %
            if ( rank(D, Options.abs_tol) > 0 ),
                %
                % critical region is lower dimensional facet of another region
                %
                cr.type = 0;
                return;
            end
        end
        %
	% check the rank of matrix and strict complementarity
	%
	if ( length(compl_lambdas) < length(idx_active) ),
	    %
	    % strict complementarity is violated => DUAL DEGENERACY
	    %
	    cr.type = -1;
	elseif ( rank(U,Options.abs_tol) < nu ),
	    %
	    % optimizer not unique => DUAL DEGENERACY
	    %
	    cr.type = -1;
        else
	    % 
	    % we have full-dimensional critical region
	    %
            cr.type = 1;
        end
    end

    if ( cr.type == -1 ),
        %
        %
        %
        %-----------------
        % dual degeneracy
        %-----------------
        % solution proposed in:
        % 
        % J. Spjotvold, P.Tondel, T.A.Johansen
        % "A Method for Obtaining Continuous Solutions to
        % Multiparametric Linear Programs", submitted to CDC04
        %
        % formulate QP which minimizes optimizer norm along the
        % active constraints
        %-------------------------------------------------------
	
	% take only those constraints which satisfy strict
        % complementarity
	%
	idx_active = idx_active(compl_lambdas);
        idx_inactive = setdiff((1:size(G,1)),idx_active);
        bQP = W + S * x0;
        Aqp_ineq = G(idx_inactive,:); 
        Bqp_ineq = bQP(idx_inactive,:);
        Aqp_eq   = G(idx_active,:); 
        Bqp_eq   = bQP(idx_active,:);
        %
        %------------
        % SOLVE QP
        %------------   
        Hess = eye(nu); f = zeros(nu,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Tondel, example
        %
% $$$         Hess = eye(3); f = zeros(3,1);
% $$$         G = [1 0 -1;-1 0 -1;0 1 -1;0 -1 -1];
% $$$         S = [1 0;-1 0;0 -1;0 1];
% $$$         W = -[1 1 1 1]';
% $$$         x0 = [0 0.5]';
% $$$         Aqp_ineq = G; bQP = W + S * x0;
% $$$         Bqp_ineq = bQP; Aqp_eq = []; Bqp_eq = [];
% $$$         bndA = [1 0;-1 0;0 1;0 -1];
% $$$         bndb = [2 2 2 2]';
% $$$         nu = 3;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [zoptQP,lambda,how,eflag,objqp] = ...
            mpt_solveQP (Hess, f, Aqp_ineq, Bqp_ineq, ...
                         Aqp_eq, Bqp_eq, z0, Options.qpsolver);
        idx_activeQP = find ( abs(G*zoptQP-bQP) < Options.abs_tol );
        Ga = G(idx_activeQP,:);
        Wa = W(idx_activeQP,:);
        Sa = S(idx_activeQP,:);
        rankGa = rank(Ga, Options.abs_tol);
        %
        if ( rankGa < length(idx_activeQP) ),
            %
            % primal degenerate QP: we'll use the projection
            % here
            %
            % calculate the optimizer from the reduced set of
            % equations
            %
            
            if Options.verbose>1,
                disp('===> Primal degenerate QP');
            end
            
            % obtain a reduced set of equations satisfying LICQ
            %
            [Qt,Rt,permT] = qr(Ga');
            [Q,R,perm]    = qr(Ga);
            reduxGa = permT' * Ga;
            reduxGa = reduxGa(1:rankGa,:);
            idx_activeQPredux = permT' * idx_activeQP;
            idx_activeQPredux = idx_activeQPredux(1:rankGa);
            %
            % calculate the optimizer
            %
            FiGi = reduxGa' * ((reduxGa * reduxGa') \ [ ...
                S(idx_activeQPredux,:) W(idx_activeQPredux)]);
            Fi = FiGi(:,1:nx); 
            Gi = FiGi(:,nx+1);
            %
            % determine the critical region: the solution is
            % obtained using null-space method. Matrix P1 is chosen
            % such that matrix (Ga'*P1) is well conditioned.
            %
            idx_inactiveQP = setdiff((1:size(G,1)),idx_activeQP);
            [dummy,idx_activediffQP] = setdiff(idx_activeQP, ...
                                               idx_active);
            nactive = length(idx_activeQP);
            Gna = G(idx_inactiveQP,:);
            Sna = S(idx_inactiveQP,:);
            Wna = W(idx_inactiveQP,:);
            P1  = Q(:,1:nu) / (R(1:nu,1:nu)');
            P2  = null(Ga');
            L   = P1 * ((Ga' * P1) \ FiGi); % this defines P1*lambda1
            H1  = [L(idx_activediffQP,1:nx), -P2(idx_activediffQP,:)];
            K1  = -L(idx_activediffQP,nx+1);
            H2  = [Gna * Fi - Sna; bndA];
            K2  = [Wna - Gna * Gi; bndb];
            H   = [H1; H2,zeros(size(H2,1),size(P2,2))];
            K   = [K1;K2];
            LXpoly = polytope(H,K);
            if ~isfulldim(LXpoly), 
                % 
                %  can this really happen? Why?
                %
                cr.type = 0;
                return;
            end
	    %
	    % old version of MPT used this syntax
            % cr.P = projection(LXpoly,[nx+1:size(H,2)]);
	    %
            cr.P = projection(LXpoly,(1:nx),Options);
	    % this should happen only due to the bug in projection
            % algorithm
	    %
            if ~isfulldim(cr.P),
                cr.type = 0;
            end
        else
            %
            % nondegenerate QP
            %
            idx_inactiveQP = setdiff((1:size(G,1)),idx_activeQP);
            %
            % compute the boundaries of the sub- region and the
            % optimizer
            %
            Gna = G(idx_inactiveQP,:);
            Wna = W(idx_inactiveQP,:);
            Sna = S(idx_inactiveQP,:);
            lambdaHK = -(Ga*Ga') \ [Sa Wa];
            Fi = -Ga' * lambdaHK(:,1:nx);
            Gi = -Ga' * lambdaHK(:,nx+1);
            [dummy,idx_activediffQP] = setdiff(idx_activeQP,idx_active);
            H  = [ Gna * Fi - Sna; -lambdaHK(idx_activediffQP,1:nx); bndA ];
            K  = [ Wna - Gna * Gi;  lambdaHK(idx_activediffQP,nx+1); ...
                   bndb ];
            cr.P  = polytope(H,K);
            if ~isfulldim(cr.P),
                cr.type = 0;
            end
        end
    else
        BigB0 = Q \ Wa;
        Afbf  = U \ [P BigB0(1:nu)];
        Fi    = -Afbf(:,1:nx);
        Gi    =  Afbf(:,nx+1);
        H     = [bndA; Gna * Fi - Sna];
        K     = [bndb; Wna - Gna * Gi];
        cr.P  = polytope(H,K);

        % %++++++++++++++++++++>>>>>>>>>>> MODIFIED/CORRECTED/ADDED BY MATO BAOTIC 
        if ~isfulldim(cr.P),
            cr.type = 0;
            %disp('LO AND BEHOLD UNBELIEVABLE CAN HAPPEN!!!');
        end
        % %++++++++++++++++++++>>>>>>>>>>> MODIFIED/CORRECTED/ADDED BY MATO BAOTIC 
    end
    %
    cr.Fi   = Fi; 
    cr.Gi   = Gi;
    
    return
%
%============= END getCR ===========================    
    


%--------------------Raphael-------------------------
% change the function to:
%  1. compute all remaining extreme points of the bounding box
%  2. change the procedure to handle arbitrary dimensions
%========================================
function q = sub_whichquadrant(P,center)
%========================================
%    
Options.noPolyOutput = 1;
[R,l,u] = bounding_box(P,Options);

nx = dimension(P);
vert = ones(nx,2^nx);
binOne=dec2bin(1);
for i=1:2^nx
    decVec=dec2bin(i-1,nx);
    for j=1:nx
        if(decVec(j)==binOne)
            vert(j,i)=l(j);
        else
            vert(j,i)=u(j);
        end
    end
end

%-------------------My-------------------
q1 = zeros(1,2^nx);
for i=1:2^nx
    q2 = sub_whichquadrant_point(vert(:,i),center);
    q1(q2) = 1;
end
q = find(q1==1);

return
%
%========= END sub_whichquadrant ===============


%===================================================
%
function q = sub_whichquadrant_point(xBeyond,center)
%    
%===================================================
q  = [];
q1 = [];

nx = length(xBeyond);
%center = zeros(nx,1);       % center per default
twoQuadrant = find(xBeyond==0);

if isempty(twoQuadrant),
    %One = eye(nx);
    vi = (xBeyond - center) > 0;
    q1 = 1;
    for i=0:nx-1
        q1 = q1 + (vi(i+1))*2^i;    
    end   
else    
    for j=1:length(twoQuadrant)
        xBeyond(twoQuadrant(j))=1;
        %One = eye(nx);
        vi = (xBeyond - center) > 0;
        q = 1;
        for i=0:nx-1
            q = q + (vi(i+1))*2^i;    
        end 
        
        q1 = [q1 q];
        xBeyond(twoQuadrant(j))=-1;
        %One = eye(nx);
        vi = (xBeyond - center) > 0;
        q = 1;
        for i=0:nx-1
            q = q + (vi(i+1))*2^i;    
        end 
        xBeyond(twoQuadrant(j))=0;
        q1 = [q1 q]; 
    end      
end
q = q1;

return
%
% ============ END sub_whichquadrant_point =========


