function [U,feasible,region,cost,details]=mpt_getInput(ctrl,x0,Options)
%MPT_GETINPUT For a given state, extracts the (optimal) output from a controller structure
%
% [U,feasible,region,cost,details] = mpt_getInput(ctrlStruct,x0,Options)
%
% ------------------------------------------------------------------------------
% DESCRIPTION
% ------------------------------------------------------------------------------
% For the given state x0, this function extracts the optimal output from a
% controller given by means of the controller structure ctrlStruct. If the
% controller partition is overlapping in the X space, the input U will be picked
% up such that an associated cost is minimized (i.e. if for x0 there are 2 or
% more associated control laws, only the one which minimizes a given criterion
% is returned. The criterion is either value of the objective function for
% optimal solution or minimum time for the time-optimal solution).
%
% ------------------------------------------------------------------------------
% INPUT
% ------------------------------------------------------------------------------
% ctrl              - MPT controller
% x0                - initial state
% Options.openloop  - 0 by default. If set to 1, the full optimizer as obtained
%                     as a solution to the finite-time optimal control problem
%                     is returned, i.e. U = [u_0 u_1 ... u_N] where N is the
%                     prediction horizon
% Options.recover   - If set to 1 and there is no region associated to the
%                     current state x0, but the state itself lies in the
%                     feasible set of the controller, control law of the nearest
%                     neighbour is used. Default is 0.
% Options.abs_tol   - absolute tolerance
% Options.verbose   - Level of verbosity
% Options.useXU     - if 1, use a control input based on an XUset istead of the 
%                     usual (optimization) based control action. (default 0)
% Options.nlsolver  - which solver to use for nonlinear optimization.
%                       'global' - can be VERY slow, but gives results which are
%                                  close to global optimum (default)
%                       'local'  - faster than global, but does not
%                                  guarantee global optimality
% Options.nliter    - If global nonlinear solver is chosen, this parameter
%                     specifies number of iterations (default is 10)
% Options.MaxSQPIter - Maximum number of iterations for fmincon. If your
%                      nonlinear optimization seems to be stuck, try to set this
%                      number to something like 1e3.
% Options.lowersolver - Lower bound solver for YALMIP (e.g. 'cdd', 'glpk', 'clp')
% Options.uppersolver - Upper bound solver for YALMIP (e.g. 'fmincon', 'none')
%                       Set Options.uppersolver='none' if you don't have fmincon
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ------------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ------------------------------------------------------------------------------
% U         - control input computed as U=F*x0 + G (or via other strategies),
% region    - index of a region which contains the optimal control input
%             associated to the given state x0 
% cost      - value of the associated cost function
%             NOTE: The cost is not necessarily the cost of the objective
%                   function obtained for state x0, it can be also distance to
%                   the invariant set(s) in case of time-optimal controller!
% feasible  - returns 1 if  the there is at least one control law associated to
%             a given state x0, 0 otherwise
% details
%  .inwhich - vector of indicies of regions which contain state x0
%  .fullopt - full optimizer associated to the state x0
%  .runtime - runtime of the on-line MPC
%  .nops    - number of numerical operations (multiplications, summations,
%             comparisons) needed to identify and compute a control action
%             associated to a given x0 (only for explicit controllers)
%
% see also MPT_COMPUTETRAJECTORY, MPT_PLOTTIMETRAJECTORY
%

% Copyright is with the following author(s):
%
% (c) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (c) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (c) 2004 Arne Linder, Faculty of Electrical, Information and Media Engineering, 
%          Wuppertal University, alinder@uni-wuppertal.de

% ------------------------------------------------------------------------------
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
% ------------------------------------------------------------------------------

error(nargchk(2,3,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if ~isa(ctrl, 'mptctrl') & ~isstruct(ctrl)
    error('MPT_GETINPUT: first argument must be an MPTCTRL object!');
end

if ~isa(x0, 'double')
    error('MPT_GETINPUT: second argument must be vector!');
end

if nargin==3 & ~(isstruct(Options) | isempty(Options))
    error('MPT_GETINPUT: third argument must be an Options structure!');
end

if nargin<3,
    Options=[];
end

Options = mpt_defaultOptions(Options, ...
    'openloop', 0, ...
    'abs_tol', mptOptions.abs_tol, ...
    'convertconvex', 0, ...    
    'verbose', mptOptions.verbose, ...
    'recover', 0, ...
    'useXU', 0, ...
    'nlsolver', 'global', ...
    'lowersolver', '', ...
    'uppersolver', '', ...
    'nliter', 10, ...
    'MaxSQPIter', 'auto', ...
    'lpsolver', mptOptions.lpsolver, ...
    'qpsolver', mptOptions.qpsolver, ...
    'milpsolver', mptOptions.milpsolver, ...
    'miqpsolver', mptOptions.miqpsolver);

% according to Johan Loefberg, 'convertconvex' should be set to 0 for
% general nonlinear optimization. However then fmincon just stalls on the
% nonlin2.m demo (matlab R14 works fine, though)

if Options.useXU
    if isfield(ctrl.details,'XU')
        if ~isfield(ctrl.details.XU, 'XUset') | ~isfield(ctrl.details.XU, 'Idx_orig'),
            error('MPT_GETINPUT: provide XUset and/or Idx_orig in ctrl.details.XU')
        end
    else
        error('MPT_GETINPUT: provide XUset and/or Idx_orig in ctrl.details.XU')
    end
end


x0=x0(:);
nops.multiplications = 0;
nops.summations = 0;
nops.comparisons = 0;
details = struct('inwhich', [], 'fullopt', [], 'runtime', [], 'nops', nops);

if isa(ctrl, 'mptctrl') & ~isexplicit(ctrl)
    % solve an QP/LP/MIQP/MILP for on-line controllers

    % NOTE! NOTE! NOTE! from now on, we convert the "ctrl" object into a
    % structure to get faster access to internal fields
    ctrl = struct(ctrl);
    
    sysStruct = ctrl.sysStruct;
    probStruct = ctrl.probStruct;
    
    if isfield(ctrl.details, 'yalmipData'),
        % MPC problem for non-linear systems; use solvesdp()
        F = ctrl.details.yalmipData.F;
        obj = ctrl.details.yalmipData.obj;
        vars = ctrl.details.yalmipData.vars;
        if isfield(vars, 'uprev'),
            % assign value to previous control input
            if ~isfield(Options, 'Uprev'),
                fprintf('WARNING: deltaU constraints can only be satisfied if Options.Uprev is given.\n');
            end
            uprev = mpt_defaultField(Options, 'Uprev', zeros(ctrl.details.dims.nu, 1));
            uprev = uprev(:);
            if length(uprev) ~= ctrl.details.yalmipData.uprev_length,
                error(sprintf('Wrong size of "Options.uprev". Expecting %d elements, got %d.', ...
                    ctrl.details.yalmipData.uprev_length, length(uprev)));
            end
            F = F + set(vars.uprev == uprev);
        end
        if isfield(vars, 'ref')
            % assign value to reference
            if ~isfield(Options, 'reference'),
                fprintf('WARNING: Options.reference not specified, assuming zero.\n');
            end
            ref = mpt_defaultField(Options, 'reference');
            if isempty(ref),
                if isfield(ctrl.probStruct, 'Qy'),
                    % reference has the dimension of output
                    ref = zeros(ctrl.details.dims.ny, 1);
                else
                    % reference has the dimension of state
                    ref = zeros(ctrl.details.dims.nx, 1);
                end
            end
            ref = ref(:);
            if length(ref) ~= ctrl.details.yalmipData.reference_length,
                error(sprintf('Wrong size of "Options.reference". Expecting %d elements, got %d.', ...
                    ctrl.details.yalmipData.reference_length, length(ref)));
            end
            F = F + set(vars.ref == ref);
        end

        % assign value to x0
        F = F + set(vars.x{1} == x0);
        
        % now solve the non-linear problem
        yalmipOptions = mptOptions.sdpsettings;
        if isequal(Options.nlsolver, 'global'),
            % switch to a global solver. WARNING! global solvers are extremely
            % slow!
            yalmipOptions.solver = 'bmibnb';
        end
        if isequal(Options.MaxSQPIter, 'auto'),
            % automatically set this option based on the precise number of
            % variables and constraints. similar setting is used in optimset()
            % under Matlab R14.
            nvars = length(sdpvar(F));
            nineqs = length(unique([depends(F) depends(obj)]));
            Options.MaxSQPIter = 5*max(nvars, nineqs);
        end
        yalmipOptions.convertconvexquad = Options.convertconvex;
        yalmipOptions.verbose = Options.verbose;
        yalmipOptions.fmincon.MaxSQPIter = Options.MaxSQPIter;
        yalmipOptions.bmibnb.lowersolver = Options.lowersolver;
        yalmipOptions.bmibnb.uppersolver = Options.uppersolver;
        yalmipOptions.bmibnb.maxiter = Options.nliter;
        if isfield(ctrl.details.yalmipData, 'model_expanded'),
            % the model was already expanded in mptctrl()
            yalmipOptions.expand = 0;
        end
        
        starttime = cputime;
        diagnost = solvesdp(F, obj, yalmipOptions);
        runtime = cputime - starttime;
        
        feasible = (diagnost.problem==0);
        fullopt = [];
        U = double(cat(1, vars.u{:}));
        cost = double(obj);
        
    elseif isfield(ctrl.details, 'yalmipMatrices'),
        % we do have an optimization problem constructed by mpt_yalmipcftoc(),
        % solve the corresponding LP/QP/MILP/MIQP problem
        M = ctrl.details.yalmipMatrices;
        mi_problem = ~isempty(M.binary_var_index);
        qp_problem = M.qp;
        
        %==================================================================
        % set default values of Options.Uprev and/or Options.reference
        if M.uprev_length > 0 & (length(x0)~=size(M.E, 2)),
            % include previous input u(k-1) into the state vector if needed
            if ~isfield(Options, 'Uprev'),
                fprintf('WARNING: deltaU constraints can only be satisfied if Options.Uprev is given.\n');
            end
            uprev = mpt_defaultField(Options, 'Uprev', zeros(ctrl.details.dims.nu, 1));
            uprev = uprev(:);
            if length(uprev) ~= M.uprev_length,
                error(sprintf('Wrong size of "Options.uprev". Expecting %d elements, got %d.', ...
                    M.uprev_length, length(uprev)));
            end
            x0 = [x0; uprev];
        end
        if M.reference_length > 0 & (length(x0)~=size(M.E, 2)),
            % include reference into the state vector if needed
            if ~isfield(Options, 'reference'),
                fprintf('WARNING: Options.reference not specified, assuming zero.\n');
            end
            ref = mpt_defaultField(Options, 'reference');
            if isempty(ref),
                if isfield(ctrl.probStruct, 'Qy'),
                    % reference has the dimension of output
                    ref = zeros(ctrl.details.dims.ny, 1);
                else
                    % reference has the dimension of state
                    ref = zeros(ctrl.details.dims.nx, 1);
                end
            end
            ref = ref(:);
            if length(ref) ~= M.reference_length,
                error(sprintf('Wrong size of "Options.reference". Expecting %d elements, got %d.', ...
                    M.reference_length, length(ref)));
            end
            x0 = [x0; ref];
        end

        %==================================================================
        % prepare constraints, plug in the parametric variables
        if length(x0(:)) ~= size(M.E, 2),
            error(sprintf('Wrong dimension of x0 (expecting %d elements)', size(M.E, 2)));
        end
        nvars = size(M.G, 2);
        nparams = length(M.param_var);
        A = M.G;
        B = M.W + M.E*x0;
        if isempty(M.Aeq),
            % happens if we used mpt_remove_equalities in mptctrl.m
            Aeq = M.Aeq;
            Beq = M.beq;
        else
            Aeq = M.Aeq;
            Beq = M.beq - M.Beq*x0;
        end
        lb = M.lb(1:end-nparams);  % bounds on parametric variables are always last!
        ub = M.ub(1:end-nparams);


        %==================================================================
        % denote respective variables as binary, enforce bounds as constraints
        % if necessary (for mpt_solveLP and mpt_solveQP)
        if mi_problem,
            % denote selected variables as binary
            vartype = repmat('C', nvars, 1);
            vartype(M.binary_var_index) = 'B';
            
        elseif ~isempty(lb)
            % mpt_solveLP() and mpt_solveQP() currently do not support
            % lower/upper bound on variables, therefore we have to convert them
            % into dummy constraints
            A = [A; eye(nvars); -eye(nvars)];
            B = [B; ub; -lb];
            lb = []; ub = [];
        end
        
        %==================================================================
        % solve the problem
        starttime = cputime;
        PARAM = [];
        xinit = mpt_defaultField(Options, 'usex0', []);
        if qp_problem,
            % QP/MIQP problem
            if mi_problem,
                [xopt, cost, how, exitflag] = mpt_solveMIQP(...
                    M.H, M.F'*x0+M.Cf', A, B, Aeq, Beq, lb, ub, vartype, ...
                    PARAM, Options, Options.miqpsolver);
                
            else
                [xopt, lambda, how, exitflag, cost] = mpt_solveQP(...
                    M.H, M.F'*x0+M.Cf', A, B, Aeq, Beq, xinit, ...
                    Options.qpsolver);
                
            end
            cost = cost + x0'*M.Y*x0 + M.Cx*x0 + M.Cc;
            
        else
            % LP/MILP problem
            if mi_problem,
                [xopt, cost, how, exitflag] = mpt_solveMILP(...
                    M.H, A, B, Aeq, Beq, lb, ub, vartype, PARAM, Options, ...
                    Options.milpsolver);
                
            else
                [xopt, cost, lambda, exitflag, how] = mpt_solveLP(...
                    M.H, A, B, Aeq, Beq, xinit, Options.lpsolver);
                
            end
            
        end
        runtime = cputime - starttime;
        feasible = strcmp(how, 'ok');


        if ~isempty(M.getback),
            % if mpt_project_on_equality() was used (done in mptctrl.m),
            % a different expressions must be used to obtain the extract
            % the correct optimizer
            xopt = M.getback.S1*xopt + M.getback.S2*x0 + M.getback.S3;
        else
            % otherwise we take the whole optimizer
            xopt = xopt(:);
        end
        
        U = xopt(M.requested_variables);
        fullopt = xopt;
        
    elseif iscell(sysStruct.A),
        % PWA system
        
        if ctrl.details.dims.nx~=length(x0),
            error(sprintf('MPT_GETINPUT: state x0 should be a column vector with %d elements!',ctrl.details.dims.nx));
        end
        
        if Options.verbose<=1,
            Options.verbose=0;
        end
        if ~isfield(sysStruct.data,'MLD')
            error('MPT_GETINPUT: no MLD model available in system structure! Did you use ''mpt_sys(hysdelmodel)''?');
        end
        
        ny = sysStruct.data.MLD.ny;
        nx = sysStruct.data.MLD.nx;
        
        % define weights
        if isfield(probStruct, 'Qy')
            weights.Qy = probStruct.Qy;
        else
            weights.Qx = probStruct.Q;
        end
        if isfield(probStruct, 'Qz'),
            % penalty on "z" variables
            weights.Qz = probStruct.Qz;
        end
        if isfield(probStruct, 'Qd'),
            % penalty on "d" variables
            weights.Qd = probStruct.Qd;
        end
        weights.Qu = probStruct.R;

        % defined options
        if isfield(ctrl.details, 'Matrices'),
            % use pre-calculated problem matrices stored in the controller
            if ~isempty(ctrl.details.Matrices),
                if isfield(ctrl.details.Matrices, 'F1') & isfield(ctrl.details.Matrices, 'IntIndex'),
                    % double-check that the matrices contain appropriate fields
                    Options.problemmatrices = ctrl.details.Matrices;
                end
            end
        end
        Options.norm = probStruct.norm;
        Options.eps2 = 1e3; % tolerance on fulfillment of terminal state constraint
        Options.umin = sysStruct.umin;
        Options.umax = sysStruct.umax;
        if any(~isinf(sysStruct.dumax)) | any(~isinf(sysStruct.dumin))
            Options.dumax = sysStruct.dumax;
            Options.dumin = sysStruct.dumin;
        end
        if isfield(sysStruct, 'xmax') & isfield(sysStruct, 'xmin'),
            Options.xmax = sysStruct.xmax;
            Options.xmin = sysStruct.xmin;
        end
        if isfield(sysStruct, 'ymax') & isfield(sysStruct, 'ymin'),
            Options.ymax = sysStruct.ymax;
            Options.ymin = sysStruct.ymin;
        end
        
        Options.TerminalConstraint = 0;
        if probStruct.Tconstraint==2 & isfield(probStruct, 'Tset'),
            if isfulldim(probStruct.Tset),
                % tell mpc_mip to include terminal set constraint
                Options.Tset = probStruct.Tset;
            end
        end
        
        % define references
        if isfield(Options, 'reference') & probStruct.tracking>0,
            if isfield(probStruct, 'Qy')
                % output reference
                ref.y = Options.reference;
            else
                % state reference
                ref.x = Options.reference;
            end
        else
            if isfield(probStruct, 'Qy'),
                if isfield(probStruct, 'yref'),
                    ref.y = probStruct.yref;
                else
                    ref.y = zeros(ny,1);
                end
            else
                if isfield(probStruct, 'xref'),
                    ref.x = probStruct.xref;
                    if isfield(probStruct, 'uref')
                        ref.u = probStruct.uref;
                    end
                else
                    ref.x = zeros(nx,1);
                end
            end
        end
        
        if isfield(probStruct, 'xN'),
            clear ref
            ref.x = probStruct.xN;
            % tell mpc_mip to include terminal set constraint
            Options.TerminalConstraint = 1;
            % zero tolerance on satisfaction of terminal set constraint
            Options.eps2 = Options.abs_tol;
        end

        % in case of time-varying penalties, mpc_mip expects them to be in a
        % cell array
        %   weights{1}.Qx (.Qy, .Qu, .Qd, .Qz)
        %   ...
        %   weights{N}.Qx (.Qy, .Qu, .Qd, .Qz)
        % we make the conversion in a subfunction...
        weights = sub_fixweights(weights);
        
        % if time-varying penalties is used, mpc_mip wants to have one MLD
        % model per each sampling time, so we make it happy...
        S = {};
        for ii = 1:probStruct.N,
            S{ii} = sysStruct.data.MLD;
        end

        % the prediction horizon in mpc_mip is determined as sum(horizon)...
        horizon = repmat(1, 1, probStruct.N);
        
        [ut,dt,zt,Eflag] = mpc_mip(S, x0, ref, weights, horizon, horizon, Options);
        
        feasible = strcmp(Eflag.solverflag, 'ok');
        cost = Eflag.fopt;
        U = Eflag.u(:);
        details.fullopt = Eflag.full_xopt;
        details.runtime = Eflag.runtime;
    else
        % LTI system
        isycost = isfield(probStruct, 'Qy');
        if isycost,
            reflength = ctrl.details.dims.ny;
        else
            reflength = ctrl.details.dims.nx;
        end
        if ctrl.details.dims.nx~=length(x0),
            if ctrl.probStruct.tracking==1,
                if length(x0) ~= (ctrl.details.dims.nx + ctrl.details.dims.nu + reflength),
                    error('For tracking==1, the state vector x0 must be in form [x; u; ref] !');
                end
            elseif ctrl.probStruct.tracking==2,
                if length(x0) ~= (ctrl.details.dims.nx + reflength)
                    error('For tracking==2, the state vector x0 must be in form [x; ref] !');
                end
            else
                error(sprintf('MPT_GETINPUT: state x0 should be a column vector with %d elements!',ctrl.details.dims.nx));
            end
        end
        
        Matrices = ctrl.details.Matrices;
        [U, feasible, cost] = mpt_solveMPC(x0, sysStruct, probStruct, Matrices);
    end
    
    nu = ctrl.details.dims.nu;
    if ~feasible,
        U = repmat(NaN, nu*probStruct.N, 1);
    end
    % if control horizon was used, make sure that we have the full open-loop
    % optimizer at disposal (if mpt_yalmipcftoc() was used to do move blocking,
    % the optimizer would just consist of free control moves, therefore we add
    % the remaining manipulated variables at the end)
    if isfield(probStruct, 'Nc') & (length(U(:)) < probStruct.N*nu),
        ulast = U(end-nu+1:end);
        U = [U; repmat(ulast, probStruct.N-probStruct.Nc, 1)];
    end
    if Options.openloop==0,
        U = U(1:nu);
    end

    region = [];
    return
    
    
else
    % NOTE! NOTE! NOTE! from now on, we convert the "ctrl" object into a
    % structure to get faster access to internal fields
    ctrlStruct = struct(ctrl);
    
end


% evaluate explicit control law for state x0
U=[];
cost=-Inf;
feasible=0;
region=0;

if Options.useXU
    [U, feasible, region, XUreg] = mpt_getInputXU(ctrlStruct, x0, Options);
    cost = 0;
    return
end

PA = ctrlStruct.Pn;

sysStruct = ctrlStruct.sysStruct;

if dimension(PA)~=length(x0),
    if ctrlStruct.probStruct.tracking==1,
        disp('For tracking==1, the state vector x0 must be in form [x; u; ref] !');
    elseif ctrlStruct.probStruct.tracking==2,
        disp('For tracking==2, the state vector x0 must be in form [x; ref] !');
    end
    error(sprintf('MPT_GETINPUT: state x0 should be a column vector with %d elements!',dimension(PA)));
end

%nu = ctrl.details.dims.nu;
if iscell(sysStruct.B),
    ispwa = 1;
    nu = size(sysStruct.B{1},2);
else
    ispwa = 0;
    nu = size(sysStruct.B,2);
end


U = []; inwhich = [];
cost=-Inf;
feasible=0;
region=0;

if Options.recover,
    % if a state lies in the feasible state-space but no region is associated to
    % it, we take the control law of the nearest neighbour.
    
    % check if the given point lies inside of the feasible state-space. use
    % fastbreak=1 to tell isinside() to abort quickly once it finds a region
    Options.fastbreak = 1;
    isin = isinside(ctrlStruct.Pfinal, x0, Options);
    if ~isin,
        % the state does not lie in the feasible state-space. we handle this
        % case later
    else
        % x0 does lie in the feasible state space, but maybe in a hole?
        % check that. use fastbreak=0 such that we can choose region with lowest
        % cost.
        locOpt.fastbreak = 0;
        [isin, inwhich, closest] = isinside(ctrlStruct.Pn,x0,Options);
        if ~isin,
            % indeed, the current state x0 lies in a hole, use feedback law of
            % neighbouring region
            if Options.verbose>0,
                fprintf('MPT_GETINPUT: State x = [%s] lies in a hole, taking the nearest neighbour\n',num2str(x0'));
            end
            isin = 1;
            inwhich = closest;
        end
    end
else
    % this is a safe setting, we _could_ eventually abort the search for active
    % region once we find at least one such region, but:
    %  * in PWA / CFTOC / linear norms case it is very likely that the state ends up
    %    on a boundary of multiple regions, therefore we have to compare costs and
    %    pick up region where the cost is minimal
    %  * in minimum-time case we could break after the first found region, because
    %    all regions are already sorted in a "good" way, but what if somebody
    %    modifies the controller with modify() and includes his own regions?
    %
    % bottom line, we don't loose much if we always search through _all_ regions
    locOpt.fastbreak = 0;
    
    [isin, inwhich] = isinside(ctrlStruct.Pn,x0,Options);
end
details.inwhich = inwhich;

if ~isin
    % no associated control law
    feasible=0;
    if Options.verbose>0,
        if isfield(ctrlStruct.sysStruct, 'xmax'),
            if (any(x0 > ctrlStruct.sysStruct.xmax) | any(x0 < ctrlStruct.sysStruct.xmin)),
                fprintf('\nState x = %s violates state constraints!\n', mat2str(x0));
            end
        end
        fprintf('MPT_GETINPUT: NO REGION FOUND FOR STATE x = %s\n', mat2str(x0));
    end
    return
elseif isfield(ctrlStruct.details,'searchTree'),
    % use search tree for fast region identification
    % contributed by Arne Linder
    searchTree = ctrlStruct.details.searchTree;
    [lenST colST]=size(searchTree);
    node=1; niter = 0;
    while node>0  % node<0 means node is now number of control law
        niter = niter + 1;
        H=searchTree(node,1:colST-3);
        K=searchTree(node,colST-2);
        if H*x0-K<0
            node=searchTree(node,colST);    % x0 on plus-side of the hyperplane
        else
            node=searchTree(node,colST-1);  % x0 on minus-side of the hyperplane
        end
    end

    node = -round(node);
    
    U = ctrlStruct.Fi{node}*x0+ctrlStruct.Gi{node};
    cost = x0'*ctrlStruct.Ai{node}*x0 + ctrlStruct.Bi{node}(:)'*x0 + ctrlStruct.Ci{node};
    region = node;
    feasible = 1;
    if ~Options.openloop
        % return just U at time 0
        U = U(1:nu);
    end

    %====================================================================
    % compute number of operations needed to find the optimal control law
    
    % H*x0, where H is a row vector needs "nx" multiplications
    details.nops.multiplications = niter*length(x0);
        
    % H*x0   requires "nx" summations
    % H*x0-K requires 1 additional summation (adding -K)
    details.nops.summations = niter*(length(x0) + 1);
        
    % comparing (H*x0 - K)<0 needs 1 comparison, since H is a row vector
    details.nops.comparisons = niter;
        
    % F*x0 + G
    details.nops.multiplications = details.nops.multiplications + nu*length(x0);
    details.nops.summations = details.nops.summations + nu + nu*length(x0);
    
    return

elseif length(inwhich)==1,
    % x0 lies only in one region, return the control law and cost evaluated at x0
    U = ctrlStruct.Fi{inwhich}*x0 + ctrlStruct.Gi{inwhich};
    cost = x0'*ctrlStruct.Ai{inwhich}*x0 + x0'*ctrlStruct.Bi{inwhich}(:) + ctrlStruct.Ci{inwhich};
    feasible = 1;
    region=inwhich;
else
    % x0 belongs to several regions, go through all of them and isolate the control law with least associated cost
    mincost = Inf;      % minimal cost
    mincostregion = 0;  % index of a region in which cost associated to the state x0 is minimal
    for ii=length(inwhich):-1:1
        region = inwhich(ii);
        cost = x0'*ctrlStruct.Ai{region}*x0 + ctrlStruct.Bi{region}(:)'*x0 + ctrlStruct.Ci{region};
        if cost <= mincost,
            mincost = cost;
            mincostregion = region;
        end
    end
    U = ctrlStruct.Fi{mincostregion}*x0 + ctrlStruct.Gi{mincostregion};
    cost = mincost;
    feasible = 1;
    region = mincostregion;
end
if ctrlStruct.probStruct.feedback,
    % in case of pre-stabilization, actual control move is:
    % U = (Fi + K)*x0 + Gi
    U_nofeedback=U;
    x0_i=x0;
    for i=1:(length(U)/nu)
        U(((i-1)*nu+1):i*nu) = U(((i-1)*nu+1):i*nu) + ctrlStruct.probStruct.FBgain*x0_i;
        x0_i = ctrlStruct.sysStruct.A*x0_i + ctrlStruct.sysStruct.B*U(((i-1)*nu+1):i*nu);
    end
end

if isfield(ctrlStruct.probStruct,'inputblocking'),
    if(~isempty(ctrlStruct.probStruct.inputblocking))
        inputblocking=ctrlStruct.probStruct.inputblocking;
        u=[];
        x0_i=x0;
        for i=1:length(inputblocking)
            for j=1:inputblocking(i)
                if ~ctrlStruct.probStruct.feedback
                    u=[u; U(((i-1)*nu+1):i*nu)];    
                elseif ctrlStruct.probStruct.feedback
                    u = [u; U_nofeedback(((i-1)*nu+1):i*nu) + ctrlStruct.probStruct.FBgain*x0_i];
                    x0_i = ctrlStruct.sysStruct.A*x0_i + ctrlStruct.sysStruct.B*u(end+1-nu:end);    
                end
            end
        end
        U=u;
    end
end
if isfield(ctrlStruct.probStruct,'deltablocking'),
    if(~isempty(ctrlStruct.probStruct.deltablocking))
        deltablocking=ctrlStruct.probStruct.deltablocking;
        u=[];
        x0_i=x0;
        if ctrlStruct.probStruct.feedback
            U=U_nofeedback;
        end
        u=U(1:nu);
        for i=2:length(deltablocking)
            u(((deltablocking(i)-1)*nu+1):deltablocking(i)*nu)=U(((i-1)*nu+1):i*nu);
            for j=(deltablocking(i-1)+1):(deltablocking(i)-1)
                u(((j-1)*nu+1):j*nu)=u((((j-1)*nu+1):j*nu)-nu)+(U(((i-1)*nu+1):i*nu)-U((((i-1)*nu+1):i*nu)-nu))/(deltablocking(i)-deltablocking(i-1));
            end
        end
        if ctrlStruct.probStruct.feedback
            for i=1:(length(u)/nu)
                u(((i-1)*nu+1):i*nu) = u(((i-1)*nu+1):i*nu) + ctrlStruct.probStruct.FBgain*x0_i;
                x0_i = ctrlStruct.sysStruct.A*x0_i + ctrlStruct.sysStruct.B*u(((i-1)*nu+1):i*nu);
            end
        end
        
        U=u;
    end
end

% if control horizon was used, make sure that we have the full open-loop
% optimizer at disposal (if mpt_yalmipcftoc() was used to do move blocking,
% the optimizer would just consist of free control moves, therefore we add
% the remaining manipulated variables at the end)
if isfield(ctrl.probStruct, 'Nc') & (length(U(:)) < ctrl.probStruct.N*nu),
    ulast = U(end-nu+1:end);
    U = [U; repmat(ulast, ctrl.probStruct.N-ctrl.probStruct.Nc, 1)];
end

if ~Options.openloop
    % return just U at time 0
    U = U(1:nu);
end

if nargout<5,
    clear details
end
return


%---------------------------------------------------
function W = sub_fixweights(weights),
% in case of time-varying penalties, mpc_mip expects them to be in a
% cell array
%   weights{1}.Qx (.Qy, .Qu, .Qd, .Qz)
%   ...
%   weights{N}.Qx (.Qy, .Qu, .Qd, .Qz)
% we make the conversion here

Qx_cell = 0;
Qy_cell = 0;
Qu_cell = 0;
Qz_cell = 0;
Qd_cell = 0;
haveQy  = 0;
haveQz  = 0;
haveQd  = 0;
Qx_length = 0;
Qy_length = 0;
Qu_length = 0;
Qz_length = 0;
Qd_length = 0;

if isfield(weights, 'Qx'),
    if iscell(weights.Qx),
        Qx_cell = 1;
        Qx_length = length(weights.Qx);
    end
end
if isfield(weights, 'Qu'),
    if iscell(weights.Qu),
        Qu_cell = 1;
        Qu_length = length(weights.Qu);
    end
end
if isfield(weights, 'Qy'),
    haveQy = 1;
    if iscell(weights.Qy),
        Qy_cell = 1;
        Qy_length = length(weights.Qy);
    end
end
if isfield(weights, 'Qz'),
    haveQz = 1;
    if iscell(weights.Qz),
        Qz_cell = 1;
        Qz_length = length(weights.Qz);
    end
end
if isfield(weights, 'Qd'),
    haveQd = 1;
    if iscell(weights.Qd),
        Qd_cell = 1;
        Qd_length = length(weights.Qd);
    end
end

if Qx_cell | Qu_cell | Qy_cell | Qz_cell | Qd_cell,
    N = max([Qx_length Qu_length Qy_length Qz_length Qd_length]);
    W = {};
    for ii = 1:N,
        if Qx_cell,
            W{ii}.Qx = weights.Qx{ii};
        else
            W{ii}.Qx = weights.Qx;
        end
        if Qu_cell,
            W{ii}.Qu = weights.Qu{ii};
        else
            W{ii}.Qu = weights.Qu;
        end
        if haveQy,
            if Qy_cell,
                W{ii}.Qy = weights.Qy{ii};
            else
                W{ii}.Qy = weights.Qy;
            end
        end
        if haveQz,
            if Qz_cell,
                W{ii}.Qz = weights.Qz{ii};
            else
                W{ii}.Qz = weights.Qz;
            end
        end
        if haveQd,
            if Qd_cell,
                W{ii}.Qd = weights.Qd{ii};
            else
                W{ii}.Qd = weights.Qd;
            end
        end
    end
else
    W = weights;
end
