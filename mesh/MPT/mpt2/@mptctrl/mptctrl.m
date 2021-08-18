function ctrl = mptctrl(varargin)
%MPTCTRL Constructor for the MPT controller object
%
% ctrl = mptctrl(ctrlStruct)
% ctrl = mptctrl(sysStruct, probStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Creates an MPT controller object. The controller can represent both an
% explicit control law (i.e. a Piecewise-Affine feedback law defined over a
% polyhedral partition of a state-space), or an implicit controller (i.e. an MPC
% controller where the optimization problem is being solved at every time step
% for current measurements of the states).
%
% USAGE:
%
%   ctrl = mptctrl(ctrlStruct)
%
%     If called with one input parameter, the input is assumed to be a
%     controller structure, as defined in the MPT manual. The output will then
%     be an explicit controller.
%
%   ctrl = mptctrl(sysStruct, probStruct)
%
%     If called with two input parameters, creates an on-line controller.
%
% ---------------------------------------------------------------------------
% REPRESENTATION OF CONTROLLERS
% ---------------------------------------------------------------------------
%
% Each MPTCTRL object is internally represented as a structure with following
% fields:
%
%   .type        - type of the controller {'explicit' | 'online'}
%   .sysStruct   - system structure
%   .probStruct  - problem structure
%   .Pn          - polyhedral partition over which the explicit controller is
%                  defined (POLYTOPE object)
%   .Pfinal      - feasible set of the controller, i.e. set of states for which
%                  there exists a control law (POLYTOPE object) 
%
%   .Fi          - if the given state "x" is in region Pn(i), the associated
%   .Gi          /    control law is given by U = Fi{i}*x + Gi{i}
%
%   .Ai          \
%   .Bi          - value of the objective function at point 'x' is given by
%   .Ci          /    J = x'*Ai{i}*x + Bi{i}*x + Ci{i}
%
%   .dynamics    - a vector which associates dynamics of a PWA system to a given
%                     region of controller parition. 
%   .details     - additional details about the solution
%   .overlaps    - binary flag, takes a "true" value if regions of the
%                     controller partition overlap, "false" otherwise
%   .simplified  - binary flag, takes a "true" value if regions of the
%                     controller partition have been simplified using greedy or
%                     optimal merging, "false" otherwise 
%
%
% ---------------------------------------------------------------------------
% ACCESING INTERNAL DATA OF A CONTROLLER OBJECT
% ---------------------------------------------------------------------------
%
% Even though the MPT controller is represented as an object, it's fields can
% still be accessed as it were a structure. That means that you can do, for
% instance, the following:
%
%   ctrl.Pn
%   ctrl.details
%
% and so on. There are also couple of short-cut functions defined to quickly
% access important fields:
%
% nR = length(ctrl)    - returns number of regions of the controller partition
% time = runtime(ctrl) - returns runtime
% isexplicit(ctrl)     - returns 1 if the controller is an explicit one, 0 if it
%                          is an on-line controller
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% either: 
%   ctrlStruct  - valid controller structure (see manual for more details)
%
% or:
%   sysStruct   - valid system structure
%   probStruct  - valid problem structure
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrl          - an MPTCTRL object
%
%
% see also MPTCTRL/ANALYZE, MPTCTRL/ISEXPLICIT, MPTCTRL/LENGTH, MPTCTRL/PLOT
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

global mptOptions
if ~isstruct(mptOptions),
    mpt_error;
end

if nargin > 0,
    % if input is an MPTCTRL object, do nothing
    if isa(varargin{1}, 'mptctrl')
        ctrl = varargin{1};
        return
    end
end

% =========================================================================
% only these fields are allowed to be present in a ctrlStruct structure
allowedfields = {'sysStruct', 'probStruct', 'Pn', 'Pfinal', 'Fi', 'Gi', 'Ai', 'Bi', ...
        'Ci', 'dynamics', 'details', 'overlaps', 'simplified', 'type'};

% =========================================================================
% set default values
ctrl.sysStruct = [];
ctrl.probStruct = [];
ctrl.Pn = polytope;
ctrl.Pfinal = polytope;
ctrl.Fi = {};
ctrl.Gi = {};
ctrl.Ai = {};
ctrl.Bi = {};
ctrl.Ci = {};
ctrl.dynamics = [];
ctrl.details = [];
ctrl.overlaps = 0;
ctrl.simplified = 0;
ctrl.type = [];

if nargin==0,
    ctrl = class(ctrl, 'mptctrl');
    return
end

% =========================================================================
% first input has always to be a structure or a cell array (if we have
% multi-model dynamics)
if ~(isstruct(varargin{1}) | iscell(varargin{1}))
    error('MPTCTRL: first input must be a structure!');
end


if isfield(varargin{1}, 'Pn') & isfield(varargin{1}, 'Fi')
    
    % =========================================================================
    % this structure is a controller structure
    
    ctrlStruct = varargin{1};
    
    if mpt_isValidCS(ctrlStruct, struct('nowarnings',1)),
        % copy only allowed fields from ctrlStruct to ctrl, any remaining fields
        % will be stored to ctrl.details.additional
        fnames = fields(ctrlStruct);
        for ii=1:length(fnames)
            fval = getfield(ctrlStruct, fnames{ii});
            ctrl = setfield(ctrl, fnames{ii}, fval);
        end
        ctrl = copyvalidfields(ctrl, allowedfields);
        
    else
        error('MPTCTRL: first input is not a valid controller structure!');
    end

    if isfield(varargin{1}, 'type'),
        ctrl.type = varargin{1}.type;
    else
        ctrl.type = 'explicit';
    end
    
elseif nargin==1 & mpt_issysstruct(varargin{1}),
    
    % =========================================================================
    % input is a system structure

    sysStruct = mpt_verifySysStruct(varargin{1});
    
    if ~isstruct(sysStruct.A),
        sysStruct = mpt_lti2pwa(sysStruct);
    end
    if nnz([sysStruct.B{:}]) > 0 | nnz(cat(1, sysStruct.guardU{:})) > 0,
        error('MPTCTRL: only autonomous systems can be converted.');
    end
    if isfield(sysStruct, 'nonlinhandle'),
        error('MPTCTRL: nonlinear systems cannot be converted.');
    end
    [nx, nu, ny, ndyn, nbool, ubool, intInfo] = mpt_sysStructInfo(sysStruct);
    
    ctrl.sysStruct = sysStruct;
    probStruct = struct('Q', eye(nx), 'R', eye(nu), 'N', 1, 'norm', 1, 'subopt_lev', 0);
    ctrl.probStruct = mpt_verifyProbStruct(probStruct);
    ctrl.details.origSysStruct = sysStruct;
    ctrl.details.origProbStruct = probStruct;
    
    ctrl.Pn = polytope(intInfo.PdynX);
    ctrl.Pfinal = ctrl.Pn;
    ctrl.Fi = cell(1, ndyn); [ctrl.Fi{:}] = deal(zeros(nu, nx));
    ctrl.Gi = cell(1, ndyn); [ctrl.Gi{:}] = deal(zeros(nu, 1));
    ctrl.Ai = cell(1, ndyn); [ctrl.Ai{:}] = deal(zeros(nx, nx));
    ctrl.Bi = cell(1, ndyn); [ctrl.Bi{:}] = deal(zeros(nx, 1));
    ctrl.Ci = cell(1, ndyn); [ctrl.Ci{:}] = deal(0);
    ctrl.dynamics = zeros(1, ndyn);
    ctrl.details.runTime = 0;
    ctrl.type = 'explicit';
    
elseif nargin==2 | nargin==3
    
    % =========================================================================
    % two inputs, possibly sysStruct + probStruct case, verify that
    sysStruct = varargin{1};
    probStruct = varargin{2};
    userSysStruct = sysStruct;    
    if iscell(userSysStruct),
        sysStruct = userSysStruct{1};
    end
    if ~isstruct(sysStruct) | ~isstruct(probStruct)
        error('MPTCTRL: both input arguments must be structures!');
    end

    % =========================================================================
    % verify sysStruct and probStruct, catch any errors
    if ~isfield(sysStruct, 'verified') | ~isfield(probStruct, 'verified'),
        try
            [sysStruct, probStruct] = mpt_verifySysProb(sysStruct, probStruct);
        catch
            rethrow(lasterror);
        end
    end
    
    ctrl.sysStruct = sysStruct;
    ctrl.probStruct = probStruct;
    ispwa = iscell(sysStruct.A);
    
    % =========================================================================
    % check if on-line solution can be computed
    if mpt_isnoise(sysStruct.noise) & iscell(sysStruct.A)
        error('MPTCTRL: no on-line solution for PWA systems with additive noise!');
    elseif isfield(sysStruct, 'Aunc') & iscell(sysStruct.A)
        error('MPTCTRL: no on-line solution for PWA systems with parametric uncertainty!');
    elseif probStruct.subopt_lev==1,
        error('MPTCTRL: no on-line solution to minimum-time problems! Please set ''probStruct.subopt_lev=0''.');
    elseif probStruct.subopt_lev==2,
        error('MPTCTRL: no on-line solution to low-complexity problems! Please set ''probStruct.subopt_lev=0''.');
    elseif isinf(probStruct.N)
        error('MPTCTRL: no on-line solution to infinite-horizon problems! Please set ''probStruct.N=finite''.');
    end
    
    ctrl.details.origSysStruct = sysStruct;
    ctrl.details.origProbStruct = probStruct;
    
    % =========================================================================
    % some additional tests for conditions under which an on-line controller
    % cannot be computed
    if probStruct.feedback 
        if iscell(sysStruct.A),
            error('MPTCTRL: Feedback Pre-Stabilization not supported for PWA systems!');
        end
        if ~isfield(probStruct,'FBgain')
            % compute the pre-stabilization feedback for LTI system if it is not given
            [FB,S,E] = mpt_dlqr(sysStruct.A,sysStruct.B,probStruct.Q,probStruct.R);
            probStruct.FBgain = -FB;
        end
    end

    if nargin==3,
        Options = varargin{end};
    else
        Options = [];
    end
    Options = mpt_defaultOptions(Options, ...
        'verbose', mptOptions.verbose, ...
        'yalmip_data', [], ...
        'remove_equalities_always', 0, ...
        'remove_equalities_lti', 1, ...
        'force_mpt', 0 );

    if Options.force_mpt == 0 | ~isempty(Options.yalmip_data) | iscell(userSysStruct),
        % we can use mpt_yalmipcftoc(), do so
        
        if isempty(Options.yalmip_data),
            if Options.verbose > 1,
                fprintf('Using mpt_yalmipcftoc() to construct the optimization problem.\n');
            end
            % first construct constraints and optimization objective
            Options.yalmip_online = 1;
            Options.dont_solve = 1;
            Options.dp = 0;
            [dummy, F, obj, vars] = mpt_yalmipcftoc(userSysStruct, probStruct, Options);
            
            % mpt_yalmipcftoc() might have ran mpt_yalmipTracking() or
            % mpt_yalmipDU() which changes sysStruct/probStruct
            sysStruct = dummy.sysStruct;
            probStruct = dummy.probStruct;
            
        else
            % use data provided from mpt_owncost()
            F = Options.yalmip_data.constraints;
            obj = Options.yalmip_data.objective;
            vars = Options.yalmip_data.variables;
            dummy = [];
            
        end

        % now convert the optimization problem into MPTs native form
        % vars.x{1} corresponds to x0
        % vars.uprev corresponds to previous input (for deltaU constraints)
        uprev_length = 0; reference_length = 0;
        parametric_vars = vars.x{1};
        fprintf('\n');
        if isfield(vars, 'uprev'),
            fprintf('For deltaU constraints to hold, you will need to provide\n');
            fprintf('"Options.Uprev" when calling "u=ctrl(x0, Options)".\n');
            parametric_vars = [parametric_vars; vars.uprev];
            uprev_length = length(vars.uprev);
   
        elseif isfield(dummy, 'uprev_in_x0'),
            % u(k-1) is already part of the state vector, hence parametric_vars
            % don't need to be updated
            fprintf(['The controller will return deltaU=u(k)-u(k-1) instead of '...
                    'u(k) when evaluated!\n']);
            uprev_length = dummy.uprev_in_x0;
            ctrl.details.uprev_in_x0 = 1;
            
        elseif isfield(vars, 'uprev_in_x0'),
            ctrl.details.uprev_in_x0 = 1;
            
        end
        
        if isfield(vars, 'ref'),
            % we have the reference as a parametric variable
            fprintf('For tracking to work, you will need to provide\n');
            fprintf('"Options.reference" when calling "u=ctrl(x0, Options)".\n');
            parametric_vars = [parametric_vars; vars.ref];
            reference_length = length(vars.ref);
            
        elseif isfield(dummy, 'reference_in_x0'),
            % reference is already part of the state vector, hence parametric_vars
            % don't need to be updated
            reference_length = dummy.reference_in_x0;
            ctrl.details.reference_in_x0 = 1;
            
        elseif isfield(vars, 'reference_in_x0'),
            ctrl.details.reference_in_x0 = 1;
            
        end
        
        % define which variables should be optimization variables
        requested_vars = cat(1, vars.u{:});

        Fnonlin = any(~is(F, 'linear'));  % constraints are not linear
        if isfield(sysStruct, 'nonlinhandle') | Fnonlin,
            % we have a non-linear model, we must not convert constraints into
            % MPT matrices, because it wouldn't work. instead, we keep the
            % constraints, objective and variables as they are and later solve
            % the problem using solvesdp() in mpt_getInput().
            if ~isfield(sysStruct, 'nonlinhandle'),
                fprintf('WARNING: Setup contains nonlinear constraints.\n');
            end
            
            % expand the model. two reasons for doing that:
            %  1. expansion is expensive and hence we save time by not doing it
            %     at every step in mpt_getInput()
            %  2. YALMIP as of 16-May-2006 cannot save and load models which
            %     involve operators like iff() or implies(). however we can
            %     save/load the expanded model
            [Fexp, failure, cause] = expandmodel(F, obj, mptOptions.sdpsettings);
            if failure,
                fprintf('\n%s\n\n', cause);
                error('Cannot deal with given setup, see message above.');
            end

            ctrl.details.yalmipData = struct('F', Fexp, 'obj', obj, 'vars', vars, ...
                'uprev_length', uprev_length, 'reference_length', reference_length, ...
                'model_expanded', 1);
            
        else
            % expand the YALMIP model into MPT format
            M = mpt_yalmip2mpt(F, obj, parametric_vars, requested_vars);
            M.uprev_length = uprev_length;
            M.reference_length = reference_length;

            % remove equality constraints if needed (default for LTI systems to
            % improve performance of LP/QP solvers). however, according to
            % Johan, we should only do so if we have no binary variables
            justLTI = isempty(M.binary_var_index);
            
            if (justLTI & Options.remove_equalities_lti) | ...
                    Options.remove_equalities_always,
                % remove equalities
                if Options.verbose > 1,
                    fprintf('Projecting on equalities....\n');
                end
                M = mpt_collect_equalities(M, []);
                M = mpt_remove_equalities(M, []);
                M = mpt_project_on_equality(M);
            end
            
            ctrl.details.yalmipMatrices = M;
        end
         
    end
    
    if iscell(sysStruct.A) & Options.force_mpt == 1,
    
        % =========================================================================
        % we have a PWA system, check if MLD model is available in sysStruct
        
        haveMLD = 0;
        if isfield(sysStruct, 'data')
            if isfield(sysStruct.data, 'MLD')
                haveMLD = 1;
            end
        end
        if ~haveMLD
            % =========================================================================
            % no MLD model available, obtain it by doing a PWA2MLD conversion
            if Options.verbose > -1,
                disp('======================================================================');
                disp('WARNING: No MLD model available in system structure! Trying to convert');
                disp('         PWA model to MLD representation. This conversion can be very');
                disp('         inefficient! You should always opt for HYSDEL model instead.');
                disp('======================================================================');
                fprintf('\n'); 
            end
            if Options.verbose > 0,
                fprintf('Converting PWA system into MLD representation...\n');
            end
            if ~isfield(sysStruct, 'xmax') | ~isfield(sysStruct, 'xmin')
                % xmax and xmin must be defined for pwa2mld translation
                error('State constraints must be defined!');
            end
            
            % do the conversion...
            S = mpt_pwa2mld(sysStruct);
            
            % store MLD matrices to sysStruct...
            sysStruct.data.MLD = S;
            ctrl.sysStruct.data.MLD = S;
        end
        
        if isfield(probStruct, 'Rdu'),
            fprintf('Warning: Penalties on deltaU not supported for this type of controllers.\n');
        end
        [nx,nu,ny] = mpt_sysStructInfo(sysStruct);
        sysStruct.dims = struct('nx', nx, 'nu', nu, 'ny', ny);
        
        % try to construct matrices of the MPC problem
        % NOTE! not possible if probStruct.tracking > 0, in this case the
        % subfunction will return an empty matrix, which is fine (see
        % mpt_getInput and mpt_mip)
        if Options.verbose > 0,
            disp('Constructing data for on-line computation...');
        end
        ctrl.details.Matrices = sub_getMLDmatrices(sysStruct, probStruct, Options);
        
    elseif Options.force_mpt == 1,
        % =========================================================================
        % we have an LTI system, construct matrices

        % first augment matrices for tracking if necessary
        if probStruct.tracking==2 & (any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax)))
            % if we have deltaU constraints and tracking without deltaU formulation is
            % enabled, we use tracking with deltaU formulation which satisfies both
            % needs
            probStruct.tracking = 1;
        end
        
        if (probStruct.tracking | isfield(probStruct,'xref') | isfield(probStruct,'uref')) & ...
                ~iscell(sysStruct.A)
            % if tracking is requested, augment system and problem matrices
            [sysStruct, probStruct] = mpt_prepareTracking(sysStruct, probStruct);
        end
        
        useDUmode = 0;
        if isfield(probStruct, 'Rdu'),
            if any(probStruct.Rdu~=0),
                useDUmode = 1;
            end
        end
        if (any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax)))
            useDUmode = 1;
        end
        if useDUmode & ~probStruct.tracking,
            % augment state vector for deltaU constraints case to guarantee fullfilment
            % of those constraints in closed loop
            [sysStruct, probStruct] = mpt_prepareDU(sysStruct, probStruct);
        end

        if nargin==3,
            Options = varargin{end};
        else
            Options = [];
        end
        if ~isfield(Options, 'noConstraintReduction'),
            Options.noConstraintReduction = 1;
        end
        if ~isfield(Options, 'verbose'),
            Options.verbose = mptOptions.verbose;
        end
        if Options.verbose > -1,
            disp('Constructing data for on-line computation...');
        end
        if (isfield(probStruct,'inputblocking') | isfield(probStruct,'deltablocking'))
            opt = Options;
            opt.noConstraintReduction = 1;
            % do not reduce constraints yet, it will be done (if requested) in
            % mpt_blockingMatrices
            ctrl.details.Matrices = mpt_constructMatrices(sysStruct, probStruct, opt);
            ctrl.details.Matrices = mpt_blockingMatrices(ctrl.details.Matrices,sysStruct,probStruct, Options);
        else
            ctrl.details.Matrices = mpt_constructMatrices(sysStruct, probStruct, Options);
        end

    end
    ctrl.sysStruct = sysStruct;
    ctrl.probStruct = probStruct;

    ctrl.type = 'online';
else

    error('MPTCTRL: unknown type of input argument(s)!');
end


% =========================================================================
% store dimensions
if ~isfield(ctrl.details, 'dims'),
    [nx,nu,ny] = mpt_sysStructInfo(ctrl.sysStruct);
    ctrl.details.dims = struct('nx', nx, 'nu', nu, 'ny', ny);
end


% =========================================================================
% set format of the state vector
if ~isfield(ctrl.details, 'x0format'),
    ctrl = sub_setX0format(ctrl);
end


% =========================================================================
% create the object
ctrl = class(ctrl, 'mptctrl');



%--------------------------------------------------------------------------
function Matrices = sub_getMLDmatrices(sysStruct, probStruct, Options),

global mptOptions

verb = Options.verbose;
Options.verbose=0;
if probStruct.tracking > 0,
    % cannot pre-compute matrices of the MPC problem for time-varying references
    Matrices = [];
    return;
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

if isfield(probStruct, 'xN'),
    clear ref
    ref.x = probStruct.xN;
    % tell mpc_mip to include terminal set constraint
    Options.TerminalConstraint = 1;
    % zero tolerance on satisfaction of terminal set constraint
    Options.eps2 = mptOptions.abs_tol;
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

% tell mpc_mip() not to solve the given problem, rather just to return matrices
% of the MPC problem
Options.returnproblem = 1;

% dummy initial state, we are not solving the problem, so it does not matter
% what the value is
x0 = zeros(nx, 0);

Matrices = mpc_mip(S, x0, ref, weights, horizon, horizon, Options);
Matrices.F1eq = [];
Matrices.F2eq = [];
Matrices.F3eq = [];

Matrices.info = [];
if isfield(Options, 'convert2eq'),
    if Options.convert2eq==1,
        Matrices.info.eqconverted = 1;
        [ne, nx_f1] = size(Matrices.F1);
        [Ain,Bin,Aeq,Beq] = mpt_ineq2eq([Matrices.F1 -Matrices.F3], Matrices.F2);
        if ~isempty(Aeq),
            Matrices.F1 = Ain(:, 1:nx_f1);
            Matrices.F2 = Bin;
            Matrices.F3 = -Ain(:, nx_f1+1:end);
            Matrices.F1eq = Aeq(:, 1:nx_f1);
            Matrices.F2eq = Beq;
            Matrices.F3eq = -Aeq(:, nx_f1+1:end);
            if verb > 0,
                fprintf('%d inequalities originally / %d inequalities replaced by %d equalities / %d inequalities now\n', ...
                    ne, 2*length(Beq), length(Beq), length(Bin));
            end
            Matrices.info.ineqorig = ne;
            Matrices.info.neweq = length(Beq);
            Matrices.info.newineq = size(Ain, 1);
        end
    end
end

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


%---------------------------------------------------
function ctrl = sub_setX0format(ctrl)
% determines whether the state vector was augmented to deal with tracking/deltaU
% formulation. sets "ctrl.details.x0format" to a structure with following
% fields:
%
%  x0format.required   - number of states required
%  x0format.reference  - number of references (0 if the reference is not part of
%                        the state vector)
%  x0format.uprev      - number of inputs (0 if the previous control input is
%                        not part of the state vector)
%
% This data are then later used in extendx0()

sysStruct = ctrl.sysStruct;
probStruct = ctrl.probStruct;
dims = ctrl.details.dims;
uprev_in_x0 = isfield(ctrl.details, 'uprev_in_x0');
reference_in_x0 = isfield(ctrl.details, 'reference_in_x0');

% is u(k-1) included in the state vector?
include_uprev = (probStruct.tracking==0 & isfield(sysStruct, 'dumode')) | ...
    (probStruct.tracking==1 & (isfield(probStruct, 'tracking_augmented') | ...
    isfulldim(ctrl.Pn))) | uprev_in_x0;

% is the reference included in the state vector?
include_ref = reference_in_x0 | ((probStruct.tracking > 0) & ...
    (isfield(probStruct, 'tracking_augmented') | isfulldim(ctrl.Pn)));

[nx, nu, ny] = mpt_sysStructInfo(ctrl.details.origSysStruct);

if include_ref==0,
    nref = 0;
elseif isfield(ctrl.probStruct, 'Qy'),
    nref = ny;
else
    nref = nx;
end

if include_uprev==0,
    nu = 0;
end

ctrl.details.x0format = struct('required', nx + nu + nref, ...
    'reference', nref, ...
    'uprev', nu );

% remove fields which are no longer needed
if uprev_in_x0,
    ctrl.details = rmfield(ctrl.details, 'uprev_in_x0');
end
if reference_in_x0,
    ctrl.details = rmfield(ctrl.details, 'reference_in_x0');
end
