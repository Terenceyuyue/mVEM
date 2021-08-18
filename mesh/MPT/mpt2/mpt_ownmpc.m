function varargout = mpt_ownmpc(varargin)
% MPT_OWNMPC The "Design your own MPC" function
%
% [CON, OBJ, VAR] = mpt_ownmpc(sysStruct, probStruct)
% ctrl = mpt_ownmpc(sysStruct, probStruct, CON, OBJ, VAR)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% "Design Your Own MPC problem". This function operates in two modes:
% 
% 1. Problem construction phase:
%      In this step, matrices defining constraints and objective of a given MPC
%      problem are formulated. The matrices, together with variables which
%      define them, are stored to your workspace as CON (constraints), OBJ
%      (objective) and VAR (variables) objects.
%
% 2. Computation phase:
%      In this step, a control law is calculated according to constraints and
%      objective provided.
%
% Example:
%   We would like to impose polytopic constraints on all predicted states, but
%   not on the initial condition x0:
%
%     [H, K] = double(unitbox(2, 2));   % polytopic constraints
%     Double_Integrator
%     [CON, OBJ, VAR] = mpt_ownmpc(sysStruct, probStruct);
%     for k = 2:length(VAR.x)
%       % k==1 corresponds to x0, k==2 to x1 etc.
%       CON = CON + set(H*VAR.x{k} <= K);
%     end
%     ctrl = mpt_ownmpc(sysStruct, probStruct, CON, OBJ, VAR)
%
% Note!
%   To design an on-line MPC controller, you must append an additional 'online'
%   flag when calling this function, i.e.:
%
%     [C, O, V] = mpt_ownmpc(sysStruct, probStruct, 'online')
%     ctrl = mpt_ownmpc(sysStruct, probStruct, C, O, V, 'online')
%

% Copyright is with the following author(s):
%
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end
Options = [];
Options = mpt_defaultOptions(Options, ...
    'verbose', mptOptions.verbose, ...
    'dont_solve', 1, ...
    'dp', 0);

% by default we compute an explicit controller
ctrltype = '';
designphase = 1;
sysStruct = varargin{1};
probStruct = varargin{2};
if nargin > 2,
    if ischar(varargin{3}),
        if nargin ~= 3,
            error('Wrong number of input arguments.');
        end
        ctrltype = varargin{3};
    else
        if nargin < 5,
            error('Wrong number of input arguments.');
        end
        F = varargin{3};
        obj = varargin{4};
        vars = varargin{5};
        designphase = 0;
    end
end
if nargin==6,
    ctrltype = varargin{6};
    designphase = 0;
end

ctrl = [];
if designphase,
    % generate constraints, objective and variables

    if ~(iscell(sysStruct) | isstruct(sysStruct)),
        error('First input must be a structure.');
    end
    if ~isstruct(probStruct),
        error('Second input must be a structure.');
    end

    if isequal(lower(ctrltype), 'online') | isequal(lower(ctrltype), 'on-line'),
        % we need the following special setting for on-line controllers
        Options.yalmip_online = 1;
        Options.dont_solve = 1;
        Options.dp = 0;
    end

    Options.ownmpc = 1;
    try
        [dummy, C, O, V] = mpt_yalmipcftoc(sysStruct, probStruct, Options);
    catch
        error(lasterr);
    end
    
    if nargout==0,
        % assign varibales in caller's workspace
        
        % check for existence of certain variables in caller's workspace
        checkvars = {'CON', 'OBJ', 'VAR'};
        for ic = 1:length(checkvars),
            evalcmd = sprintf('exist(''%s'', ''var'')', checkvars{ic});
            exists = evalin('caller', evalcmd);
            if exists,
                fprintf('Variable "%s" was overwritten...\n', checkvars{ic});
            end
        end
        
        assignin('caller', 'CON', C);
        assignin('caller', 'OBJ', O);
        assignin('caller', 'VAR', V);
        fprintf('Variables CON, OBJ and VAR are now stored in your workspace.\n');
        
    elseif nargout==3,
        % return constraints, objective and variables
        varargout{1} = C;
        varargout{2} = O;
        varargout{3} = V;
        
    else
        error('Wrong number of output arguments.');
        
    end
    
else
    % compute controller based on custom constraints/objective
    
    if ~(isa(F, 'lmi') | isa(F, 'set')),
        error('Third input must be a set of constraints.');
    end
    if ~(isa(obj, 'sdpvar') | isa(obj, 'double'))
        error('Fourth input must be an optimization objective.');
    end
    if ~isstruct(vars),
        error('Fifth input must be a structure.');
    end
    if ~isfield(vars, 'x') | ~isfield(vars, 'u') | ~isfield(vars, 'y'),
        error('Wrong type of fifth input argument.');
    end
    if any(is(F, 'sigmonial')),
        error('Sigmonial constraints (e.g. 1/x <= a) not supported.');
    end
    if isa(obj, 'sdpvar') & is(obj, 'sigmonial'),
        error('Sigmonial objectives not supported.');
    end
    
    
    % =================================================================
    % verify sysStruct and probStruct
    verOpt.ybounds_optional = 1;
    verOpt.verbose = -1;    
    verOpt.useyalmip = 1; % to tell mpt_verifyProbStruct we can deal with move blocking
    userSysStruct = sysStruct;
    if iscell(userSysStruct),
        sysStruct = userSysStruct{1};
    end
    [sysStruct, probStruct] = mpt_verifySysProb(sysStruct, probStruct, verOpt);
    origSysStruct = sysStruct;
    origProbStruct = probStruct;

    
    % =================================================================
    % look whether we have to generate an on-line controller
    if isempty(ctrltype),
        ctrltype = mpt_defaultField(vars, 'type', 'explicit');
    else
        ctrltype = lower(ctrltype);
    end
    if isequal(lower(ctrltype), 'online') | isequal(lower(ctrltype), 'on-line'),
        % we need the following special setting for on-line controllers
        Options.yalmip_data.constraints = F;
        Options.yalmip_data.objective = obj;
        Options.yalmip_data.variables = vars;
        ctrl = mptctrl(userSysStruct, probStruct, Options);
        varargout{1} = ctrl;
        return
    end
    
    
    % ============================================================================
    % no explicit solution for non-linear system
    if isfield(sysStruct, 'nonlinhandle'),
        error('No explicit solution available for non-linear systems.');
    end

    
    % =================================================================
    % augment system to deal with tracking (on-line controllers use Uprev
    % and the reference directly as additional variables, hence we must not do
    % this augmentation for such case)
    if probStruct.tracking > 0 & ~isfield(probStruct, 'tracking_augmented')
        [sysStruct, probStruct] = mpt_yalmipTracking(sysStruct, probStruct, verOpt);
    end

    
    % =================================================================
    % display information about a given model
    if Options.verbose > 0,
        % expand the model to see how many binary variables we have
        Matrices = mpt_yalmip2mpt(F, obj, vars.x{1}, cat(1, vars.u{:}));
        nbinary = length(mpt_defaultField(Matrices, 'binary_var_index', []));
        nparam = length(mpt_defaultField(Matrices, 'param_var', []));
        noptim = length(mpt_defaultField(Matrices, 'requested_variables', []));
        isqp = mpt_defaultField(Matrices, 'qp', 0);
        
        sol_type = 'mp';
        if nbinary > 0,
            sol_type = [sol_type, 'MI'];
        end
        if isqp,
            sol_type = [sol_type 'QP'];
        else
            sol_type = [sol_type 'LP'];
        end
        fprintf('Solving an %s problem with %d parametric, %d optimization and %d binary variables\n\n', ...
                sol_type, nparam, noptim, nbinary);
    end
    
    
    % =================================================================
    % compute the controller
    yalmipOptions = mptOptions.sdpsettings;
    yalmipOptions.debug = 1;
    yalmipOptions.verbose = 1;
    yalmipOptions.mp.algorithm = 3;
    starttime = cputime;
    [sol, diagnost, Uz] = solvemp(F, obj, yalmipOptions, vars.x{1}, cat(1, vars.u{:}));
    if diagnost.problem ~= 0,
        fprintf('%s\n', diagnost.info);
        error('mpt_ownmpc: an error has occurred, see message above.');
    end
    
    if length(sol)==1,
        ctrl = sol{1};
        overlaps = 0;
        
    else
        if isempty(sol),
            ctrl = mptctrl;
            fprintf('The problem is infeasible:\n');
            diagnost.info
        end
        
        if probStruct.norm==2,
            ctrl = mpt_mergeCS(sol);
            overlaps = 1;
        else
            ctrl = mpt_removeOverlaps(sol);
            overlaps = 0;
        end
    end

    if isempty(sol),
        error('Problem is infeasible or an error occurred.');
    end
    if iscell(sol),
        if isempty(sol{1}),
            error('Problem is infeasible or an error occurred.');
        end
    end
    
    % =================================================================
    % set necessary fields of the controller structure
    ctrl.details.runTime = cputime - starttime;
    nx = length(ctrl.Bi{1});
    empty_Ais = find(cellfun('isempty', ctrl.Ai));
    if ~isempty(empty_Ais)
        [ctrl.Ai{empty_Ais}] = deal(zeros(nx));
    end

    ctrl.overlaps = overlaps;
    ctrl.details.origSysStruct = origSysStruct;
    ctrl.details.origProbStruct = origProbStruct;
    ctrl.sysStruct = sysStruct;
    ctrl.probStruct = probStruct;
    
    % copy data if tracking was used
    if isfield(vars, 'uprev_in_x0'),
        ctrl.details.uprev_in_x0 = 1;
    end
    if isfield(vars, 'reference_in_x0'),
        ctrl.details.reference_in_x0 = 1;
    end

    try
        ctrl = mptctrl(ctrl);
    catch
        fprintf('Error while creating controller object, result returned as a structure\n');
    end

    varargout{1} = ctrl;
    
end
