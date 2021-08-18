function [ctrl, F, obj, variables] = mpt_yalmipcftoc(sysStruct, probStruct, Options)
% MPT_YALMIPCFTOC CFTOC of PWA systems
%
% ctrl = mpt_yalmipcftoc(sysStruct, probStruct)
% ctrl = mpt_yalmipcftoc(sysStruct, probStruct, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the explicit solution of a given CFTOC problem with either linear
% or quadratic cost function. Supports PWA and LTI systems with
% boolean/integer/alphabet inputs (sysStruct.Uset), time-varying penalties, soft
% constraints (probStruct.{S|Sx|Su|Sy}), move blocking (probStruct.Nc), terminal
% state constraints (probStruct.xN), multiple target sets (probStruct.Tset),
% computes stabilizing sets for 2-norm problems.
%
% There a two approaches:
%  1. one-shot formulation
%  2. dynamic programming (DP) formulation
%
% You can choose the DP formulation by setting Options.dp=1
%
% By default we use the DP formulation if the cost function is linear, otherwise
% we use the one-shot formulation.
%
% Following problem formulations are supported:
%  * terminal state constraints - when probStruct.xN is given, terminal state
%                                 constraint is added (x(N)==probStruct.xN)
%  * time-varying penalties     - when probStruct.Q, probStruct.R, probStruct.Qy
%                                 are given as cell arrays
%  * multi-model dynamics       - allows to specify one model per prediction
%                                 step. to use this feature, specify sysStruct
%                                 as a cell array of system structures.
%  * move blocking              - if probStruct.Nc is given, first Nc control
%                                 moves are considered free (i.e.
%                                 u_0...u_(Nc-1)), while all  subsequent control
%                                 moves (u_Nc...u_(N-1)) are equal to u_(Nc-1).
%                                 works also for PWA and MLD systems!
%  * soft constraints - if "probStruct.Sx" is specified, uses soft state
%                       constraints.
%                     - probStruct.sxmax - maximum allowed violation of each
%                       state constraint. e.g. probStruct.sxmax=[10;0] allows to
%                       violate 1st state constraint by 10, while keeping the
%                       second constraint hard (i.e. no violation are allowed)
%                     - probStruct.Sy, probStruct.symax - softening of output
%                       constraints
%                     - probStruct.Su, probStruct.sumax - softening of input
%                       constraints
%  * penalize switching of PWA dynamics - probStruct.Qdyn, probStruct.Qswitch
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct         - System structure in the sysStruct format
% probStruct        - Problem structure in the probStruct format
% Options           - Additional options
%    .verbose       - level of verbosity
%    .mp.algorithm  - which enumeration algorithm to use
%                     (Options.mp_algorithm=3 uses different enumeration)
%    .mp.presolve   - perform pre-solving if this option is true (default)
%    .dp            - if set to 1, use the DP formulation
%    .force_mld     - force MLD dynamics (mainly for explicit controllers)
%    .force_pwa     - force PWA dynamics (mainly for on-line controllers)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% ctrl              - MPTCTRL object 
%
% see also MPT_CONTROL, MPT_OPTCONTROLPWA
%

% Copyright is with the following author(s):
%
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2006 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%          loefberg@control.ee.ethz.ch

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

% This function is the most general control function in the whole history of
% MPT. All the credit goes to Johan and his YALMIP, which makes it very easy to
% solve the CFTOC problem for LTI and PWA systems with both continuous and
% discrete inputs. The problem can be solved either in the dynamic programming
% fashion (default for linear cost functions), or in the one-shot formulation
% (default for quadratic cost functions). The user always force one of the
% formulation by setting "Options.dp" to a proper value.
%
% NOTE! if we have probStruct.xref/probStruct.uref, mpt_control() already shifts
% the dynamics such that the reference becomes the new origin. Therefore we
% always set the state/input reference to zero in the sequel. However, if you
% intend to call this function not via mpt_control, look for appropriate lines
% in this file and follow the instructions.

% NOTE! if sysStruct contains both the PWA and MLD representations, choice of
% which dynamics will be used to create the CFTOC model depends on which control
% strategy is used. For on-line controllers we currently prefere MLD models
% (they are much faster to solve than PWA models), while for explicit
% controllers we prefer to use the PWA representation. However, it is possible
% to force using MLD or PWA representation by setting either "Options.force_mld"
% or "Options.force_pwa" to 1.

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3,
    Options = [];
end
Options = mpt_defaultOptions(Options, ...
    'verbose', mptOptions.verbose, ...
    'includeLQRset', 1, ...
    'pwa_index', [], ...
    'dont_solve', 0, ...
    'details', 0, ...
    'force_pwa', 0, ...
    'force_mld', 0, ...
    'ownmpc', 0, ...
    'yalmip_online', 0);


%===============================================================================
% set options which we pass to solvemp()
yalmipOptions = mptOptions.sdpsettings;
yalmipOptions.verbose = 0;  % to keep intemediate mpLPs silent
f = fields(Options);
for ii = 1:length(f),
    yalmipOptions = setfield(yalmipOptions, f{ii}, getfield(Options, f{ii}));
end

% we need to check that prediction horizon is not infinity already at this stage
if ~isfield(probStruct, 'N'),
    error('Prediction horizon must be given.');
end
if isinf(probStruct.N),
    error('Prediction horizon must not be infinity.');
end


%===============================================================================
% handle multi-model systems (one sysStruct per one prediction step)
% part 1.
if iscell(sysStruct),
    multi_model = 1;
    if (length(sysStruct) ~= probStruct.N),
        error('Number of models must be equal to the prediction horizon.');
    else
        SST = sysStruct;        
    end
    SST{end+1} = SST{end};
    origSysStruct = sysStruct{1};
else
    multi_model = 0;
    origSysStruct = sysStruct;
    SST = sub_sysStruct2cell(sysStruct, probStruct.N+1);
end


%===============================================================================
% verify sysStruct and probStruct structures

% to keep mpt_verify*() silent
verOpt.verbose = -1;
% to tell mpt_verifySysStruct() that output constraints are optional
verOpt.ybounds_optional = 1;
% to tell mpt_verifyProbStruct() we can deal with move blocking
verOpt.useyalmip = 1; 

[origSysStruct, origProbStruct] = mpt_verifySysProb(origSysStruct, ...
    probStruct, verOpt);
pst = origProbStruct;

for ii = 1:length(SST),
    if ~isfield(SST{ii}, 'verified'),
        [SST{ii}, pst] = mpt_verifySysProb(SST{ii}, origProbStruct, verOpt);
    end
end

sysStruct = SST{1};


%===============================================================================
% augment the state vector to cope with tracking/deltaU formulation
[SST, probStruct, tracking_data] = sub_augmentx0(SST, pst, Options.yalmip_online);
sysStruct = SST{1};
pst = probStruct;
x0_format = 'x0';
if tracking_data.uprev_in_x0,
    x0_format = [x0_format '; u(k-1)'];
end
if tracking_data.reference_in_x0,
    x0_format = [x0_format '; reference'];
end
if Options.ownmpc & Options.verbose > 0 & ~isequal(x0_format, 'x0'),
    fprintf('State vector will have the format x = [%s].\n', x0_format);
end


%===============================================================================
% check system dimensions
[nx,nu,ny,nPWA,nbool,ubool] = mpt_sysStructInfo(sysStruct);


%===============================================================================
% do we have an LTI, PWA or an MLD dynamics?
[dynamics_type, nPWA, MLD, SST, sysStruct] = sub_checkmodels(SST, sysStruct, Options);

% do we have at least one MLD/PWA/LTI system?
haveMLD = ~isempty(findstr([dynamics_type{:}], 'mld'));
haveMLD3 = ~isempty(findstr([dynamics_type{:}], 'ml3'));
havePWA = ~isempty(findstr([dynamics_type{:}], 'pwa'));
haveLTI = ~isempty(findstr([dynamics_type{:}], 'lti'));
haveNONLIN = ~isempty(findstr([dynamics_type{:}], 'nonlin'));
havePNONLIN = ~isempty(findstr([dynamics_type{:}], 'pnonlin'));

% only on-line non-linear MPC is allowed
if haveNONLIN & Options.yalmip_online==0,
    error('Only on-line controllers can be designed for non-linear plants.');
end

% when formulating the MPC problem for PWA and piecewise nonlinear systems, we
% use the implies() operator, which requires all variables to be bounded. this
% flag therefore indicates later that we should introduce non-infty bounds on
% x,u and y variables to prevent warnings from YALMIP.
bigM_used = havePWA | havePNONLIN;


%===============================================================================
% check whether one can use the DP formulation
if ~isfield(Options, 'dp'),
    % use the DP formulation for linear cost, otherwise use the one-shot
    % formulation
    if probStruct.norm==2,
        Options.dp = 0;
    else
        Options.dp = 1;
    end
end
if Options.dp,
    if Options.dont_solve & probStruct.N>1,
        fprintf('WARNING: switching to one-shot formulation because Options.dont_solve is 1.\n');
        Options.dp = 0;
    elseif isfield(probStruct, 'Qswitch'),
        fprintf('WARNING: switching to one-shot formulation because probStruct.Qswitch is given.\n');
        Options.dp = 0;
    elseif Options.details == 1,
        fprintf('WARNING: switching to one-shot formulation because Options.details=1\n');
        Options.dp = 0;
    elseif haveNONLIN,
        fprintf('WARNING: switching to one-shot formulation because we have non-linear dynamics\n');
        Options.dp = 0;
    end
end


% set pwa_index
if ~isempty(Options.pwa_index),
    if ~(havePWA | haveLTI),
        error('Options.pwa_index can only be used for PWA systems.');
        
    else
        pwa_index = cell(1, length(SST));
        [pwa_index{:}] = deal(Options.pwa_index);
    end
    
else
    pwa_index = cell(1, length(nPWA));
    for im = 1:length(nPWA),
        pwa_index{im} = 1:nPWA(im);
    end
end


%===============================================================================
% reject certain problem descriptions
if probStruct.subopt_lev > 0,
    error('Suboptimal strategies not supported by this function.');
end
if isinf(probStruct.N),
    error('mpt_yalmipcftoc: Prediction horizon must be finite!');
end
if mpt_isnoise(sysStruct.noise),
    error('mpt_yalmipcftoc: additive noise not supported.');
end
if (nbool>0 & probStruct.tracking==1),
    error('mpt_yalmipcftoc: tracking cannot be used for systems with integer inputs. Set probStruct.tracking=2.');
end
if isfield(probStruct, 'Aunc') | isfield(probStruct, 'Bunc'),
    error('mpt_yalmipcftoc: parametric uncertainties are not supported.');
end


%===============================================================================
% handle control horizons
if isfield(probStruct, 'Nc'),
    Nc = probStruct.Nc;
    if Nc<probStruct.N & Options.dp==1,
        fprintf('Switching to one-shot formulation because control horizon is given.\n');
        Options.dp = 0;
    end
else
    Nc = probStruct.N;
end


%===============================================================================
% compute stabilizing target set and associated LQR cost-to-go for 2-norm
% problems

% we can't compute the stabilizing set if we have time-varying penalties,
% becuase mpt_computePWATset does not support them.
timevar_penalties = iscell(probStruct.Q) | iscell(probStruct.R);

% now try to compute the stabilizing target set. we can't do that if we have
% tracking, discrete inputs, penalties on otputs and/or own target set is
% already specified
if probStruct.norm==2 & probStruct.Tconstraint==1 & ...
        probStruct.tracking==0 & nbool == 0 & ~isfield(probStruct, 'Qy') & ...
        ~isfulldim(probStruct.Tset) & (havePWA | haveLTI) & multi_model==0,

    if timevar_penalties,
        fprintf('WARNING: No stabilizing target can be computed with time-varying penalties.\n');
        
    else
        if Options.verbose >= 0,
            fprintf('Adding stabilizing target set constraint...\n');
        end
        [Tset, P_N, origin_in_dyn] = sub_stabilizingset(SST{1}, pst, nx, nPWA(1));
        if isempty(origin_in_dyn),
            fprintf(['WARNING: No dynamics contains the origin in it''s interior, '...
                    'cannot compute a stabilizing target set.\n']);
            
        elseif isfulldim(Tset),
            probStruct.Tset = Tset;
            probStruct.P_N = P_N;
            
        else
            fprintf('WARNING: No stabilizing target set found!\n');
            
        end
    end
end


%===============================================================================
% create a cell array of probStruct.[R,Q,Qy] if they are not already provided as
% a cell array. this will be used later for simpler implementation of
% time-varying penalties
%
% Note that we must do this only after the stabilizing target set is computed,
% because mpt_computePWATset does not support time-varying penalties
%
probStruct = sub_timevarpenalties(probStruct);


%===============================================================================
% find out whether sysStruct.C, sysStruct.D and sysStruct.g all contain the same
% elements for all dynamics. if so, we can simplify things a bit
if multi_model | haveMLD | haveMLD3
    % always consider outputs as separate variables for the multi-model
    % approach, it makes things less messy coding-wise. do the same if we have
    % at least one MLD model.
    YeqSame = 0;
else
    YeqSame = sub_allYeqsEqual(SST);
end


%===============================================================================
% reject certain objectives if we have MLD/LTI systems
if ~havePWA,
    if isfield(probStruct, 'Qswitch'),
        error('probStruct.Qswitch can only be used for PWA systems.');
    elseif isfield(probStruct, 'Qdyn'),
        error('probStruct.Qdyn can only be used for PWA systems.');
    end
end


%===============================================================================
% check whether constraints are present
%
% haveXbounds(i) = 0    - no constraints
% haveXbounds(i) = 1    - constraints present (can be +/-Inf)
haveXbounds = zeros(1, length(SST));
haveYbounds = zeros(1, length(SST));
haveUbounds = zeros(1, length(SST));
dUconstraints = zeros(1, length(SST));
for im = 1:length(SST),
    haveXbounds(im) = isfield(SST{im}, 'xmax') & isfield(SST{im}, 'xmin');
    haveYbounds(im) = isfield(SST{im}, 'ymax') & isfield(SST{im}, 'ymin');
    haveUbounds(im) = isfield(SST{im}, 'umax') & isfield(SST{im}, 'umin');
    dUconstraints(im) = any(~isinf(SST{im}.dumin)) | any(~isinf(SST{im}.dumax)) | ...
        isfield(probStruct, 'Rdu');
end

% if bigM technique is used for PWA and Piecewise Nonlinear systems, we do
% require state and input constraints to be present in order to be able to make
% tight relaxations. If such constraints are not present, we give a warning and
% set them to +/- mptOptions.infbox (which is 1e4 by default)
noXbounds = find(haveXbounds==0);
noYbounds = find(haveYbounds==0);
noUbounds = find(haveUbounds==0);
if bigM_used & ~isempty(noXbounds),
    fprintf('WARNING: state constraints not given but required, setting them to +/- %d\n', ...
        mptOptions.infbox);
    SST = sub_force_bounds(SST, 'x', nx, noXbounds, mptOptions.infbox);
    haveXbounds(noXbounds) = 1;
end
if bigM_used & ~isempty(noYbounds),
    fprintf('WARNING: output constraints not given but required, setting them to +/- %d\n', ...
        mptOptions.infbox);
    SST = sub_force_bounds(SST, 'y', ny, noYbounds, mptOptions.infbox);
    haveYbounds(noYbounds) = 1;
end
if bigM_used & ~isempty(noUbounds),
    fprintf('WARNING: input constraints not given but required, setting them to +/- %d\n', ...
        mptOptions.infbox);
    SST = sub_force_bounds(SST, 'u', nu, noUbounds, mptOptions.infbox);
    haveUbounds(noUbounds) = 1;
end


%===============================================================================
% prepare references
yref = mpt_defaultField(probStruct, 'yref', zeros(ny, 1));
xref = mpt_defaultField(probStruct, 'xref', zeros(nx, 1));
uref = mpt_defaultField(probStruct, 'uref', zeros(nu, 1));
if haveMLD,
    % references for "d" and "z" variables of an MLD model
    nd = MLD{1}.nd; nz = MLD{1}.nz;
    dref = mpt_defaultField(probStruct, 'dref', zeros(nd, 1));
    zref = mpt_defaultField(probStruct, 'zref', zeros(nz, 1));
end
if haveMLD3
    nw = MLD{1}.nw;
    wref = mpt_defaultField(probStruct, 'wref', zeros(nw, 1));
end
if isfield(probStruct, 'xref_augmented') | probStruct.tracking>0,
    % this flag is set in mpt_prepareTracking and indicates that the system was
    % already augmented such that xref/uref becomes the new origin. in this case
    % we MUST set xref/uref to zero
    xref = 0*xref;
    uref = 0*uref;
end


%===============================================================================
if Options.yalmip_online,
    if Options.verbose > 0,
        if multi_model,
            fprintf('Using a time-varying model...\n');
        elseif haveMLD,
            fprintf('Using an MLD model...\n');
        elseif haveMLD3,
            fprintf('Using an MLD3 model...\n');            
        elseif havePWA,
            fprintf('Using a PWA model...\n');
        elseif haveNONLIN,
            fprintf('Using a nonlinear model...\n');
        else
            fprintf('Using an LTI model...\n');
        end
    end
end


%===============================================================================
% prepare necessary variables
N = probStruct.N + 1;

% States x(k), ..., x(k+N)
x = sub_cellsdpvar(nx, N);

% Inputs u(k), ..., u(k+N-1)
u = sub_cellsdpvar(nu, N-1);

% Outputs y(k), ..., y(k+N-1)
y = sub_cellsdpvar(ny, N-1);

d = cell(1, N-1);
z = cell(1, N-1);
for im = 1:length(SST)-1,
    if isequal(dynamics_type{im}, 'mld'),
        % prepare "d" and "z" variables for MLD models
        nd = MLD{im}.nd; nz = MLD{im}.nz;
        d{im} = binvar(nd, 1);
        z{im} = sdpvar(nz, 1);
        
    elseif isequal(dynamics_type{im}, 'ml3'),
        % prepare auxiliary variables
        w{im} = sdpvar(MLD{im}.nw, 1);
        
    else
        % we just need binary variables for PWA selection
        d{im} = binvar(nPWA(im), 1);
        
    end
end

Fbounds = set([]);
reference_used = 0;
uprev_used = 0;
if Options.yalmip_online & tracking_data.need_uprev
    % we have either tracking or deltaU constraints/penalties. therefore we need
    % to introduce a new variable to denote the previous control input
    uprev = sdpvar(nu, 1);
    % add bounds on this variable for better scaling
    [umin, umax] = sub_getMaxBounds(SST, 'u', nu);
    tag = 'umin < u(k-1) < umax';
    Fbounds = Fbounds + set(umin < uprev < umax, tag);
    uprev_used = 1;
end
if Options.yalmip_online & tracking_data.need_reference,
    % we have a varying reference in tracking problems, introduce a new
    % variable to denote the reference
    if isfield(probStruct, 'Qy'),
        % reference has the dimension of outputs
        reference = sdpvar(ny, 1);
        % add bounds on this variable for better scaling
        [ymin, ymax] = sub_getMaxBounds(SST, 'y', ny);
        tag = 'ymin < reference < ymax';
        Fbounds = Fbounds + set(ymin < reference < ymax, tag);
        yref = reference;
        
    else
        % reference has the dimension of states
        reference = sdpvar(nx, 1);
        % add bounds on this variable for better scaling
        [xmin, xmax] = sub_getMaxBounds(SST, 'x', nx);
        tag = 'xmin < reference < xmax';
        Fbounds = Fbounds + set(xmin < reference < xmax, tag);
        xref = reference;
    end
    reference_used = 1;
end


if Options.yalmip_online,
    if ~isfield(SST{1}, 'ymax'),
        fprintf('WARNING: no output constraints defined, problem can be badly scaled...\n');
    end
    if ~isfield(SST{1}, 'xmax'),
        fprintf('WARNING: no state constraints defined, problem can be badly scaled...\n');
    end
end


%===============================================================================
% introduce slack variables for soft constraints if necessary
dims = struct('nx', nx, 'nu', nu, 'ny', ny);
[soften, slacks, smax, sweights] = sub_prepareslacks(probStruct, N, haveXbounds, haveYbounds, dims);


%===============================================================================
% initialize the solution
F = set([]);
if uprev_used | reference_used,
    % add any constraints defined previously on uprev and the reference
    F = F + Fbounds;
end
J{N} = 0;  % initial cost
obj = 0;
starttime = cputime;


%===============================================================================
% now formulate constraints and objective of the CFTOC problem. use either the
% DP approach or one-shot formulation (depending on "Options.dp")
for k = N-1:-1:1    

    
    %===============================================================================
    if Options.dp & Options.verbose > -1,
        fprintf('--- Step %d (out of %d) ---\n', N-k, probStruct.N);
    end
    if k > Nc,
        % block u_k, k=Nc+1..N to be equal to u_Nc
        ku = Nc;
    else
        % otherwise u_k is a free variable
        ku = k;
    end
    iN = k - 1; iNu = ku - 1;

    
    %===============================================================================
    if Options.dp,
        % in the DP formulation we formulate constraints for each step
        % separatelly. In the one-shot formulation, we add constraints for all
        % steps.
        F = set([]);
        if uprev_used | reference_used,
            % add any constraints defined previously on uprev and the reference
            F = F + Fbounds;
        end
        
        % in the DP formulation we start a new objective for each step, in the
        % one-shot formulation we sum up cost for each horizon.
        obj = 0;
    end

    
    %===============================================================================
    % impose bounds on x_0 if necessary
    if isfulldim(SST{k}.Pbnd) & k==1 & Options.yalmip_online==0,
        % only include Pbnd constraint for explicit solutions
        tag = sprintf('x_%d in Pbnd', iN);        
        F = F + set(ismember(x{k}, SST{k}.Pbnd), tag);
    end
    
    
    %===============================================================================
    % input constraints
    if k <= Nc,
        umin = SST{k}.umin - slacks.all{k} - slacks.u{k};
        umax = SST{k}.umax + slacks.all{k} + slacks.u{k};
        tag = sprintf('umin < u_%d < umax', iNu);
        if soften.u,
            tag = [tag ' (soft)'];
        end
        F = F + set(umin < u{ku} < umax, tag);
    end
   
    % input slew rate constraints
    if dUconstraints(k) & k <= Nc - 1,
        dumin = SST{k}.dumin - slacks.all{k} - slacks.u{k};
        dumax = SST{k}.dumax + slacks.all{k} + slacks.u{k};

        if uprev_used & k==1,
            % add constraint dumin < u(0)-u(-1) < dumax, such that we can
            % enforce deltaU constraints knowing the previous input (which is a
            % parametric variable)
            tag = 'dumin < u_0 - u_prev < dumax';
            du_var = u{k} - uprev;
            uprev_used = 1;
            
        else
            % add a constraint dumin < u{k+1} - u{k} < dumax
            tag = sprintf('dumin < u_%d - u_%d < dumax', iNu+1, iNu);
            du_var = u{ku+1} - u{ku};
            
        end
        F = F + set(dumin < du_var < dumax, tag);
    end
    
    % some inputs can be boolean or from finite alphabet
    for iu = 1:nu,
        [t, inputset] = sub_inputtype(SST{k}, iu);
        tag = sprintf('u_%d(%d) in %s', iNu, iu, mat2str(inputset));
        if t == 'B',
            % this input is boolean
            F = F + set(binary(u{ku}(iu)), tag);
        elseif t == 'A',
            % this input is from finite alphabet
            F = F + set(ismember(u{ku}(iu), inputset), tag);
        end
    end

    
    %===============================================================================
    % state constraints
    if haveXbounds(k),
        xmin = SST{k}.xmin - slacks.all{k} - slacks.x{k};
        xmax = SST{k}.xmax + slacks.all{k} + slacks.x{k};
        tag = sprintf('xmin < x_%d < xmax', iN);
        if soften.x,
            tag = [tag ' (soft)'];
        end        
        F = F + set(xmin < x{k} < xmax, tag);
    end
    if (k==N-1 | (Options.dp==1 & bigM_used)) & haveXbounds(k+1),
        % add state constraints on x_N or on x{k+1} if needed
        %
        % notice that if the dynamic programming approach is used and bigM
        % technique has to be used by means of the implies() operator, we must
        % also add constraints on x{k+1}, otherwise YALMIP would scream that the
        % problem is badly scaled
        
        xmin = SST{k+1}.xmin - slacks.all{k+1} - slacks.x{k+1};
        xmax = SST{k+1}.xmax + slacks.all{k+1} + slacks.x{k+1};
        tag = sprintf('xmin < x_%d < xmax', iN+1);
        if soften.x,
            tag = [tag ' (soft)'];
        end        
        F = F + set(xmin < x{k+1} < xmax, tag);
    end
        
    
    %===============================================================================
    % output constraints
    if haveYbounds(k) & (k > 1 | probStruct.y0bounds==1),
        % soften constraints if necessary
        ymin = SST{k}.ymin - slacks.all{k} - slacks.y{k};
        ymax = SST{k}.ymax + slacks.all{k} + slacks.y{k};
        tag = sprintf('ymin < y_%d < ymax', iN);
        if soften.y,
            tag = [tag ' (soft)'];
        end
        F = F + set(ymin < y{k} < ymax, tag);
    end
    
    
    %===============================================================================
    % add constraints on slacks (slacks have to be positive)
    if soften.all,
        tag = sprintf('0 < sa_%d < smax', iN);
        F = F + set(0 <= slacks.all{k} <= smax.all, tag);
    end
    if soften.x,
        tag = sprintf('0 < sx_%d < sxmax', iN);
        F = F + set(0 <= slacks.x{k} <= smax.x, tag);
    end
    if soften.u,
        tag = sprintf('0 < su_%d < sumax', iN);
        F = F + set(0 <= slacks.u{k} <= smax.u, tag);
    end
    if soften.y,
        tag = sprintf('0 < sy_%d < symax', iN);
        F = F + set(0 <= slacks.y{k} <= smax.y, tag);
    end
    
        
    %===============================================================================
    % Dynamics
    switch lower(dynamics_type{k}),
        case 'lti',
            % LTI dynamics. note that we have converted LTI systems into PWA form
            % with one dynamics
            tag = sprintf('x_%d == A*x_%d + B*u_%d', iN+1, iN, iNu);
            F = F + set(x{k+1} == SST{k}.A{1}*x{k} + SST{k}.B{1}*u{ku} + SST{k}.f{1}, tag);
            
            tag = sprintf('y_%d == C*x_%d + D*u_%d', iN, iN, iNu);
            F = F + set(y{k} == SST{k}.C{1}*x{k} + SST{k}.D{1}*u{ku}, tag);

        case 'nonlin',
            % non-linear time-invariant dynamics
            tag = sprintf('x_%d == %s(''state'', x_%d, u_%d)', ...
                iN+1, func2str(SST{k}.nonlinhandle), iN, iNu);
            xeq = feval(SST{k}.nonlinhandle, 'state', x{k}, u{ku});
            error(sub_check_nonlinearity(xeq));
            F = F + set(x{k+1} == xeq, tag);
            
            tag = sprintf('y_%d == %s(''output'', x_%d, u_%d)', ...
                iN, func2str(SST{k}.nonlinhandle), iN, iNu);
            yeq = feval(SST{k}.nonlinhandle, 'output', x{k}, u{ku});
            error(sub_check_nonlinearity(yeq));
            F = F + set(y{k} == yeq, tag);
            
        case 'pnonlin',
            % piecewise non-linear time-invariant dynamics
            for i = pwa_index{k}
                tag = sprintf('d_%d(%d) => x_%d == %s(''state'', x_%d, u_%d, %d)', ...
                    iN, i, iN+1, func2str(SST{k}.nonlinhandle), iN, iNu, i);
                xeq = feval(SST{k}.nonlinhandle, 'state', x{k}, u{ku}, i);
                error(sub_check_nonlinearity(xeq));
                F = F + set(implies(d{k}(i), x{k+1} == xeq), tag);
                
                tag = sprintf('d_%d(%d) => y_%d == %s(''output'', x_%d, u_%d, %d)', ...
                    iN, i, iN, func2str(SST{k}.nonlinhandle), iN, iNu, i);
                yeq = feval(SST{k}.nonlinhandle, 'output', x{k}, u{ku}, i);
                error(sub_check_nonlinearity(yeq));
                F = F + set(implies(d{k}(i), y{k} == yeq), tag);
                
                tag = sprintf('d_%d(%d) => %s(''guards'', x_%d, u_%d, %d)', ...
                    iN, i, func2str(SST{k}.nonlinhandle), iN, iNu, i);
                guardeq = feval(SST{k}.nonlinhandle, 'guards', x{k}, u{ku}, i);
                error(sub_check_nonlinearity(guardeq));
                
                % YALMIP would give a nasty warning if the guard expression is
                % not bounded by a aconstraint. however, it is currently not
                % clear how to come up with good bounds. just imagine our guard
                % is defined as follows:
                %   guardeq = (x(1)^2 + x(2)^2 <= 1)
                % even though all elements of the state vector are bounded,
                % YALMIP currently does not propagate such bounds into nonlinear
                % expressions. if we just add dummy bounds like
                %   set( -1e4 < guardeq < 1e4)
                % then it is NOT imposing +/-1e4 bounds on individual variables,
                % but instead on the whole nonlinear expression (which can be a
                % wrong thing to do). 
                %
                % therefore for now we don't assign any bounds, which will cause
                % a warning to be displayed by YALMIP. we use it as a reminder
                % that this should be improved
                F = F + set(implies(d{k}(i), guardeq), tag);
            end
            
            tag = sprintf('sum d_%d = 1', iN);
            F = F + set(sum(d{k}(pwa_index{k})) == 1, tag);
            
        case 'pwa',
            % PWA dynamics
            for i = pwa_index{k}
                tag = sprintf('d_%d(%d) => x_%d == A{%d}*x_%d + B{%d}*u_%d + f{%d}', ...
                    iN, i, iN+1, i, iN, i, iNu, i);
                F = F + set(implies(d{k}(i), x{k+1} == SST{k}.A{i}*x{k} + ...
                    SST{k}.B{i}*u{ku} + SST{k}.f{i}), tag);
                
                if ~YeqSame,
                    % we need to have the output as a variable
                    tag = sprintf('d_%d(%d) => y_%d == C{%d}*x_%d + D{%d}*u_%d + g{%d}', ...
                        iN, i, iN, i, iN, i, iNu, i);
                    F = F + set(implies(d{k}(i), y{k} == SST{k}.C{i}*x{k} + ...
                        SST{k}.D{i}*u{ku} + SST{k}.g{i}), tag);
                end
                
                tag = sprintf('d_%d(%d) => guardX{%d}*x_%d + guardU{%d}*u_%d <= guardC{%d}', ...
                    iN, i, i, iN, i, iNu, i);
                F = F + set(implies(d{k}(i),SST{k}.guardX{i}*x{k} + ...
                    SST{k}.guardU{i}*u{ku} <= SST{k}.guardC{i}), tag);
            end
            
            tag = sprintf('sum d_%d = 1', iN);
            F = F + set(sum(d{k}(pwa_index{k})) == 1, tag);
            
            if YeqSame,
                % all output equations are the same for all dynamics, it is
                % enough to define one output variable
                tag = sprintf('y_%d == C*x_%d + D*u_%d + g', iN, iN, iNu);
                F = F + set(y{k} == SST{k}.C{1}*x{k} + ...
                    SST{k}.D{1}*u{ku} + SST{k}.g{1}, tag);
            end
            
        case 'mld',
            % MLD dynamics
            %        x(k+1) = A x + B1 u + B2 d + B3 z
            %         y(k)  = C x + D1 u + D2 d + D3 z
            %  E2 d + E3 z <= E1 u + E4 x + E5
            tag = sprintf('x_%d == A*x_%d + B1*u_%d + B2*d_%d + B3*z_%d + B5', ...
                iN+1, iN, iNu, iN, iN);
            F = F + set(x{k+1} == MLD{k}.A*x{k} + MLD{k}.B1*u{ku} + ...
                MLD{k}.B2*d{k} + MLD{k}.B3*z{k} + MLD{k}.B5, tag);
            
            tag = sprintf('y_%d == C*x_%d + D1*u_%d + D2*d_%d + D3*z_%d + D5', ...
                iN, iN, iNu, iN, iN);
            F = F + set(  y{k} == MLD{k}.C*x{k} + MLD{k}.D1*u{ku} + ...
                MLD{k}.D2*d{k} + MLD{k}.D3*z{k} + MLD{k}.D5, tag);
            
            tag = sprintf('E2*d_%d + E3*z_%d <= E1*u_%d + E4*z_%d + E_5', ...
                iN, iN, iNu, iN);
            F = F + set(MLD{k}.E2*d{k} + MLD{k}.E3*z{k} <= MLD{k}.E1*u{ku} + ...
                MLD{k}.E4*x{k} + MLD{k}.E5, tag);
            
            % some states, inputs and/or outputs can be boolean:
            % states x(MLD.nx-MLD.nxb+1:MLD.nx) are boolean,
            % inputs u(MLD.nu-MLD.nub+1:MLD.nu) are boolean,
            % outputs y(MLD.ny-MLD.nyb+1:MLD.ny) are boolean
            ib1 = MLD{k}.nx-MLD{k}.nxb+1; ib2 = MLD{k}.nx; ib = ib1:ib2;
            if ~isempty(ib),
                tag = sprintf('x_%d(%d:%d) binary', iN, ib1, ib2);
                F = F + set(binary(x{k}(ib)), tag);
            end
            ib1 = MLD{k}.nu-MLD{k}.nub+1; ib2 = MLD{k}.nu; ib = ib1:ib2;
            if ~isempty(ib),
                tag = sprintf('u_%d(%d:%d) binary', iN, ib1, ib2);
                F = F + set(binary(u{k}(ib)), tag);
            end
            ib1 = MLD{k}.ny-MLD{k}.nyb+1; ib2 = MLD{k}.ny; ib = ib1:ib2;
            if ~isempty(ib),
                tag = sprintf('y_%d(%d:%d) binary', iN, ib1, ib2);
                F = F + set(binary(y{k}(ib)), tag);
            end
            
            % add bounds on auxiliary "z" variables
            tag = sprintf('MLD.zl < z_%d < MLD.zu', iN);
            F = F + set(MLD{k}.zl <= z{k} <= MLD{k}.zu, tag);

        case 'ml3',
            % MLD model by HYSDEL3
            %        x(k+1) = A x + Bu u + Baux w + Baff
            %         y(k)  = C x + Du u + Daux w + Daff
            %  Ex x + Eu u + Eaux w <= Eaff
            tag = sprintf('x_%d == A*x_%d + Bu*u_%d + Baux*w_%d + Baff', ...
                iN+1, iN, iNu, iN);
            F = F + set(x{k+1} == MLD{k}.A*x{k} + MLD{k}.Bu*u{ku} + ...
                MLD{k}.Baux*w{k} + MLD{k}.Baff, tag);
            
            tag = sprintf('y_%d == C*x_%d + Du*u_%d + Daux*w_%d + Daff', ...
                iN, iN, iNu, iN);
            F = F + set(  y{k} == MLD{k}.C*x{k} + MLD{k}.Du*u{ku} + ...
                MLD{k}.Daux*w{k} + MLD{k}.Daff, tag);

            ineq_idx = MLD{k}.j.ineq;
            if ~isempty(ineq_idx)
                tag = sprintf('Ex*x_%d + Eu*u_%d + Eaux*w_%d <= Eaff', ...
                    iN, iNu, iN);
                F = F + set(MLD{k}.Ex(ineq_idx, :)*x{k} + ...
                    MLD{k}.Eu(ineq_idx, :)*u{ku} + ...
                    MLD{k}.Eaux(ineq_idx, :)*w{k} <= ...
                    MLD{k}.Eaff(ineq_idx, :), tag);
            end

            eq_idx = MLD{k}.j.eq;
            if ~isempty(eq_idx)
                tag = sprintf('Ex*x_%d + Eu*u_%d + Eaux*w_%d == Eaff', ...
                    iN, iNu, iN);
                F = F + set(MLD{k}.Ex(eq_idx, :)*x{k} + ...
                    MLD{k}.Eu(eq_idx, :)*u{ku} + ...
                    MLD{k}.Eaux(eq_idx, :)*w{k} == ...
                    MLD{k}.Eaff(eq_idx, :), tag);
            end

            % some states, inputs and/or outputs can be boolean:
            % states x(MLD.j.xb) are boolean,
            % inputs u(MLD.j.ub) are boolean,
            % outputs y(MLD.j.yb) are boolean
            % auxiliaries w(MLD.j.d) are boolean
            bidx = MLD{k}.j.xb;
            if ~isempty(bidx)
                tag = sprintf('x_%d(%s) binary', iN, mat2str(bidx));
                F = F + set(binary(x{k}(bidx)), tag);
            end
            bidx = MLD{k}.j.ub;
            if ~isempty(bidx)
                tag = sprintf('u_%d(%s) binary', iNu, mat2str(bidx));
                F = F + set(binary(u{ku}(bidx)), tag);
            end
            bidx = MLD{k}.j.d;
            if ~isempty(bidx)
                tag = sprintf('w_%d(%s) binary', iN, mat2str(bidx));
                F = F + set(binary(w{k}(bidx)), tag);
            end
            ridx = MLD{k}.j.z;
            if ~isempty(ridx)
                % add bounds on auxiliary "z" variables
                tag = sprintf('MLD.wl < w_%d(%s) < MLD.wu', iN, mat2str(ridx));
                F = F + set(MLD{k}.wl(ridx) <= w{k}(ridx) <= MLD{k}.wu(ridx), tag);
            end

        otherwise
            error('Unknown type of system dynamics.');
            
    end
    

    %===============================================================================
    % add target set constraint
    if k==N-1 & isfulldim(probStruct.Tset),
        % ismember automatically handles cases where Tset is a polytope array
        tag = sprintf('x_%d in Tset', iN+1);
        F = F + set(ismember(x{k+1}, probStruct.Tset), tag);
    end


    %===============================================================================
    % add terminal state constraint
    if k==N-1 & isfield(probStruct, 'xN'),
        tag = sprintf('x_%d == probStruct.xN', iN+1);
        F = F + set(x{k+1} == probStruct.xN, tag);
    end
    
    
    %===============================================================================
    % objective function
    if isfield(probStruct, 'Qy'),
        % add penalty on outputs
        obj = obj + sub_norm((y{k} - yref), probStruct.Qy{k}, probStruct.norm);
        
    else
        % add penalty on states
        obj = obj + sub_norm((x{k} - xref), probStruct.Q{k}, probStruct.norm);
        
        % add terminal weight if specified
        if k==N-1 & isfield(probStruct, 'P_N'),
            obj = obj + sub_norm((x{k+1} - xref), probStruct.P_N, probStruct.norm);
        end
        
    end
    
    
    %===============================================================================
    % add penalty on inputs
    if 1 | (k <= Nc),
        % we could penalize just the inputs from k=0..Nc-1, i.e.
        %
        % min (\sum_{k=0}^{Nc-1} ||Q*x_k||_p + ||R*u_k||) + (\sum_{k=Nc}^{N-1} ||Q*x_k||_p)
        %
        % but currently we penalize
        % min \sum_{k=0}^{N-1} ||Q*x_k||_p + ||R*u_k||_p
        obj = obj + sub_norm((u{ku} - uref), probStruct.R{k}, probStruct.norm);
    end

    
    %===============================================================================
    % add penalty on deltaU
    if isfield(probStruct, 'Rdu') & k <= Nc - 1,
        obj = obj + sub_norm(du_var, probStruct.Rdu{ku}, probStruct.norm);
    end


    %===============================================================================
    % penalize delta and z variables from the MLD model
    if isequal(dynamics_type{k}, 'mld'),
        if isfield(probStruct, 'Qz'),
            obj = obj + sub_norm((z{k} - zref), probStruct.Qz{k}, probStruct.norm);
        end
        if isfield(probStruct, 'Qd'),
            obj = obj + sub_norm((d{k} - dref), probStruct.Qd{k}, probStruct.norm);
        end
    elseif isequal(dynamics_type{k}, 'ml3'),
        if isfield(probStruct, 'Qw')
            obj = obj + sub_norm((w{k} - wref), probStruct.Qw{k}, probStruct.norm);
        end
    end

    
    %===============================================================================
    % penalize switching between dynamics of a PWA system    
    if isequal(dynamics_type{k}, 'pwa') & isfield(probStruct, 'Qswitch'),
        obj = obj + sub_norm(d{k+1} - d{k}, probStruct.Qswitch{k}, probStruct.norm);
    end

    
    %===============================================================================
    % penalize "d" variables which denote to which dynamics does a state belong
    % to (only for PWA systems!)
    if isequal(dynamics_type{k}, 'pwa') & isfield(probStruct, 'Qdyn'),
        obj = obj + sub_norm(d{k}, probStruct.Qdyn{k}, probStruct.norm);
    end
    
    
    %===============================================================================
    % penalize slacks
    if soften.all,
        obj = obj + sub_norm(slacks.all{k}, sweights.all, probStruct.norm);
    end    
    if soften.x,
        obj = obj + sub_norm(slacks.x{k}, sweights.x, probStruct.norm);
    end    
    if soften.y,
        obj = obj + sub_norm(slacks.y{k}, sweights.y, probStruct.norm);
    end
    if soften.u,
        obj = obj + sub_norm(slacks.u{k}, sweights.u, probStruct.norm);
    end
        

    %===============================================================================
    if Options.dp,
        % Solve one-step problem, add cost-to-go J{k+1}
        
        % decrease verbosity value. by default we use Options.verbose=1, but
        % this makes mpt_mplp/mpt_mpqp display number of regions obtained in
        % each call, which can be rather disturbing.
        %yalmipOptions.verbose = Options.verbose - 1;
        [sol{k}, diagnost{k}, Uz{k}, J{k}] = solvemp(F, obj + J{k+1}, ...
            yalmipOptions, x{k}, u{k});
        if Options.verbose > -1 & ~isempty(sol{k}),
            nr = 0;
            for ii = 1:length(sol{k}),
                if ~isempty(sol{k}{ii}),
                    nr = nr+length(sol{k}{ii}.Pn);
                end
            end
            fprintf('Generated %d regions\n', nr);
        end
    end
    
end

variables.x = x;
variables.u = u(1:Nc);
variables.y = y;
if havePWA | haveMLD | havePNONLIN,
    variables.d = d;
end
if haveMLD3
    variables.w = w;
end
if haveMLD,
    variables.z = z;
end
if soften.x,
    variables.sx = slacks.x;
end
if soften.y,
    variables.sy = slacks.y;
end
if soften.u,
    variables.su = slacks.u(1:Nc);
end
if Options.yalmip_online,
    if reference_used,
        variables.ref = reference;
    end
    if uprev_used,
        variables.uprev = uprev;
    end
    % remember that this setup corresponds to an on-line controller    
    variables.type = 'online';
else
    variables.type = 'explicit';
end
if tracking_data.reference_in_x0,
    variables.reference_in_x0 = 1;
end
if tracking_data.uprev_in_x0,
    variables.uprev_in_x0 = 1;
end

if Options.dont_solve,
    % just return constraints, objectives and variables
    ctrl = [];
    ctrl.sysStruct = sysStruct;
    ctrl.probStruct = pst;
    if tracking_data.uprev_in_x0,
        % we did augment the state vector by the previous input
        ctrl.uprev_in_x0 = nu;
    end
    if tracking_data.reference_in_x0,
        % we did augment the state vector by the reference
        ctrl.reference_in_x0 = tracking_data.nref;
    end
    ctrl.dynamics_type = dynamics_type;
    return
end

if Options.dp & Options.verbose > -1 & probStruct.norm~=2,
    % make an empty line to separate output from mpt_removeOverlaps
    fprintf('\n');
end

%===============================================================================
% solve the one-shot problem as mpMIQP/mpMILP if necessary
if Options.dp==0,
    % obtain the open-loop solution:
    [sol{k}, diagnost{k}, Uz{k}] = solvemp(F, obj, yalmipOptions, x{k}, cat(1, u{:}));
end


%===============================================================================
% collect overlapping partitions together, remove overlaps if cost is linear
overlaps = 1;
if length(sol{k})==1,
    ctrl = sol{k}{1};
    overlaps = 0;
    
else
    if isempty(sol{k}),
        fprintf('\n%s\n', diagnost{k}.info);
        error('mpt_yalmipcftoc: an error has occurred, see message above.');
    end
        
    if probStruct.norm==2,
        ctrl = mpt_mergeCS(sol{k});
    else
        ctrl = mpt_removeOverlaps(sol{k});
        overlaps = 0;
    end
    
end


%===============================================================================
if isempty(ctrl),
    % problem either infeasible or some problem occured
    ctrl = mptctrl;
    return
end

if probStruct.norm~=2,
    for ii = 1:length(ctrl.Ai),
        ctrl.Ai{ii} = zeros(nx);
    end
end

% add necessary fields and create an MPTCTRL object
ctrl.overlaps = overlaps;
% if probStruct.norm==2,
%     % overlaps are possible with quadratic cost
%     ctrl.overlaps = 1;
% else
%     % overlaps have been removed for linear cost
%     ctrl.overlaps = 0;
% end

ctrl.sysStruct = SST{1};
ctrl.probStruct = pst;
ctrl.details.origSysStruct = origSysStruct;
ctrl.details.origProbStruct = origProbStruct;
ctrl.details.runTime = cputime - starttime;
try
    ctrl = mptctrl(ctrl);
end



%-----------------------------------------------------------------
function [t, s] = sub_inputtype(sysStruct, ind)
% t='C', s=[]                     if input 'ind' is continuous
% t='B', s=[0 1]                  if input 'ind' is boolean
% t='A', s=sysStruct.Uset{ind}    if input 'ind' is from an alphabet

if isfield(sysStruct, 'Uset'),
    if isequal(sort(sysStruct.Uset{ind}), [0 1]),
        % input is boolean
        t = 'B';
        s = [0 1];
        
    elseif all(isinf(sysStruct.Uset{ind})),
        % input is continuous
        t = 'C';
        s = [];
        
    else
        % input is from finite alphabet
        t = 'A';
        s = sysStruct.Uset{ind};
    end
    
else
    % all inputs are continuous
    t = 'C';
    s = [];
end


%-----------------------------------------------------------------
function yesno = sub_allYeqsEqual(sysStruct)
% returns true if sysStruct.C, sysStruct.D and sysStruct.g are identical for all
% dynamics

yesno = 1;
if ~iscell(sysStruct),
    sysStruct = {sysStruct};
end
C = sysStruct{1}.C{1};
D = sysStruct{1}.D{1};
g = sysStruct{1}.g{1};
for ii = 1:length(sysStruct)
    for jj = 1:length(sysStruct{ii}.C),
        if ~isequal(C, sysStruct{ii}.C{jj}) | ...
                ~isequal(D, sysStruct{ii}.D{jj}) | ~isequal(g, sysStruct{ii}.g{jj}),
            yesno = 0;
            return
        end
    end
end


%-----------------------------------------------------------------
function probStruct = sub_timevarpenalties(probStruct)
% if probStruct.[R,Q,Qy] are single matrices, we convert them to a cell array of
% length probStruct.N

penalties = {'R', 'Q', 'Rdu', 'Qy', 'Qd', 'Qz', 'Qdyn', 'Qswitch'};

for ii = 1:length(penalties),
    if isfield(probStruct, penalties{ii}),
        P = getfield(probStruct, penalties{ii});
        if ~iscell(P),
            Pcell = cell(1, probStruct.N);
            [Pcell{:}] = deal(P);
            probStruct = setfield(probStruct, penalties{ii}, Pcell);
        end
    end
end


%-----------------------------------------------------------------
function [Tset, P_N, origin_in_dynamics] = sub_stabilizingset(sysStruct, probStruct, nx, nPWA)
% computes a stabilizing set for a given PWA system (LTI systems are
% converted to PWA systems anyhow in the main function)

global mptOptions

Tset = polytope;
P_N = [];

% first find out which dynamics contain the origin in their interior
%
% here we could take xref and uref from probStruct.xref and probStruct.uref,
% respectively. but mpt_control already takes over this and performes change of
% coordinates by calling mpt_prepareTracking, therefore we must use zero
% reference in the sequel 
%
origin = zeros(nx, 1);

origin_in_dynamics = [];

for ii=1:nPWA
    if(all(sysStruct.f{ii}==0)) | isfield(probStruct,'xref') | isfield(probStruct,'uref')
        %.... othewise the origin would not be an equilibrium point
        if all(sysStruct.guardU{ii}==0) & ...
                max(sysStruct.guardX{ii}*origin - sysStruct.guardC{ii}) <= mptOptions.abs_tol,
            origin_in_dynamics = [origin_in_dynamics ii];
            
        elseif(any(sysStruct.guardU{ii}~=0))
            tempP=polytope(sysStruct.guardU{ii}, ...
                -sysStruct.guardX{ii}*origin + sysStruct.guardC{ii});
            
            if(isfulldim(tempP))
                origin_in_dynamics = [origin_in_dynamics ii];
            end
            
        end
    end
end
if isempty(origin_in_dynamics),
    return
end
Options.verbose = -1;
[Tset,Fi,Gi,dynamics,pst] = mpt_computePWATset(sysStruct, probStruct, origin_in_dynamics, Options);
P_N = pst.P_N;


%-----------------------------------------------------------------
function o = sub_norm(x, P, n)

if n == 2,
    o = x'*P*x;
else
    o = norm(P*x, n);
end


%-----------------------------------------------------------------
function SST = sub_sysStruct2cell(sysStruct, N)
% convert single sysStruct into a cell array

SST = cell(1, N);
[SST{:}] = deal(sysStruct);


%-----------------------------------------------------------------
function type = sub_sysStruct_type(sysStruct)
% this subfunction should be moved into mpt_sysStructInfo()

type = struct('lti', 0, 'pwa', 0, 'mld', 0, 'ml3', 0, 'nonlin', 0);
type.lti = ~iscell(sysStruct.A);
type.pwa = iscell(sysStruct.A);
if isfield(sysStruct, 'data')
    if isfield(sysStruct.data.MLD, 'nw')
        % HYSDEL 3 model
        type.ml3 = 1;
    else
        type.mld = 1;
    end
end
if type.ml3 & type.mld
    type.mld = 0;
end
if type.mld | type.ml3,
    if isfield(sysStruct.data, 'onlymld'),
        % the {A,B,C,D} touple is just fake, we in fact do have only MLD
        % representation available
        type.pwa = 0;
    end
end
type.nonlin = isfield(sysStruct, 'nonlinhandle');


%-----------------------------------------------------------------
function [dynamics_type, nPWA, MLD, SST, sysStruct] = sub_checkmodels(SST, sysStruct, Options)
% checks whether we have MLD, PWA or LTI systems

dynamics_type = cell(1, length(SST));
nPWA = zeros(1, length(SST));
MLD = cell(1, length(SST));

for im = 1:length(SST),
    
    sys_type = sub_sysStruct_type(SST{im});
    
    if Options.force_mld & ~sys_type.mld,
        error('Cannot force MLD model because it is not available.');
    end
    if Options.force_pwa & ~sys_type.pwa,
        error('Cannot force PWA model because it is not available.');
    end
    
    if Options.force_mld,
        dynamics_type{im} = 'mld';
        
    elseif Options.force_pwa,
        dynamics_type{im} = 'pwa';

    elseif sys_type.mld & Options.yalmip_online,
        % for on-line controllers we prefer MLD models
        dynamics_type{im} = 'mld';
        
    elseif sys_type.pwa,
        dynamics_type{im} = 'pwa';
        
    elseif sys_type.mld,
        dynamics_type{im} = 'mld';

    elseif sys_type.ml3,
        dynamics_type{im} = 'ml3';

    elseif sys_type.lti,
        dynamics_type{im} = 'lti';
        
    else
        error('Unknown type of dynamical model.');
    end
    if sys_type.nonlin,
        dynamics_type{im} = 'nonlin';
    end
    
    switch dynamics_type{im},
        case {'mld', 'ml3'},
            MLD{im} = SST{im}.data.MLD;
            if isequal(dynamics_type{im}, 'mld')
                % initialize B5 and D5 matrices to zero if they are not given
                if ~isfield(MLD{im}, 'B5'),
                    MLD{im}.B5 = zeros(size(MLD{im}.A, 1), 1);
                end
                if ~isfield(MLD{im}, 'D5'),
                    MLD{im}.D5 = zeros(size(MLD{im}.C, 1), 1);
                end
            end

        case 'pwa',
            nPWA(im) = length(SST{im}.A);
            
        case 'lti',
            % convert the LTI system into a PWA form, it will allow simplier code later
            SST{im} = mpt_lti2pwa(SST{im});
            sysStruct = mpt_lti2pwa(sysStruct);
            nPWA(im) = 1;
            
        case 'nonlin',
            % find out whether we have a piecewise non-linear system or just
            % a pure non-linear plant
            if iscell(SST{im}.A),
                % piecewise non-linear system
                nPWA(im) = length(SST{im}.A);
                dynamics_type{im} = 'pnonlin';
                
            else
                % pure non-linear system
                nPWA(im) = 1;                
                % convert it to a dummy PWA system with one dynamics
                SST{im} = mpt_lti2pwa(SST{im});
                sysStruct = mpt_lti2pwa(sysStruct);
                dynamics_type{im} = 'nonlin';
                
            end
            
    end
    
end


%-----------------------------------------------------------------
function [soften, slacks, smax, sweights] = sub_prepareslacks(probStruct, N, haveXbounds, haveYbounds, dims)
% introduce slack variables for soft constraints if necessary

global mptOptions

soften = struct('all', 0, 'x', 0, 'u', 0, 'y', 0);
smax = struct('all', 0, 'x', 0, 'u', 0, 'y', 0);
sweights = struct('all', 0, 'x', 0, 'u', 0, 'y', 0);
slacks.all = cell(1, N); [slacks.all{:}] = deal(0);
slacks.x = cell(1, N); [slacks.x{:}] = deal(zeros(dims.nx, 1));
slacks.y = cell(1, N); [slacks.y{:}] = deal(zeros(dims.ny, 1));
slacks.u = cell(1, N); [slacks.u{:}] = deal(zeros(dims.nu, 1));

if isfield(probStruct, 'S') | isfield(probStruct, 'smax'),
    % only one slack which softens state, input and output constraints
    % simultaneously
    slacks.all = sub_cellsdpvar(1, N);
    soften.all = 1;
    % upper bound on this slack
    smax.all = mpt_defaultField(probStruct, 'smax', mptOptions.infbox);
    % penalty on slacks
    sweights.all = mpt_defaultField(probStruct, 'S', 1e3);
    
else
    if isfield(probStruct, 'Sx') | isfield(probStruct, 'sxmax'),
        % softening of state constraints
        if any(haveXbounds==0),
            fprintf('WARNING: no state constraints given, cannot soften them.\n');
        else
            slacks.x = sub_cellsdpvar(dims.nx, N);
            soften.x = 1;
            % upper bound on this slack
            smax.x = mpt_defaultField(probStruct, 'sxmax', mptOptions.infbox);
            % penalty on slacks
            sweights.x = mpt_defaultField(probStruct, 'Sx', 1e3);
        end
    end
    if isfield(probStruct, 'Sy') | isfield(probStruct, 'symax'),
        % softening of output constraints
        if any(haveYbounds==0),
            fprintf('WARNING: no output constraints given, cannot soften them.\n');
        else
            slacks.y = sub_cellsdpvar(dims.ny, N-1);
            soften.y = 1;
            % upper bound on this slack
            smax.y = mpt_defaultField(probStruct, 'symax', mptOptions.infbox);
            % penalty on slacks
            sweights.y = mpt_defaultField(probStruct, 'Sy', 1e3);
        end
    end
    if isfield(probStruct, 'Su')  | isfield(probStruct, 'sumax'),
        % softening of input constraints
        slacks.u = sub_cellsdpvar(dims.nu, N-1);
        soften.u = 1;
        % upper bound on this slack
        smax.u = mpt_defaultField(probStruct, 'sumax', mptOptions.infbox);
        % penalty on slacks
        sweights.u = mpt_defaultField(probStruct, 'Su', 1e3);
    end
end

smax.all = smax.all(:); smax.x = smax.x(:); smax.y = smax.y(:); smax.u = smax.u(:);


%---------------------------------------------------------
function [minb, maxb] = sub_getMaxBounds(SST, flag, rows)
% from a cell array of sysStructs, extracts maximum bounds on a given variable
%
% "flag" can be either 'u', 'y', or 'x'
% "rows" denotes the dimension of a respective variable

global mptOptions

minb = repmat(Inf, rows, 1);
maxb = repmat(-Inf, rows, 1);
for ii = 1:length(SST),
    if isfield(SST{ii}, [flag 'max']),
        maxv = getfield(SST{ii}, [flag 'max']);
        maxb = max(maxv, maxb);
    end
    if isfield(SST{ii}, [flag 'min']),
        minv = getfield(SST{ii}, [flag 'min']);
        minb = min(minv, minb);
    end
end
% replace infinite bounds by +/-1e4
maxb(find(isinf(maxb))) = mptOptions.infbox;
minb(find(isinf(minb))) = -mptOptions.infbox;


%---------------------------------------------------------
function out = sub_cellsdpvar(dim, count)
% creates a cell array of sdpvars of dimension "dim"

out = cell(1, count);
for i = 1:count,
    out{i} = sdpvar(dim, 1);
end


%---------------------------------------------------------
function out = sub_check_nonlinearity(expression)
% checks whether a given non-linear function can be solved via YALMIP
% we currently reject sigmonial expressions, such as x/u

out = '';
if any(is(set(expression), 'sigmonial')),
    out = 'Only polynomial nonlinearities are allowed.';
end


%-------------------------------------------------------------------------
function [SST, pst, data] = sub_augmentx0(SST, probStruct, yalmip_online)
% augment system/problem to cope with tracking/deltaU constraints. 

% some parameters
do_tracking   = zeros(1, length(SST));
do_deltaU     = zeros(1, length(SST));
need_tracking = zeros(1, length(SST));
need_deltaU   = zeros(1, length(SST));
nref = 0;
[nx, nu, ny] = mpt_sysStructInfo(SST{1});
dims = struct('nx', nx, 'nu', nu, 'ny', ny);

% keep the first call verbose
verOpt.verbose = 1;

% tracking=1 cannot be used for MLD/nonlinear systems (because we cannot
% augmente the description for deltaU formulation). therefore in such case we
% switch to tracking=2
%
% notice that we don't have to care whether the controller to be computed is
% explicit or on-line, since we do require PWA/LTI models for explicit control
% anyhow (this is checked by mpt_control)
sys_type = cell(1, length(SST));
for i = 1:length(SST),
    sys_type{i} = sub_sysStruct_type(SST{i});
end

if probStruct.tracking==1,
    unsupported_type = '';
    for ii = 1:length(SST),
        if sys_type{ii}.mld & ~sys_type{ii}.pwa,
            unsupported_type = 'MLD';
            break
        elseif sys_type{ii}.nonlin,
            unsupported_type = 'nonlinear';
            break
        end
    end
    
    if ~isempty(unsupported_type);,
        fprintf('tracking=1 cannot be used for %s systems, switching to tracking=2\n', ...
            unsupported_type);
        probStruct.tracking = 2;
    end
end
if probStruct.tracking > 0,
    for ii = 1:length(SST),
        if sys_type{ii}.mld & sys_type{ii}.pwa,
            % we can only deal with tracking if we switch to PWA
            % formulation. that involves clearing the MLD description from
            % the system structure.
            SST{ii} = rmfield(SST{ii}, 'data');
        end
    end
end

pst = probStruct;

for ii = 1:length(SST),
    [do_tracking(ii), do_deltaU(ii), need_tracking(ii), need_deltaU(ii)] = ...
        sub_can_do_tracking(SST{ii}, probStruct);
    % do_tracking=1    - we are able to deal with tracking=1 by
    %                    directly augmenting the state vector
    % do_tracking=2    - we are able to deal with tracking=2 by
    %                    directly augmenting the state vector and/or by
    %                    introducing a new variable to denote the reference
    % do_deltaU=1      - we are able to deal with deltaU formulation by directly
    %                    augmenting the state vector
    % need_tracking=1  - tracking is implied by sysStruct/probStruct
    % need_deltaU=1    - deltaU formulation is implied by sysStruct/probStruct
    
    if do_tracking(ii),
        % augment the system to deal with tracking
        [SST{ii}, pst] = mpt_yalmipTracking(SST{ii}, probStruct, verOpt);
        nref = SST{ii}.dims.nref;

    elseif do_deltaU(ii),
        % augment the system to deal with deltaU constraints
        [SST{ii}, pst] = mpt_yalmipDU(SST{ii}, probStruct, verOpt);
        
    elseif probStruct.tracking > 0,
        % it can be that tracking is enabled but an MLD system is used. in such
        % case we must set sysStruct.dims manually (would otherwise be done by
        % mpt_yalmipTracking)
        SST{ii}.dims = dims;
        
    end
    
    % keep all other calls silent
    verOpt.verbose = -1;
end

% check that we either do tracking/deltaU formulation for all systems or for
% none of them
if (any(need_tracking) & ~all(need_tracking)) | ...
        (any(do_tracking) & ~all(do_tracking)),
    error('Tracking not enabled or possible for some specified systems.');
end
if (any(need_deltaU) & ~all(need_deltaU)) | ...
        (any(do_deltaU) & ~all(do_deltaU)), 
    error('deltaU constraints/penalties not enabled or possible for some specified systems.');
end
if any(do_tracking==1) & ~all(do_tracking==1),
    error('probStruct.tracking=1 not possible for some specified systems.');
end
if any(do_tracking==2) & ~all(do_tracking==2),
    error('probStruct.tracking=2 not possible for some specified systems.');
end
if any(need_tracking) & ~all(do_tracking) & yalmip_online==0,
    % we are able to handle tracking for on-line controllers by introducing a
    % variable for the reference. for explicit controllers we can't do it
    error('Tracking not possible for some specified systems.');
end
if any(need_deltaU) & ~(all(do_deltaU) | all(do_tracking==1)) & yalmip_online==0,
    % we are able to handle deltaU formulation for on-line controllers by
    % introducing a variable for the previous input. for explicit controllers we
    % can't do it 
    error('deltaU constraints/penalties not possible for this setup.');
end
if any(need_deltaU) & any(do_tracking==2) & yalmip_online==0,
    error('tracking=2 cannot be used with deltaU constraints/penalties. Use tracking=1 instead.');
end

data.uprev_in_x0 = (any(need_deltaU) & all(do_deltaU)) | any(do_tracking==1);
data.need_uprev = (any(need_deltaU) & ~(all(do_deltaU) | all(do_tracking==1))) | ...
    (any(need_tracking==1) & ~all(do_tracking==1));

data.reference_in_x0 = any(need_tracking) & all(do_tracking);
data.need_reference = any(need_tracking) & ~all(do_tracking);

data.nref = nref;


%-------------------------------------------------------------------------
function [do_tracking, do_deltaU, need_tracking, need_dU] = ...
    sub_can_do_tracking(sysStruct, probStruct)
% this function checks whether we can _augment_ the state vector to cope with
% tracking (mpt_yalmipTracking) or deltaU formulation (mpt_yalmipDU). if we
% can't we still can use on-line MPC and create new variables to denote the
% previous control action and the reference.

sys_type = sub_sysStruct_type(sysStruct);
[nx,nu,ny,ndyn,nbool] = mpt_sysStructInfo(sysStruct);

need_tracking = 0;
if ~isfield(probStruct, 'tracking_augmented'),
    need_tracking = probStruct.tracking;
end
need_dU = isfield(probStruct, 'Rdu') | ~(all(isinf(sysStruct.dumax)) & ...
    all(isinf(sysStruct.dumin)));

if need_tracking & probStruct.tracking==1 & need_dU,
    % deltaU formulation is a subset of tracking=1
    need_dU = 0;
end

if (sys_type.mld & ~sys_type.pwa) | sys_type.nonlin,
    % no augmentation possible for MLD and nonlinear systems
    cantr = 0; 
    candU = 0;

elseif nbool > 0 & need_dU,
    % cannot augment the state by u(k-1) if we have boolean inputs
    cantr = 0;
    candU = 0;
    
elseif nbool > 0 & probStruct.tracking,
    % we can only handle tracking=2, i.e. we can augment the state only by the
    % reference, not by u(k-1)
    cantr = 2;
    candU = 0;

else
    % we can deal with tracking=1, tracking=2 and deltaU formulation for LTI and
    % PWA systems
    cantr = probStruct.tracking; 
    candU = 1;
    
end

if need_tracking & cantr > 0,
    do_tracking = cantr;
    do_deltaU = 0;    % deltaU formulation is implied by tracking=1
    
elseif need_dU & candU,
    do_tracking = 0;
    do_deltaU = 1;

else 
    do_tracking = 0;
    do_deltaU = 0;
    
end


%-------------------------------------------------------------------------
function SST = sub_force_bounds(SST, flag, dim, nobounds, infbox)
% sets missing state/input/output constraints to +/- infbox

for i = 1:length(nobounds),
    S = SST{nobounds(i)};
    S = setfield(S, [flag 'max'], repmat(infbox, dim, 1));
    S = setfield(S, [flag 'min'], repmat(-infbox, dim, 1));
    SST{nobounds(i)} = S;
end
