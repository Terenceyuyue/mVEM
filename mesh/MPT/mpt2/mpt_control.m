function ctrl=mpt_control(sysStruct,probStruct,ctrltype,Options)
%MPT_CONTROL Main control routine. Computes explicit controller for a given problem
%
% ctrl=mpt_control(sysStruct,probStruct)
% ctrl=mpt_control(sysStruct,probStruct, 'on-line')
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Mother of all control functions. Main control routine of the MPT toolbox.
% Based on the problem definition it calls the appropriate multiparametric 
% controller function to compute explicit controller for both LTI
% and PWA systems. This should be the only function the user needs to call. 
%
% ---------------------------------------------------------------------------
% USAGE
% ---------------------------------------------------------------------------
%
% To compute an explicit controller, call:
%   ctrl = mpt_control(sysStruct, probStruct)
%
% To compute an on-line controller, call:
%   ctrl = mpt_control(sysStruct, probStruct, 'on-line')
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct    - System structure in the sysStruct format
%                (see 'help mpt_sysStruct' for more details)
% probStruct   - Problem structure in the probStruct format
%                (see 'help mpt_probStruct' for more details)
% Options      - Optional: User options to be passed to individual control functions
% Options.autoTracking
%           If set to 0, system and problem matrices are assumed to be
%           augmented for tracking purpose by the user. By default,
%           matrices will be extended automatically to guarantee tracking.
% Options.noExtreme
%           If set to 0 (default), extreme points of each polytope in the
%           controller partition will be computed for faster plotting. Set
%           it to 1 if you don't want to compute the extreme points.
%
%	See the MPT Manual for additional details on the structure format or 
%   consult one of the example systems (e.g. Double_Integator) which were
%	provided with this package.                          
%                        
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrl   - Controller object (see 'help mptctrl' for more details)
%
% see also MPT_OPTCONTROL, MPT_OPTINFCONTROL, MPT_ITERATIVE, MPT_ITERATIVEPWA,
%          MPT_OPTCONTOLPWA, MPT_YALMIPCFTOC
%

% Copyright is with the following author(s):
%
%(C) 2003-2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%              kvasnica@control.ee.ethz.ch

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

error(nargchk(2,4,nargin));

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

if ~(isstruct(sysStruct) | iscell(sysStruct))
    error('MPT_CONTROL: First input argument must be a system structure!');
end
if ~isstruct(probStruct)
    error('MPT_CONTROL: First input argument must be a problem structure!');
end

if nargin==2,
    ctrltype = 'explicit';
    Options = [];
elseif nargin==3 & ischar(ctrltype)
    Options = [];
elseif nargin==3 & (isstruct(ctrltype) | isempty(ctrltype))
    Options = ctrltype;
    ctrltype = 'explicit';
elseif nargin==3 & ~ischar(ctrltype)
    error('MPT_CONTROL: Third input argument must be a string! Allowed values are ''explicit'' and ''on-line''.');
elseif nargin==4 & ~ischar(ctrltype)
    error('MPT_CONTROL: Third input argument must be a string! Allowed values are ''explicit'' and ''on-line''.');
elseif nargin==4 & ~isstruct(Options)
    error('MPT_CONTROL: Fourth input argument must be an Options structure!');
end

% ============================================================================
% set default options
Options = mpt_defaultOptions(Options, ...
    'noExtreme', 0, ...
    'autoTracking', 1, ...
    'statusbar', 0, ...
    'qpsolver', mptOptions.qpsolver, ...
    'sysstructname', inputname(1), ...
    'probstructname', inputname(2) );

if isempty(Options.sysstructname),
    Options.sysstructname = 'sysStruct';
end
if isempty(Options.probstructname),
    Options.probstructname = 'probStruct';
end


userSysStruct = sysStruct;
if iscell(sysStruct),
    sysStruct = sysStruct{1};
end

% ============================================================================
% verify system/problem structures
origSysStruct_notverified = sysStruct;
origProbStruct_notverified = probStruct;
if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt = Options;
    verOpt.verbose=1;
    verOpt.useyalmip = 1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end


% ============================================================================
% there are certain system/problem formulations which mpt_yalmipcftoc() cannot
% handle
noyalmip_because = [];
havedymax = 0; havedymin = 0;
if isfield(sysStruct, 'dymax'),
    havedymax = any(~isinf(sysStruct.dymax));
end
if isfield(sysStruct, 'dymin'),
    havedymin = any(~isinf(sysStruct.dymin));
end
if isinf(probStruct.N),
    noyalmip_because = 'prediction horizon is infinity';
elseif probStruct.subopt_lev > 0,
    noyalmip_because = sprintf('%s.subopt_lev is not 0', Options.probstructname);
elseif probStruct.feedback==1,
    noyalmip_because = 'feedback pre-stabilization is enabled';
elseif mpt_isnoise(sysStruct.noise),
    noyalmip_because = 'system with additive noise';
elseif isfield(sysStruct, 'Aunc') | isfield(sysStruct, 'Bunc'),
    noyalmip_because = 'system with parametric uncertainties';
elseif havedymax | havedymin,
    noyalmip_because = 'dymax/dymin constraints';
end


% ============================================================================
% there are certain system/problem formulations which MPTs native functions
% cannot handle
nompt_because = [];
if isfield(probStruct, 'Sx') | isfield(probStruct, 'Sy') | isfield(probStruct, 'Su') | ...
        isfield(probStruct, 'sxmax') | isfield(probStruct, 'symax') | ...
        isfield(probStruct, 'sumax') | isfield(probStruct, 'S')
    nompt_because = 'soft constraints';
elseif isfield(probStruct, 'xN'),
    nompt_because = 'terminal state constraint';
elseif isfield(probStruct, 'Qdyn') | isfield(probStruct, 'Qswitch'),
    nompt_because = 'penalized switching';
elseif iscell(userSysStruct),
    nompt_because = 'time-varying model';
elseif isfield(sysStruct, 'nonlinhandle'),
    nompt_because = 'non-linear dynamics';
elseif iscell(sysStruct.A) & isfield(probStruct, 'Nc') & ...
        probStruct.subopt_lev~=1,
    nompt_because = 'control horizon for PWA systems';
elseif probStruct.norm~=2 & isfield(probStruct, 'Nc'),
    nompt_because = 'control horizon with linear cost';
end


% ============================================================================
% give an error if neither mpt_yalmipcftoc() nor MPTs native functions can
% handle a given system/problem setup
error(sub_cannotsolve(noyalmip_because, nompt_because));


% ============================================================================
% deal with on-line controllers
if strcmpi(ctrltype, 'on-line') | strcmpi(ctrltype, 'online')
    % compute matrices of on-line controller and return
    Options = mpt_defaultOptions(Options, ...
        'force_mpt', ~isempty(noyalmip_because) );
    if Options.force_mpt & ~isempty(nompt_because),
        fprintf('Cannot force MPT because of following problem:\n');
        fprintf(' * %s\n\n', nompt_because);
        error('Cannot continue.');
    end
    if iscell(userSysStruct),
        % handle multi-model dynamics
        ctrl = mptctrl(userSysStruct, probStruct, Options);
    else
        ctrl = mptctrl(sysStruct, probStruct, Options);
    end
    return
end


% ============================================================================
% cannot compute an explicit controller for MLD systems with no corresponding
% PWA representation
% Note: this is not fully true, since mpt_yalmipcftoc() can do that, but we
% reject it here because many other functions rely on existence of the PWA
% representation. but if you nevertheless want to compute an explicit controller
% for MLD systems, use this call:
%
%   ctrl = mpt_yalmipcftoc(sysStruct, probStruct)
%
% and it will return a valid controller object.

if isfield(sysStruct, 'data'),
    if isfield(sysStruct.data, 'onlymld'),
        if sysStruct.data.onlymld,
            % cannot compute an explicit controller if PWA model is not
            % available
            fprintf('\nPWA representation of the hybrid system must be available in order to compute explicit solutions.\n');
            fprintf('Call "%s = mpt_sys(%s.data.MLD)" to get the PWA representation.\n\n', Options.sysstructname, Options.sysstructname);
            error('Cannot compute explicit solution if PWA representation is not available in sysStruct.');
        end
    end
end


% ============================================================================
% no explicit solution for non-linear system
if isfield(sysStruct, 'nonlinhandle'),
    error('No explicit solution available for non-linear systems.');
end


origSysStruct = sysStruct; origProbStruct = probStruct;
if probStruct.tracking==2 & (isfield(probStruct, 'Rdu') | ...
        (any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax))))
    % if we have deltaU constraints and tracking without deltaU formulation is
    % enabled, we use tracking with deltaU formulation which satisfies both
    % needs
    probStruct.tracking = 1;
end
if length(sysStruct.Pbnd)>1 & ~iscell(sysStruct.A),
    % Pbnd is (possibly) non-convex and the system is LTI, convert it to PWA
    % format with sysStruct.guardX and sysStruct.guardC corresponding to each
    % element of sysStruct.Pbnd
    sysStruct = sub_expandPbnd(sysStruct);
end
[nx,nu,ny,ndyn,nbool,ubool] = mpt_sysStructInfo(sysStruct);
isPWA = iscell(sysStruct.A);

mpt_call = []; yalmip_call = []; no_overlaps = 0; check_feasibility = 0;
prefered = 'yalmip';

switch probStruct.subopt_lev,
    case 0,
        % optimal solution
        
        if isinf(probStruct.N),
            % CITOC problem (infinite-time solution)
            
            if nbool > 0,
                error('Systems with discrete inputs not supported for infinite-time solutions.');
            end
            
            if isPWA & probStruct.norm == 2,
                error('No infinite-time solution for PWA systems and quadratic cost.'); 
                
            elseif ~isPWA & probStruct.norm == 2,
                % CITOC for quadratic cost
                mpt_call = @mpt_optInfControl;
                
            elseif probStruct.norm ~= 2,
                % use CITOC for PWA systems and linear cost, this function also
                % handles LTI systems 
                mpt_call = @mpt_optInfControlPWA;
                
            else
                error('Unsupported model/problem combination.');
                
            end
            
        else
            % CFTOC problem

            if nbool > 0,
                % discrete inputs, mpt_yalmipcftoc() handles all norms and
                % PWA/LTI systems
                yalmip_call = @mpt_yalmipcftoc;
            
            elseif isPWA,
                % PWA systems with no discrete inputs
                
                if probStruct.norm == 2,
                    % mpt_yalmipcftoc() can do the DP for quadratic cost
                    yalmip_call = @mpt_yalmipcftoc;
                    
                else
                    % mpt_optControlPWA() is our fast DP procedure for PWA
                    % systems with linear cost, but we can also use
                    % mpt_yalmipcftoc() which also does DP (less efficiently,
                    % though)
                    mpt_call = @mpt_optControlPWA;
                    yalmip_call = @mpt_yalmipcftoc;
                    
                    % mpt_optControlPWA() is faster than mpt_yalmipcftoc()
                    prefered = 'mpt';
                    
                end
                
            elseif ~isPWA,
                % linear systems with no discrete inputs
                yalmip_call = @mpt_yalmipcftoc;
                
                % for linear costs we would like to use our fast DP procedure in
                % mpt_optControlPWA, but that one requires the latest mpLP
                % solver to be working. And we can't use the latest mpLP if
                % there is no QP solver available.
                %
                % Moreover, we must not use mpt_optControlPWA for systems with
                % uncertainties (additive or parametric)
                if isfield(Options, 'mplpver'),
                    mplpver = Options.mplpver;
                else
                    % choose the fastest
                    mplpver = Inf;
                end
                oneshot = 0;
                % solve the problem in one shot if:
                %  - older version of the MPLP solver is requested, or
                %  - no QP solver is available (in which case we cannot
                %    use newer version of MPLP solver)
                %  - system has uncertainty
                if isfield(sysStruct, 'Aunc') | isfield(sysStruct, 'Bunc'),
                    oneshot = 1;
                elseif mplpver < 4 | Options.qpsolver==-1,
                    oneshot = 1;
                elseif isfield(sysStruct, 'noise'),
                    if mpt_isnoise(sysStruct.noise),
                        oneshot = 1;
                    end
                end                 
                
                if probStruct.norm == 2,
                    % CFTOC for linear systems
                    mpt_call = @mpt_optControl;
                    no_overlaps = 1;
                    
                    % mpt_yalmipcftoc() is more numerically robust for this case
                    prefered = 'yalmip';
                    
                elseif oneshot,
                    % solve the one-shot formulation
                    mpt_call = @mpt_optControl;
                    no_overlaps = 1;
                    
                else
                    % use DP approach with PWA cost-to-go (faster than
                    % mpt_optControl for linear objectives)
                    mpt_call = @mpt_optControlPWA;
                    
                    % mpt_optControlPWA() is faster than mpt_yalmipcftoc()
                    prefered = 'mpt';
                
                end
                
            else
                error('Unsupported model/problem combination.');
                
            end
        end
        
    case 1,
        % minimum-time solution
        
        if probStruct.tracking,
            error('Tracking not supported for minimum-time problems.');
        end

        if nbool == nu,
            % all inputs are boolean
            mpt_call = @mpt_boolMinTime;
            
        elseif nbool > 0,
            % some inputs boolean, some continuous
            mpt_call = @mpt_mixedMinTime;
            
        elseif isPWA,
            % time-optimal solution for PWA systems
            mpt_call = @mpt_iterativePWA;
            
        elseif ~isPWA,
            % time-optimal solution for LTI systems
            mpt_call = @mpt_iterative;
            
        else
            error('Unsupported model/problem combination.');
            
        end
        
    case 2,
        % low-complexity solution
        
        if nbool > 0,
            error('probStruct.subopt_lev=2 not supported for systems with discrete inputs!');
        end
        
        if isPWA,
            % for PWA systems we can only obtain the low-complexity solution for
            % linear cost
            
            if probStruct.norm==2
                error('2-norm problems not allowed for subopt_lev=2 and PWA systems!');
                
            else
                mpt_call = @mpt_iterativePWA;
                
            end
            
        elseif ~isPWA,
            mpt_call = @mpt_oneStepCtrl;
            check_feasibility = 1;
            
        else
            error('Unsupported model/problem combination.');
            
        end
        
    otherwise,
        error('Unknown type of probStruct.subopt_lev');
        
end


% ============================================================================
% perform computation
method_used = '';
canuse_yalmip = 0; canuse_mpt = 0;

if ~isempty(yalmip_call) & isempty(noyalmip_because),
    % we can use mpt_yalmipcftoc() for this model/problem, let's do so
    function_call = yalmip_call;
    method_used = 'yalmip';
    canuse_yalmip = 1;
end
if ~isempty(mpt_call) & isempty(nompt_because),
    % we can use MPTs native function for this setup
    
    function_call = mpt_call;
    method_used = 'mpt';
    canuse_mpt = 1;
end
if canuse_mpt==0 & canuse_yalmip==0,
    % no suitable method to deal with given setup
    if isempty(noyalmip_because),
        noyalmip_because = 'No appropriate function to call with this setup.';
    end
    if isempty(nompt_because),
        nompt_because = 'No appropriate function to call with this setup.';
    end
    error(sub_cannotsolve(noyalmip_because, nompt_because));
end

if isfield(Options, 'prefered'),
    prefered = Options.prefered;
end

if canuse_yalmip & canuse_mpt,
    % we can use both approaches, decide now based on our preference
    % (mpt_yalmipcftoc() is prefered always but for CFTOC for PWA systems, where
    % mpt_optControlPWA() is faster)
    if isequal(prefered, 'mpt'),
        % prefer MPTs native functions
        method_used = 'mpt';
        function_call = mpt_call;
        
    elseif isequal(prefered, 'yalmip'),
        % prefer mpt_yalmipcftoc()
        method_used = 'yalmip';
        function_call = yalmip_call;
    end
end

if isequal(method_used, 'mpt')
    % previously we have verified user input with Options.useyalmip=1 which may
    % disable some important checks which have to be performed if MPT-based
    % function are to be used. therefore we verify the structures once again
    % with Options.useyalmip disabled.
    verOpt = Options;
    verOpt.verbose = -1;
    verOpt.useyalmip = 0;
    [sysStruct,probStruct]=mpt_verifySysProb(origSysStruct_notverified, ...
        origProbStruct_notverified, verOpt);
    
    % we need to deal with tracking or feedback pre-stabilization if we call
    % MPTs native functions
    % (note that tracking is handled directly in mpt_yalmipcftoc())
    [sysStruct, probStruct] = sub_augmentSetup(sysStruct, probStruct, Options);

end

% now call the computation
if check_feasibility,
    [ctrlStruct, feasible] = feval(function_call, sysStruct, probStruct, Options);
    if ~feasible,
        disp('No Lyapunov function was found for this system, the final result may not be a stabilizing controller!!!');
    else
        disp('Lyapunov function found, controller is stabilizing');
    end
    
else
    if iscell(userSysStruct),
        % handle multi-model dynamics (sysStruct is a cell)
        ctrlStruct = feval(function_call, userSysStruct, probStruct, Options);
    else
        ctrlStruct = feval(function_call, sysStruct, probStruct, Options);
    end
end


% ============================================================================
% check for possible errors

if ~exist('ctrlStruct', 'var'),
    error('mpt_control: Problem is infeasible...');
end
if isempty(ctrlStruct.Pn),
    error('mpt_control: Problem is infeasible...');
end
if ~isfulldim(ctrlStruct.Pn) | length(ctrlStruct.Fi)==0,
    error('mpt_control: Problem is infeasible...');
end
if isempty(ctrlStruct.Fi{1}),
    error('mpt_control: Problem is infeasible...');
end

nR = length(ctrlStruct.Fi);
if isa(ctrlStruct, 'mptctrl'),
    ctrlStruct = struct(ctrlStruct);
end

if isequal(method_used, 'mpt') & (isfield(probStruct,'xref') | isfield(probStruct,'uref')),
    % for fixed-state tracking, do the substitution back (but only for MPTs
    % native functions
    Pn = [];
    xref = probStruct.xref;
    uref = probStruct.uref;
    for reg = 1:nR,
        nu = round(size(ctrlStruct.Fi{reg}, 1) / length(uref(:)));
        ctrlStruct.Gi{reg} = ctrlStruct.Gi{reg} - ...
            ctrlStruct.Fi{reg}*xref + repmat(uref,nu,1);
        if probStruct.subopt_lev==0,
            % for cost-optimal solution, translate also the cost
            ctrlStruct.Ci{reg} = ctrlStruct.Ci{reg} - ctrlStruct.Bi{reg}*xref + xref'*ctrlStruct.Ai{reg}*xref;
            ctrlStruct.Bi{reg} = ctrlStruct.Bi{reg} - 2*xref'*ctrlStruct.Ai{reg};
            % quadratic term remains identical
        end
        % translate regions:
        [H,K] = double(ctrlStruct.Pn(reg));
        Pn = [Pn polytope(H, K + H*xref)];
    end
    
    % translate feasible set:
    Pfinal = polytope;
    for fin = 1:length(ctrlStruct.Pfinal),
        if ~isfulldim(ctrlStruct.Pfinal(fin)), continue, end
        [Hf,Kf] = double(ctrlStruct.Pfinal(fin));
        Pfinal = [Pfinal polytope(Hf, Kf + Hf*xref)];
    end
    if ~isfulldim(Pfinal),
        Pfinal = Pn;
    end
    
    % store the new data
    ctrlStruct.Pfinal = Pfinal;
    ctrlStruct.Pn = Pn;
    ctrlStruct.details.origSysStruct = origSysStruct;
    ctrlStruct.details.origProbStruct = origProbStruct;
    ctrlStruct.sysStruct = sysStruct;    
    ctrlStruct.probStruct = probStruct;
end


% ============================================================================
% compute extreme points of controller regions if asked for
if dimension(ctrlStruct.Pn)<=3 & length(ctrlStruct.Pn)<1000 & Options.noExtreme==0,
    % if dimension is <= 3, we compute and store also extreme points for faster plotting
    Pn = ctrlStruct.Pn;
    for ii=1:length(Pn),
        [V,R,Pn(ii)]=extreme(Pn(ii));
    end
    ctrlStruct.Pn = Pn;
end


% ============================================================================
% remove expression of the cost function from details (mpt_optControl() returns
% it there)
if isfield(ctrlStruct.details,'Bi'),
    ctrlStruct.details = rmfield(ctrlStruct.details,'Bi');
end
if isfield(ctrlStruct.details,'Ci'),
    ctrlStruct.details = rmfield(ctrlStruct.details,'Ci');
end
if isfield(ctrlStruct.details,'Pn'),
    ctrlStruct.details = rmfield(ctrlStruct.details,'Pn');
end

ctrlStruct.details.origSysStruct = origSysStruct;
ctrlStruct.details.origProbStruct = origProbStruct;

fprintf('\nSolution consists of %d regions\n\n',nR);

try
    ctrl = mptctrl(ctrlStruct);
catch
    fprintf('An unexpected error occured when creating an MPT controller object\n');
    fprintf('Please be so kind and send your system and problem definition to:\n');
    fprintf('   mpt@control.ee.ethz.ch\n\n');
    fprintf('We will investigate the problem.\n');
    fprintf('Meanwhile, the controller is returned as a ctrlStruct structure...\n\n');
    ctrl = ctrlStruct;
end




%--------------------------------------------------------------------------
function sysStruct = sub_expandPbnd(sysStruct)

if length(sysStruct.Pbnd)>1 & ~iscell(sysStruct.A),
    % Pbnd is (possibly) non-convex and the system is LTI, convert it to PWA
    % format with sysStruct.guardX and sysStruct.guardC corresponding to each
    % element of sysStruct.Pbnd
    warning('sysStruct.Pbnd is a polyarray, converting LTI system to PWA form...');
    sysStruct = mpt_lti2pwa(sysStruct);
    npbnd = length(sysStruct.Pbnd);
    for ii = 2:npbnd,
        sysStruct.A{ii} = sysStruct.A{1};
        sysStruct.B{ii} = sysStruct.B{1};
        sysStruct.C{ii} = sysStruct.C{1};
        sysStruct.D{ii} = sysStruct.D{1};
        sysStruct.f{ii} = sysStruct.f{1};
        sysStruct.g{ii} = sysStruct.g{1};
        sysStruct.guardX{ii} = sysStruct.guardX{1};
        sysStruct.guardU{ii} = sysStruct.guardU{1};
        sysStruct.guardC{ii} = sysStruct.guardC{1};
    end
    for ii = 1:npbnd,
        [sysStruct.guardX{ii},sysStruct.guardC{ii}] = double(sysStruct.Pbnd(ii));
    end
    sysStruct.Pbnd = hull(sysStruct.Pbnd);
end


%--------------------------------------------------------------------------
function [sysStruct, probStruct] = sub_augmentSetup(sysStruct, probStruct, Options)
% augment setup to deal with tracking or feedback pre-stabilization

if (probStruct.tracking & Options.autoTracking) | isfield(probStruct,'xref') | ...
        isfield(probStruct,'uref'),
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

if probStruct.feedback 
    if iscell(sysStruct.A),
        error('Feedback Pre-Stabilization not supported for PWA systems!');
    end
    if ~isfield(probStruct,'FBgain')
        % compute the pre-stabilization feedback for LTI system if it is not given
        [FB,S,E] = mpt_dlqr(sysStruct.A,sysStruct.B,probStruct.Q,probStruct.R);
        probStruct.FBgain = -FB;
    end
end


%--------------------------------------------------------------------------
function err=sub_cannotsolve(noyalmip_because, nompt_because)
% returns a non-empty error messages if we can't solve a given problem neither
% by mpt_yalmipcftoc() nor by MPTs native functions

err = '';
if ~isempty(noyalmip_because) & ~isempty(nompt_because),
    fprintf('\nCannot solve this setup because of conflicting objectives/constraints:\n');
    if ~isempty(noyalmip_because)
        fprintf(' * %s\n', noyalmip_because);
    end
    if ~isempty(nompt_because),
        fprintf(' * %s\n', nompt_because);
    end
    fprintf('\nPlease modify your system/problem setup to remove one of the conflicts.\n\n');
    err='Cannot handle given setup, see message above.';
end    
