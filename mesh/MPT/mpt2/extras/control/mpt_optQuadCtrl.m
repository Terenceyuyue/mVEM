function ctrl = mpt_optQuadCtrl(sysStruct, probStruct, Options)
% MPT_OPTQUADCTRL CFTOC of PWA systems with quadratic cost
%
% ctrl = mpt_optQuadCtrl(sysStruct, probStruct)
% ctrl = mpt_optQuadCtrl(sysStruct, probStruct, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the explicit solution of a given CFTOC problem with quadratic cost
% fucntion.
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
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2005 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
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

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3,
    Options = [];
end
Options = mpt_defaultOptions(Options, ...
    'verbose', mptOptions.verbose);

%===============================================================================
% set options which we pass to solvemp()
yalmipOptions = mptOptions.sdpsettings;
f = fields(Options);
for ii = 1:length(f),
    yalmipOptions = setfield(yalmipOptions, f{ii}, getfield(Options, f{ii}));
end


%===============================================================================
% verify sysStruct and probStruct structures
if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end
origSysStruct = sysStruct;
origProbStruct = probStruct;


%===============================================================================
% prepare data such that we can later create a valid MPTCTRL object
origSysStruct = sysStruct;
origProbStruct = probStruct;
if ~iscell(sysStruct.A),
    % convert LTI system to PWA format
    sysStruct = mpt_lti2pwa(sysStruct);
end


%===============================================================================
% create a cell array of probStruct.[R,Q,Qy] if they are not already provided as
% a cell array. this will be used later for simpler implementation of
% time-varying penalties
probStruct = sub_timevarpenalties(probStruct);


%===============================================================================
% check system dimensions
[nx,nu,ny,nPWA,nbool,ubool] = mpt_sysStructInfo(sysStruct);
% if nbool==nu,
%     % we can't handle cases where all inputs are fixed
%     error('mpt_optQuadCtrl: at least one system input must be continuous.');l
% end


%===============================================================================
% reject certain problem descriptions
if probStruct.norm~=2,
    error('mpt_optQuadCtrl: This function only supports quadratic cost functions.');
end
if isinf(probStruct.N),
    error('mpt_optQuadCtrl: Prediction horizon must be finite!');
end
if mpt_isnoise(sysStruct.noise),
    error('mpt_optQuadCtrl: additive noise not supported.');
end
if nbool>0 & probStruct.tracking==1,
    error('mpt_optQuadCtrl: tracking cannot be used for systems with integer inputs.');
end


%===============================================================================
% find out whether sysStruct.C, sysStruct.D and sysStruct.g all contain the same
% elements for all dynamics. if so, we can simplify things a bit
YeqSame = sub_allYeqsEqual(sysStruct);


%===============================================================================
% check constraints
haveXbounds = isfield(sysStruct, 'xmax');
if any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax)),
    dUconstraints = 1;
else
    dUconstraints = 0;
end
if dUconstraints & nbool>0,
    error('mpt_optQuadCtrl: deltaU constraints not supported for systems with discrete inputs.');
end


%===============================================================================
% when deltaU formulation is introduced, input constraints are set to +/- Inf.
% make them tighther.
sysStruct.umax(isinf(sysStruct.umax)) = 1e6;
sysStruct.umin(isinf(sysStruct.umin)) = -1e6;


starttime = cputime;

N = probStruct.N+1;

% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

if ~YeqSame,
    % Outputs y(k), ..., y(k+N) (last one not used)
    y = sdpvar(repmat(ny,1,N),repmat(1,1,N));
end

obj = 0;
F = set([]);


%===============================================================================
% extract dynamics and possible modify matrices to proper dimensions if we have
% tracking. I know this is ugly, but we have to do that to be consistent with
% mpt_constructMatrices. Basically if we have tracking for PWA systems, we
% create sysStruct.Cy and sysStruct.Dy such that when we do
%   sysStruct.Cy*x
% we obtain the difference (y - ref)
sA = sysStruct.A; sB = sysStruct.B; sF = sysStruct.f;
if isfield(sysStruct, 'Cy'),
    sC = sysStruct.Cy; sD = sysStruct.Dy; 
    ny = size(sC{1}, 1);
    sG = sysStruct.g;
    for ii = 1:nPWA,
        sG{ii} = sG{ii}(1:ny);
    end
    if isfield(probStruct, 'Qy'),
        for ii = 1:length(probStruct.Qy),
            probStruct.Qy{ii} = probStruct.Qy{ii}(1:ny, 1:ny);
        end
    end
else
    sC = sysStruct.C; sD = sysStruct.D; sG = sysStruct.g;
end


%===============================================================================
% prepare reference
if isfield(probStruct, 'yref'),
    yref = probStruct.yref;
else
    yref = zeros(ny, 1);
end

% here we could take xref and uref from probStruct.xref and probStruct.uref,
% respectively. but ,[t_control already takes over this and performes change of
% coordinates by calling mpt_prepareTracking, therefore we must use zero
% reference in the sequel 
xref = zeros(nx, 1);
uref = zeros(nu, 1);
% if isfield(probStruct, 'xref'),
%     xref = probStruct.xref;
% end
% if isfield(probStruct, 'uref'),
%     uref = probStruct.uref;
% end


%===============================================================================
% now formulate the CFTOC problem
for k = N-1:-1:1
    % Feasible region
     
    % Binary for PWA selection
    d = binvar(nPWA, 1);
    
    % set bounds on states and inputs
    if haveXbounds,
        bounds(x{k}, sysStruct.xmin, sysStruct.xmax);
        bounds(x{k+1}, sysStruct.xmin, sysStruct.xmax);
    end
    bounds(u{k}, sysStruct.umin, sysStruct.umax);
    if ~(k==1 & probStruct.y0bounds==0)
        % do not impose constraints on y0 if user does not want to
        if ~YeqSame,
            bounds(y{k}, sysStruct.ymin, sysStruct.ymax);
        end
    end
    
    % input constraints
    F = F + set(sysStruct.umin < u{k}     < sysStruct.umax);

    % some inputs can be boolean or from finite alphabet
    for iu = 1:nu,
        [t, s] = sub_inputtype(sysStruct, iu);
        if t=='B',
            % this input is boolean
            F = F + set(binary(u{k}(iu)));
        elseif t=='A',
            % this input is from finite alphabet
            F = F + set(ismember(u{k}(iu), s));
        end
    end

    % state constraints
    if haveXbounds,
        F = F + set(sysStruct.xmin < x{k}     < sysStruct.xmax);
        F = F + set(sysStruct.xmin < x{k+1}   < sysStruct.xmax);
    end
    
    % output constraints
    if ~(k==1 & probStruct.y0bounds==0)
        % do not impose constraints on y0 if user does not want to
        if YeqSame,
            % C,D,g are identical, it's enough to consider one dynamics
            F = F + set(sysStruct.ymin < sysStruct.C{1}*x{k} + sysStruct.D{1}*u{k} + sysStruct.g{1}   < sysStruct.ymax);
            F = F + set(sysStruct.ymin < sysStruct.C{1}*x{k+1} + sysStruct.D{1}*u{k+1} + sysStruct.g{1} < sysStruct.ymax);
        else
            F = F + set(sysStruct.ymin < y{k}   < sysStruct.ymax);
        end
    end

    % PWA Dynamics
    for i = 1:nPWA
        F = F + set(implies(d(i), x{k+1} == sA{i}*x{k} + sB{i}*u{k} + sF{i}));
        if ~YeqSame,
            % we need to have the output as a variable
            F = F + set(implies(d(i), y{k} == sC{i}*x{k} + sD{i}*u{k} + sG{i}));
        end
        F = F + set(implies(d(i),sysStruct.guardX{i}*x{k} + sysStruct.guardU{i}*u{k} <= sysStruct.guardC{i}));
    end
    F = F + set(sum(d) == 1);

    % add target set constraint
    if k==N-1 & isfulldim(probStruct.Tset),
        % ismember should automatically handle cases where Tset is a polytope
        % array
        F = F + set(ismember(x{k+1}, probStruct.Tset));
    end
    
    % objective function
    if isfield(probStruct, 'Qy'),
        % add penalty on outputs
        
        if YeqSame,
            % it's enough to consider just one dynamics because C,D,g are
            % identical for all dynamics
            obj = obj + (sC{1}*x{k} + sD{1}*u{k} + sG{1} - yref)' * ...
                probStruct.Qy{k} * (sC{1}*x{k} + sD{1}*u{k} + sG{1} - yref);
            
        else
            obj = obj + (y{k} - yref)' * probStruct.Qy{k} * (y{k} - yref);
            
        end

    else
        % add penalty on states
        obj = obj + (x{k} - xref)' * probStruct.Q{k} * (x{k} - xref);
        
        % add terminal weight if specified
        if k==N-1 & isfield(probStruct, 'P_N'),
            obj = obj + (x{k+1} - xref)' * probStruct.P_N * (x{k+1} - xref);
        end
        
    end
    
    % add penalty on inputs
    obj = obj + (u{k} - uref)' * probStruct.R{k} * (u{k} - uref);
    
end

%===============================================================================
% solve the mpMIQP
[mpsol{k},sol{k},Uz{k}] = solvemp(F, obj, yalmipOptions, x{k}, u{k});


%===============================================================================
% collect overlapping partitions together
if length(mpsol{k})==1,
    ctrl = mpsol{k}{1};
else
    ctrl = mpt_mergeCS(mpsol{k});
end

%===============================================================================
if isempty(ctrl),
    % problem either infeasible or some problem occured
    ctrl = mptctrl;
    return
end

% add necessary fields and create an MPTCTRL object
ctrl.overlaps = 1;
ctrl.sysStruct = origSysStruct;
ctrl.probStruct = origProbStruct;
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
    % input is continuous
    t = 'C';
    s = [];
end


%-----------------------------------------------------------------
function yesno = sub_allYeqsEqual(sysStruct)
% returns true if sysStruct.C, sysStruct.D and sysStruct.g are identical for all
% dynamics

yesno = 1;
C = sysStruct.C{1};
D = sysStruct.D{1};
g = sysStruct.g{1};
for jj = 2:length(sysStruct.C),
    if ~isequal(C, sysStruct.C{jj}) | ...
            ~isequal(D, sysStruct.D{jj}) | ~isequal(g, sysStruct.g{jj}),
        yesno = 0;
        return
    end
end

%-----------------------------------------------------------------
function probStruct = sub_timevarpenalties(probStruct)
% if probStruct.[R,Q,Qy] are single matrices, we convert them to a cell array of
% length probStruct.N

if ~iscell(probStruct.R),
    R = cell(1, probStruct.N);
    for ii = 1:probStruct.N,
        R{ii} = probStruct.R;
    end
    probStruct.R = R;
end

if ~iscell(probStruct.Q),
    Q = cell(1, probStruct.N);
    for ii = 1:probStruct.N,
        Q{ii} = probStruct.Q;
    end
    probStruct.Q = Q;
end

if isfield(probStruct, 'Qy'),
    if ~iscell(probStruct.Qy),
        Qy = cell(1, probStruct.N);
        for ii = 1:probStruct.N,
            Qy{ii} = probStruct.Qy;
        end
        probStruct.Qy = Qy;
    end
end
