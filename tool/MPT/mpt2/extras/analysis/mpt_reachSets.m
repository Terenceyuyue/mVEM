function [Psets, Vsets, PreachN, VreachN] = mpt_reachSets(structure, X0, U0, N, Options)
%MPT_REACHSETS Computes sets of reachable states for a given system / controller
%
% [Rsets, Vsets] = mpt_reachSets(structure, X0, U0, N, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates set of states which are reachable from a given polytope X0 in N
% steps. The system under consideration can be either linear (LTI), or PWA.
% Input argument can be either a system in sysStruct format, or an explicit
% controller computed by mpt_control. In that case, the autonomous system is
% defined as:
%   x(k+1) = (A+Fi)*x + (Gi+f)
%
% If the system is not autonomous, it is converted to such by setting
%   x(k+1) = A*x + fcl
% where
%   fcl = f + B*Options.U
%
% Options.U is set to zero by default.
%
% WARNING: By default, this function computes only fully dimensional reachable
% sets! Use "Options.lowdim=1" to tell the function to compute lower dimensional
% sets as well.
%
% USAGE:
%   Psets = mpt_reachSets(sysStruct, X0, N, Options)     - assumes u = 0
%   Psets = mpt_reachSets(sysStruct, X0, U0, N, Options) - assumes u \in U0
%   Psets = mpt_reachSets(ctrl, X0, N, Options)
%   [Psets, Vsets] = mpt_reachSets(sysStruct, X0, N, Options)
%   [Psets, Vsets] = mpt_reachSets(sysStruct, X0, U0, N, Options)
%   [Psets, Vsets] = mpt_reachSets(ctrl, X0, N, Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% structure       - either a sysStruct structure, or an explicit controller
% X0              - set of initial states (polytope or a polyarray)
% U0              - set of admissible inputs (polytope object)
% N               - number of steps over which to compute reachable sets
% Options.lowdim  - if set to 1, evolution of lower-dimensional sets will be
%                  computed. in such case, these lower-dimensional sets will be
%                  returned as V-represented polytopes (with vertices stored
%                  column-wise!). Default: 0
% Options.reduceV - if set to 1, V-represented lower-dimensional polytopes will
%                   be reduced to a minimal representation. Default: 0
% Options.verbose - level of verbosity (0/1/2)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Psets       - Sets of fully dimensional reachable states (polyarray)
% Vsets       - Lower-dimensional reachable sets (returned as a cell array, each
%               element represents one polytope in V-representation with
%               vertices stored row-wise).
%
% see also POLYTOPE/RANGE, MPT_REACHXU, MPT_VERIFY
%

% Copyright is with the following author(s):
%
% (C) 2004-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch

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

error(nargchk(3,5,nargin));

if nargin<5,
    Options = [];
end

if ~isa(X0, 'polytope'),
    error('mpt_reachSets: X0 must be a polytope!');
end
if isa(structure, 'mptctrl')
    if ~isexplicit(structure),
        error('This function supports only explicit controllers!');
    end
    structure = struct(structure);
end
if ~isstruct(structure),
    error('mpt_reachSets: first argument must be a structure!');
end

if isa(U0, 'double'),
    if nargin>3,
        opt = N;
    end
    N = U0;
    if nargin==4,
        if isstruct(opt),
            Options = opt;
        end
    elseif nargin==5,
        if ~isa(N, 'double'),
            error('mpt_reachSets: N must be a scalar!');
        end
    end
elseif isa(U0, 'polytope')
    if nargin==4,
        if isstruct(N),
            Options = N;
        end
    elseif nargin==5,
        if ~isa(N, 'double'),
            error('mpt_reachSets: N must be a scalar!');
        end
    end
else
    error('Unknown type of third input argument!');
end
    

if ~isa(N,'double')
    error('mpt_reachSets: third argument must be the horizon!');
end

if ~isfield(Options,'verbose')
    Options.verbose = mptOptions.verbose;
end
if ~isfield(Options, 'lowdim'),
    % if set to 1, evolution of lower-dimensional sets will be computed. In such
    % case, these lower-dimensional sets will be returned as V-represented
    % polytopes (with vertices stored column-wise!)
    Options.lowdim = 0;
end
if ~isfield(Options, 'abs_tol'),
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options, 'lpsolver'),
    Options.lpsolver = mptOptions.lpsolver;
end
if ~isfield(Options, 'reduceV'),
    % if set to 1, V-represented lower-dimensional polytopes will be reduced to
    % a minimal representation
    Options.reduceV = 0;
end
if ~isfield(Options, 'onlyVrepr'),
    % if true, all reachable sets will be returned in V-representation
    Options.onlyVrepr = 0;
end
if ~isfield(Options, 'Xf'),
    % set of final states. if not empty, the algorithm will check if reachable
    % states intersect this set
    Options.Xf = [];
end

Options.checkXf = ~isempty(Options.Xf);

if Options.lowdim | Options.reduceV | Options.checkXf,
    % test if CDD is available
    if ~ismember(3, mptOptions.solvers.extreme),
        disp('CDD is not available on your computer, lower-dimensional sets will NOT be computed!');
        Options.lowdim = 0;
        Options.reduceV = 0;
    end
end

if Options.checkXf,
    if ~isa(Options.Xf, 'polytope'),
        error('Xf must be a polytope object!');
    elseif ~isfulldim(Options.Xf),
        error('Xf must be a fully dimensional polytope!');
    end
    [Hf, Kf] = double(Options.Xf);
    Xf = Options.Xf;
    eXf = extreme(Xf);
    Options.onlyVrepr = 1;
end

abs_tol = Options.abs_tol;
lpsolver = Options.lpsolver;

if isa(U0, 'polytope'),
    infbox = mptOptions.infbox;
    sysStruct = structure;
    if ~mpt_issysstruct(sysStruct)
        error('mpt_reachSets: First input must be a system structure if U0 is a polytope!');
    end
    
    if ~isfield(sysStruct, 'verified'),
        sysStruct = mpt_verifySysStruct(sysStruct);
    end
    
    [nx, nu, ny, ndyn] = mpt_sysStructInfo(sysStruct);
    if dimension(X0)~=nx,
        error(sprintf('X0 must be a %dD polytope!', nx));
    end
    if dimension(U0)~=nu,
        error(sprintf('U0 must be a %dD polytope!', nu));
    end
    InfX = polytope([eye(nx); -eye(nx)], infbox*ones(2*nx,1));
    InfU = polytope([eye(nu); -eye(nu)], infbox*ones(2*nu,1));
    Pinf = InfX * InfU;
    [Hinf, Kinf] = double(Pinf);
    
    if ~iscell(sysStruct.A),
        sysStruct = mpt_lti2pwa(sysStruct);
    end

    xbool = [];
    if isfield(sysStruct, 'xbool'),
        xbool = sysStruct.xbool;
    elseif isfield(sysStruct, 'data'),
        if isfield(sysStruct.data, 'MLD'),
            nxr = sysStruct.data.MLD.nxr;
            nx = sysStruct.data.MLD.nx;
            xbool = nxr+1:nx;
        end
    end
    
    Pi = polytope;
    for ii = 1:ndyn,
        gX = sysStruct.guardX{ii};
        gU = sysStruct.guardU{ii};
        gC = sysStruct.guardC{ii};
        Pi = [Pi polytope([Hinf; gX gU], [Kinf; gC])];
    end

    Ei = pelemfun(@extreme, Pi);
    if isa(U0, 'polytope'),
        if length(U0)>1,
            error('This function does not support polytope arrays in U0.');
        end
        if ~isfulldim(U0),
            error('U0 must be a fully dimensional polytope.');
        end
        U0e = extreme(U0);
    end
    
    emptypoly = mptOptions.emptypoly;
    Xreach = emptypoly;
    PreachN = {};
    VreachN = {};
    Vreach = {};
    
    [Hn, Kn] = double(Pi);

    justonedim = length(Pi)==1;

    if Options.onlyVrepr,
        if isa(X0, 'polytope'),
            if ~isfulldim(X0),
                error('Set of initial states must be a fully dimensional polytope!');
            elseif length(X0)>0,
                V0 = {};
                for ii = 1:length(X0),
                    V0{end+1} = extreme(X0(ii));
                end
            end
        else
            V0 = X0;
        end
        X0 = polytope;
        Options.lowdim = 1;
        Options.reduceV = 1;
    end

    for iN = 1:N,
        
        if mod(iN,5)==0 | iN==N,
            fprintf('Horizon %d\t', iN);
            if mod(iN,25)==0,
                fprintf('\n');
            end
        end

        X0new = emptypoly;
        V0new = {};
        
        for jj = 1:length(X0),
            if justonedim,
                X0int = X0(jj)*U0;
                whichdyn = 1;
            else
                [X0int, whichdyn] = sub_intersect(X0(jj)*U0, Pi, Ei, Hn, Kn, emptypoly, abs_tol, lpsolver);
            end
            if isfulldim(X0int),
                for ii = 1:length(X0int),
                    Ai = sysStruct.A{whichdyn(ii)};
                    Bi = sysStruct.B{whichdyn(ii)};
                    fi = sysStruct.f{whichdyn(ii)};
                    [Xn,Vn] = sub_reachXU(Ai, Bi, fi, X0int(ii), xbool);
                    X0new = [X0new Xn];
                    Xreach = [Xreach Xn];
                    if ~isfulldim(Xn) & ~isempty(Vn),
                        if Options.reduceV,
                            Vn = sub_reduceV(Vn);
                        end
                        Vreach{end+1} = Vn;
                        V0new{end+1} = Vn;
                    end
                end
            end
        end
        
        if Options.lowdim,
            for jj = 1:length(V0),
                
                % first make a polytope in the extended [X U] space using
                % V-representation
                VV = V0{jj};
                VVn = [];
                for kk = 1:size(U0e, 1),
                    VVn = [VVn; VV repmat(U0e(kk,:), size(VV, 1), 1)];
                end
                V0U0 = VVn';
                
                % then continue with reachability analysis
                if justonedim,
                    X0int = {V0U0'};
                    whichdyn = 1;
                else
                    [X0int, whichdyn] = sub_intersect(V0U0, Pi, Ei, Hn, Kn, emptypoly, abs_tol, lpsolver);
                end
                if ~isempty(X0int),
                    for ii = 1:length(whichdyn),
                        Ai = sysStruct.A{whichdyn(ii)};
                        Bi = sysStruct.B{whichdyn(ii)};
                        fi = sysStruct.f{whichdyn(ii)};
                        [Xn,Vn] = sub_reachXU(Ai, Bi, fi, X0int{ii}, xbool);
                        if ~isempty(Vn),
                            if Options.checkXf,
                                dosect = sub_intersect(Vn', Xf, eXf, Hf, Kf, emptypoly, abs_tol, lpsolver);
                                if dosect,
                                    Psets = iN;
                                    Vreach{end+1} = Vn;
                                    Vsets = Vreach;
                                    fprintf('\n');
                                    return
                                end
                            end
                            if Options.reduceV,
                                Vn = sub_reduceV(Vn);
                            end
                            Vreach{end+1} = Vn;
                            V0new{end+1} = Vn;
                        end
                    end
                end
            end
        end
        
        X0 = reduceunion(X0new);
        V0 = V0new;
        PreachN{iN} = X0new;
        VreachN{iN} = V0new;
    end
    Psets = Xreach;
    Vsets = Vreach;
%     if ~isempty(Vreach),
%         fprintf('\nReachable set contains lower dimensional polytopes, V-representation returned as second argument\n');
%     end
    fprintf('\n');
    return
end

if isfield(structure, 'A'),
    % sysStruct case
    sysStruct = structure;
    if ~isfield(sysStruct, 'verified')
        sysStruct = mpt_verifySysStruct(structure);
    end
    
    [nx,nu,ny,ndyn] = mpt_sysStructInfo(sysStruct);
    
    if ~isfield(Options,'U'),
        Options.U = zeros(nu,1);
    end
    
    if ndyn<2,
        % only 1 dynamics - LTI system
        Pn = sysStruct.Pbnd;
        ACL{1} = sysStruct.A;
        if isfield(sysStruct,'f'),
            FCL{1} = sysStruct.f  + sysStruct.B*Options.U;
        else
            FCL{1} = zeros(nx,1)  + sysStruct.B*Options.U;
        end
    else
        % PWA system
        Pn = mptOptions.emptypoly;
        alloweddyn = [];
        for idyn = 1:ndyn,
            pp = polytope([sysStruct.guardX{idyn} sysStruct.guardU{idyn}*Options.U], sysStruct.guardC{idyn});
            if isfulldim(pp),
                alloweddyn = [alloweddyn idyn];
                Pdyn = projection(pp,1:nx);
                Pn = [Pn Pdyn];
            end
        end
        [Pn, kept] = reduceunion(Pn);
        if ~isfulldim(Pn),
            error('mpt_reachSets: no dynamics associated to zero input!');
        end

        keptdynamics = find(kept==1);
        
        dyns = alloweddyn(keptdynamics);
        for idyn = 1:length(dyns), 
            ACL{idyn} = sysStruct.A{dyns(idyn)};
            FCL{idyn} = sysStruct.f{dyns(idyn)} + sysStruct.B{dyns(idyn)}*Options.U;
        end
    end
    
elseif isfield(structure, 'Pn')
    % ctrlStruct case
    ctrlStruct = structure;
    isvalid = mpt_isValidCS(ctrlStruct, struct('nowarnings',1));
    if ~isvalid,
        error('mpt_reachSets: input is not a valid controller structure!');
    end
    if ctrlStruct.overlaps == 1,
        disp('remove overlaps first by calling "ctrlStruct=mpt_removeOverlaps(ctrlStruct)"');
        error('mpt_reachSets: overlapping controllers not supported.');
    end
    sysStruct = ctrlStruct.sysStruct;
    probStruct = ctrlStruct.probStruct;
    % handle feedback prestabilization:
    if isfield(probStruct, 'FBgain'),
        FBgain = probStruct.FBgain;
    else
        FBgain = zeros(size(ctrlStruct.Fi{1}));
    end
    [nx,nu,ny,ndyn] = mpt_sysStructInfo(sysStruct);
    Pn = ctrlStruct.Pn;
    for ireg = 1:length(Pn),
        dyn = 0;
        if isfield(ctrlStruct, 'dynamics'),
            dyn = ctrlStruct.dynamics(ireg);
        end
        if dyn ~= 0
            % we know a-priori which dynamics corresponds to region "ireg"
            if iscell(sysStruct.A),
                A = sysStruct.A{dyn};
                B = sysStruct.B{dyn};
                f = sysStruct.f{dyn};
            else
                A = sysStruct.A;
                B = sysStruct.B;
                if isfield(sysStruct, 'f'),
                    f = sysStruct.f;
                else
                    f = zeros(nx,1);
                end
            end
            Fi = ctrlStruct.Fi{ireg};
            Gi = ctrlStruct.Gi{ireg};
            
            % autonomous system dynamics:
            % x+ = (A+B*Fi)*x + (B*Gi + f)
            ACL{ireg} = A + B*(Fi(1:nu,:) + FBgain(1:nu,:));
            FCL{ireg} = f + B*Gi(1:nu,:);
            
        elseif ~iscell(sysStruct.A),
            % LTI system
            A = sysStruct.A;
            B = sysStruct.B;
            if isfield(sysStruct, 'f'),
                f = sysStruct.f;
            else
                f = zeros(nx,1);
            end
            Fi = ctrlStruct.Fi{ireg};
            Gi = ctrlStruct.Gi{ireg};
            
            % autonomous system dynamics:
            % x+ = (A+B*Fi)*x + (B*Gi + f)
            ACL{ireg} = A + B*(Fi(1:nu,:) + FBgain(1:nu,:));
            FCL{ireg} = f + B*Gi(1:nu,:);

        else
            % PWA system, link dynamics to regions...
            Fi = ctrlStruct.Fi{ireg};
            Gi = ctrlStruct.Gi{ireg};

            [x,R] = chebyball(Pn(ireg));            % compute center of the chebyshev's ball
            Acell = [];
            Fcell = [];
            for jj=1:length(sysStruct.A)          % go through all dynamics description
                if max(sysStruct.guardX{jj}*x+sysStruct.guardU{jj}*((Fi(1:nu,:)+FBgain(1:nu,:))*x+Gi(1:nu,:))-sysStruct.guardC{jj})<abs_tol, 
                    % check which dynamics is active in the region
                    ACL{ireg}=sysStruct.A{jj} + sysStruct.B{jj}*(Fi(1:nu,:) + FBgain(1:nu,:));
                    FCL{ireg}=sysStruct.f{jj} + sysStruct.B{jj}*Gi(1:nu,:);
                    break
                end
            end
            if isempty(ACL{ireg}) | isempty(FCL{ireg}),
                error('Faulty partition: Region could not be linked to any dynamic !!')
            end
        end
    end
        
else
    error('mpt_reachSets: unknown structure!');
end


%---------------------------------------------------------------
emptypoly = mptOptions.emptypoly;
Xreach = emptypoly;
PreachN = {};
VreachN = {};
Vreach = {};
V0 = {};

%enumerate extreme points (for faster computation of intersections of lower
%dimensional polytopes)
En = cell(1, length(Pn));
for ii=1:length(Pn),
    En{ii} = extreme(Pn(ii));
end

[Hn, Kn] = double(Pn);

if Options.onlyVrepr,
    if isa(X0, 'polytope'),
        if ~isfulldim(X0),
            error('Set of initial states must be a fully dimensional polytope!');
        elseif length(X0)>0,
            V0 = {};
            for ii = 1:length(X0),
                V0{end+1} = extreme(X0(ii));
            end
        end
    else
        V0 = X0;
    end
    X0 = polytope;
    Options.lowdim = 1;
    Options.reduceV = 1;
end

for iN = 1:N,
    
    if mod(iN,5)==0 | iN==N,
        fprintf('Horizon %d\t', iN);
        if mod(iN,25)==0,
            fprintf('\n');
        end
    end
    
    X0new = emptypoly;
    V0new = {};
    for jj = 1:length(X0),
        [X0int, whichdyn] = sub_intersect(X0(jj), Pn, En, Hn, Kn, emptypoly, abs_tol, lpsolver);
        if isfulldim(X0int),
            for ii = 1:length(X0int),
                Xi = X0int(ii);
                Ai = ACL{whichdyn(ii)};
                fi = FCL{whichdyn(ii)};
                if det(Ai) < abs_tol,
                    % mapping will be lower-dimensional
                    exi = extreme(Xi);
                    Vn = Ai*exi' + repmat(fi, 1, size(exi', 2));
                    Vn = Vn';
                    if Options.reduceV,
                        Vn = sub_reduceV(Vn);
                    end
                    Vreach{end+1} = Vn;
                    V0new{end+1} = Vn;
                else
                    Xn = range(Xi, Ai, fi);
                    X0new = [X0new Xn];
                end
            end
        end
    end

    if Options.lowdim,
        for jj = 1:length(V0),
            [V0int, whichdyn] = sub_intersect(V0{jj}', Pn, En, Hn, Kn, emptypoly, abs_tol, lpsolver);
            if ~isempty(V0int)
                if ~iscell(V0int),
                    V0int = {V0int};
                end
                for ii = 1:length(V0int),
                    exi = V0int{ii};
                    Ai = ACL{whichdyn(ii)};
                    fi = FCL{whichdyn(ii)};
                    Vn = Ai*exi + repmat(fi, 1, size(exi, 2));
                    Vn = Vn';
                    
                    if Options.checkXf,
                        dosect = sub_intersect(Vn', Xf, eXf, Hf, Kf, emptypoly, abs_tol, lpsolver);
                        if dosect,
                            Psets = iN;
                            Vreach{end+1} = Vn;
                            Vsets = Vreach;
                            fprintf('\n');
                            return
                        end
                    end
                    if Options.reduceV,
                        Vn = sub_reduceV(Vn);
                    end
                    Vreach{end+1} = Vn;
                    V0new{end+1} = Vn;
                end
            end
        end
    end
    
    % kick out identical regions
    X0 = reduceunion(X0new);
    V0 = V0new;
    
    Xreach = [Xreach X0];
    PreachN{iN} = X0new;
    VreachN{iN} = V0new;
end

Psets = Xreach;
Vsets = Vreach;

% if ~isempty(Vreach),
%     fprintf('\nReachable set contains lower dimensional polytopes, V-representation returned as second argument\n');
% end
fprintf('\n');

return
%---------------------------------------------------------------



%-----------------------------------------
function [I, which] = sub_intersect(P, Pn, En, Hn, Kn, emptypoly, abs_tol, lpsolver)

lenPn = length(Pn);

if lenPn<1,
    I = emptypoly;
    which = [];
    return
elseif lenPn == 1,
    if isa(P, 'polytope'),
        %I = P & Pn;
        [H,K] = double(P);
        [Hn,Kn] = double(Pn);
        I = polytope([H; Hn], [K; Kn]);
    else
        which = [];
        I = [];
        for ii = 1:size(P, 2);
            if all(Hn*P(:, ii) - Kn < abs_tol),
                % some points lie in the interior of Pn => intersection
                if nargout==1,
                    I = 1;
                    return
                else
                    HK = cddmex('hull', struct('V', P'));
                    [Hi, Ki] = double(Pn);
                    A = [HK.A; Hi]; B = [HK.B; Ki];
                    V = cddmex('extreme', struct('A', A, 'B', B, 'lin', HK.lin));
                    I = V.V';
                    which = 1;
                end
            end
        end
        if nargout==1,
            I = 0;
        end
        return
        
        % try to find a separating hyperplane between P and Pn
        Ve = En;
        nve = size(Ve,1);
        npe = size(P, 2);
        nx = dimension(Pn);
        f = -ones(1,nx+1);
        A = [Ve -ones(nve,1); -P' ones(npe, 1); zeros(1, nx) 1];
        B = [zeros(nve+npe, 1); 1];
        [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,A,B,[],[],[],lpsolver);

        sephp = ~all(xopt==0) & ~any(xopt==1e9);
        if sephp
            % separating hyperplane exists, polytopes do not intersect
            if nargout==1,
                I = 0;
                return
            else
                I = [];
                which = [];
                return
            end
        else
            if nargout==1,
                I = 1;
                return
            end
            % separating hyperplane does not exist, polytopes do intersect
            HK = cddmex('hull', struct('V', P'));
            [Hi, Ki] = double(Pn);
            A = [HK.A; Hi]; B = [HK.B; Ki];
            V = cddmex('extreme', struct('A', A, 'B', B, 'lin', HK.lin));
            I = V.V';
        end
    end
    which = 1;
    return
else
    if isa(P, 'polytope'),
        [H,K] = double(P);
        [Hn,Kn] = double(Pn);
        I = emptypoly;
        which = [];
        for ii=1:lenPn,
            if dointersect(P, Pn(ii));
                % intersection is fully dimensional, compute it
                Hint = [H; Hn{ii}];
                Kint = [K; Kn{ii}];
                int = polytope(Hint, Kint);
                I = [I int];
                which = [which ii];
            end
        end
        
    else
        I = {};
        which = [];
        nx = dimension(Pn);
        f = -ones(1,nx+1);
        for ii=1:lenPn,
            Ve = En{ii};
            Hi = Hn{ii};
            Ki = Kn{ii};
            
            nve = size(Ve,1);
            npe = size(P, 2);
            
            sephp = 1;
            %             for kk = 1:npe,
            %                 if all(Hi*P(:, kk) - Ki < abs_tol),
            %                     point lies in interior ot Pn(ii), intersection exists
            %                     sephp = 0;
            %                     break
            %                 end
            %             end
            if sephp,
                % we still need to check if intersection exists
                A = [Ve -ones(nve,1); -P' ones(npe, 1); zeros(1, nx) 1];
                B = [zeros(nve+npe, 1); 1];
                [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,A,B,[],[],[],lpsolver);
                sephp = ~all(xopt==0) & ~any(xopt==1e9);
            end

            if ~sephp,
                % separating hyperplane DOES NOT exist, polytopes DO intersect
                if nargout==2,
                    HK = cddmex('hull', struct('V', P'));
                    A = [HK.A; Hi]; B = [HK.B; Ki];
                    V = cddmex('extreme', struct('A', A, 'B', B, 'lin', HK.lin));
                    if ~isempty(V.V),
                        I{end+1} = V.V';
                        which = [which ii];
                    end
                else
                    I = 1;
                    return
                end
            end
        end
    end
end


%=======================================================================
function Vr = sub_reduceV(V)

Vr = [];
vert.V = V;
%vn = cddmex('v_hull_extreme', vert);
vn = cddmex('reduce_v', vert);
if isstruct(vn)
    if ~isempty(vn.V),
        Vr = vn.V;
    end
end


%=======================================================================
function [Pxn,V] = sub_reachXU(A, B, f, XU, xbool),

if nargin<5,
    xbool = [];
end

% EXPERIMENTAL
if false
    nu = size(B, 2);
    nx = size(A, 1);
    [H, K] = double(XU);
    Hx = H(:, 1:nx); z = sum(Hx.*Hx, 2); Hx = Hx(z>0, :); Kx = K(z>0, :);
    Hu = H(:, nx+1:nx+nu); z = sum(Hu.*Hu, 2); Hu = Hu(z>0, :); Ku = K(z>0, :);
    X = polytope(Hx, Kx);
    U = polytope(Hu, Ku);
    options.userange = 1;
    [Pxn, V] = mpt_reachXU(X, U, A, B, f, options);
    return
end

if isa(XU, 'polytope'),
    Exu = extreme(XU);
else
    Exu = XU;
end
if isempty(Exu)
    error('Polytope in extended space is not fully dimensional!');
end
Exn = [A B] * Exu';
Exn = Exn + repmat(f, 1, size(Exn, 2));
V = Exn';
if ~isempty(xbool),
    % handle boolean states
    nx = size(V,2);
    Vnonbool = V(:, setdiff(1:nx, xbool));
    Vnew = [];
    for jj=1:size(V,1),
        V1 = V(jj, :);
        Vnew = [Vnew; V1];
        for ii=1:length(xbool),
            V2 = V1;
            if V2(xbool(ii)) == 1,
                V2(xbool(ii))=0.99;
            else
                V2(xbool(ii)) = 0.01;
            end
            Vnew = [Vnew; V2];
        end
    end
    V = Vnew;
end
if isa(XU, 'polytope'),
    Pxn = polytope(V);
    [xcheb, rcheb] = chebyball(Pxn);
    if ~isfulldim(Pxn) | rcheb == 1e9,
        Pxn = polytope;
    end
else
    Pxn = [];
end
