function ctrlStruct = mpt_boolMinTime(sysStruct, probStruct, Options)
%MPT_BOOLMINTIME Computes minimum time controller for systems with discrete inputs
%
% ctrlStruct=mpt_boolMinTime(sysStruct,probStruct)
% ctrlStruct=mpt_boolMinTime(sysStruct,probStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the solution of a minimum-time problem for a given PWA system
%       x(k+1) = A_i x(k) + B_i u(k) + f_i
%       y(k)   = C_i x(k) + D_i u(k) + g_i
%       for i such that guardX(i) x(k) + guardU(i) u(k) <= guardC(i)  
%   s.t.
%       (ymin, ymax, umin, umax, dumin, dumax)
%
% All inputs are assumed to be discrete as defined in sysStruct.Uset
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct         - System structure in the sysStruct format
% probStruct        - Problem structure in the probStruct format
%
% Options.verbose   - Level of verbosity (see help mpt_init for more details)
% Options.details   - If set to 1, solution of each iteration is stored in the
%                       details fields of the resulting controller structure 
%                       (0 by default)
% Options.lowmem    - defines memory saving mode
%                       0 - consumes a lot of memory but runs fast (default)
%                       1 - slight memory saving 
%                       2 - heavy memory saving (slow computation)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% ctrlStruct   - Controller structure
%
% see also MPT_CONTROL, MPT_OPTCONROLPWA
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(2,3,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3,
    Options = [];
end

if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'maxiterations'),
    Options.maxiterations=100;
end

if ~isfield(sysStruct,'verified') | ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    [sysStruct,probStruct]=mpt_verifySysProb(sysStruct,probStruct,verOpt);
end

Options.noNoiseOnTset=1;

if ~iscell(sysStruct.A),
    % convert LTI system to PWA format
    sysStruct = mpt_lti2pwa(sysStruct);
end

if probStruct.tracking==1
    error('Tracking with deltaU formulation is not allowed! Please use "probStruct.tracking=2".');
end
if isfield(sysStruct, 'dumode'),
    error('Sorry, deltaU constraints are not supported!');
end
if isfield(probStruct,'Rdu')
    if ~all(probStruct.Rdu==0),
        error('mpt_optMixedCtrl: non-zero probStruct.Rdu is not supported in this version!');
    end
end

[nx,nu,ny,nPWA,nbool,ubool,intInfo] = mpt_sysStructInfo(sysStruct);

if nbool~=nu,
    error(['mpt_boolMinTime: All inputs must be integer! sysStruct.Uset has to be a cell with ' num2str(nu) ' elements.']);
end

if probStruct.Tconstraint==2 & isfulldim(probStruct.Tset),
    userTset = 1;
    inittarget = probStruct.Tset;
else
    userTset = 0;
    % use the Pbnd polytope as initial target set
    [Pn,dynamics] = sub_getInvSet(sysStruct, probStruct, zeros(nu,1),Options);
    inittarget = Pn;
end

nPWA = length(sysStruct.A);

Pdynamics = {};

%orig_gU = sysStruct.guardU;
orig_gC = sysStruct.guardC;

% determine parts of state-space where dynamics are active
Pdynamics = intInfo.Pdyn;

[bndA, bndb] = double(sysStruct.Pbnd);
pureX_gX = {}; pureX_gU = {}; pureX_gC = {};
pureU_gX = {}; pureU_gU = {}; pureU_gC = {};
for ii=1:nPWA,
    gX = sysStruct.guardX{ii};
    gC = sysStruct.guardC{ii};
    gU = sysStruct.guardU{ii};
    zerorows = [];
    for jj=1:size(gX,1)
        if all(gX(jj,:)==0),
            zerorows = [zerorows; jj];
        end
    end
    nonzerorows = setdiff(1:size(gX,1),zerorows);
    gXnz = gX(nonzerorows,:);
    gCnz = gC(nonzerorows,:);
    gUnz = gU(nonzerorows,:);
    
    gXz = gX(zerorows,:);
    gUz = gU(zerorows,:);
    gCz = gC(zerorows,:);
    pureX_gX{ii} = gXnz;
    pureX_gU{ii} = gUnz;
    pureX_gC{ii} = gCnz;
    pureU_gX{ii} = gXz;
    pureU_gU{ii} = gUz;
    pureU_gC{ii} = gCz;
end
            
if iscell(sysStruct.Uset),
    if length(sysStruct.Uset)~=nu,
        error(['mpt_boolMinTime: All inputs must be integer! sysStruct.Uset has to be a cell with ' num2str(nu) ' elements.']);
    end
elseif nu~=1,
    error('mpt_boolMinTime: All inputs must be integer! Please set sysStruct.Uset');
end

if any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax)),
    dUconstraints = 1;
else
    dUconstraints = 0;
end
if dUconstraints,
    error('mpt_boolMinTime: Constraints on delta U not supported in this version!');
end
if isfield(probStruct,'Rdu')
    if ~all(probStruct.Rdu==0),
        error('mpt_boolMinTime: non-zero probStruct.Rdu is not supported in this version!');
    end
end

Ucombs = mpt_allcombs(sysStruct.Uset);

emptypoly = polytope;
Pstore{1} = inittarget;
Ustore{1} = [];
for ii=1:length(inittarget),
    Ustore{1} = [Ustore{1}; zeros(1,nu)];
end

rdOptions = Options;
rdOptions.constructset=0;
rdOptions.reduce=0;
Options.reduce_output = 0;
Pmax = Pstore{1};
DynStore = {};
nR = 0;

starttime = cputime;

for iter = 2:Options.maxiterations+1,
    fprintf('new horizon: %d        \r', iter-1);
    Piter = emptypoly;
    Uiter = [];
    dynamics = [];
    expands=0;
    for possibleU = 1:size(Ucombs,1),
        U = Ucombs(possibleU,:)';
        DDyn=[];
        for ii=1:nPWA,
            % obtain index of dynamics associated to a given U
            if isempty(pureU_gU{ii}),
                DDyn=[DDyn ii];
            elseif all(pureU_gU{ii}*U - pureU_gC{ii}<=Options.abs_tol),
                DDyn = [DDyn ii];
            end
        end
        if isfield(sysStruct,'Uexcl'),
            if all(ismember(U(:)',sysStruct.Uexcl,'rows')),
                % exclude this U
                disp(['Excluding input ' mat2str(U')]);
                continue
            end
        end
        if isempty(DDyn),
            disp(['mpt_boolMinTime: no guardline associated to input ' num2str(U)]);
        end
        for dyn = DDyn,
            fprintf('U=%s, dynamics=%d      \r', mat2str(U'),dyn);
            
            A=sysStruct.A{dyn};
            B=sysStruct.B{dyn};
            f=sysStruct.f{dyn};
            
            % include constraints on 'y'
            C = sysStruct.C{dyn};
            D = sysStruct.D{dyn};
            g = sysStruct.g{dyn};
            
            PPfinal = Pdynamics{dyn};
            if probStruct.y0bounds,
                [Hf,Kf]=double(PPfinal);
                Hn = [Hf; C; -C];
                Kn = [Kf; sysStruct.ymax-D*U-g; -sysStruct.ymin+D*U+g];
                PPfinal = polytope(Hn,Kn);
            end
            
            for reg = 1:length(Pstore{iter-1})
                targetregion = Pstore{iter-1}(reg);
                
                deltaU = 0;
                if iter>2,
                    % control move valid in target region
                    targetU = (Ustore{iter-1}(reg,:))';
                    deltaU = U - targetU;
                    if any(deltaU<sysStruct.dumin) | any(deltaU>sysStruct.dumax),
                        % new control action violates slew constraints, continue
                        % with next sequence of inputs
                        continue
                    end
                end
                
                if ~probStruct.y0bounds,
                    [Hf,Kf]=double(targetregion);
                    Hn = [Hf; C; -C];
                    Kn = [Kf; sysStruct.ymax-D*U-g; -sysStruct.ymin+D*U+g];
                    targetregion = polytope(Hn,Kn);
                end
                
                % compute all states which can be driven from everywhere where
                % dynamics 'dyn' is active to 'targetregion' in one time step
                % using control move 'U'
                dmOptions.noReduce = 1;  % don't reduce the polytope, saves time
                [Pret,kr,feasible] = domain(targetregion,A,B*U+f,PPfinal,1,dmOptions);
                % if this set is not empty, we know there exists set of states which
                % can be driven into 'targetregion' using fixed U

                if feasible,
                    Pret = reduce(Pret);
                    
                    %                     if expands==0,
                    %                         R = regiondiff(Pret, Pmax, rdOptions);
                    %                         if isfulldim(R),
                    %                             expands=1;
                    %                         end
                    %                     end

                    %Piter = [Piter Pret];
                    Uiter = [Uiter; U'];
                    dynamics = [dynamics dyn];
                    [Piter,keep] = sub_mpt_expandlist(Pret, Piter, Options);
                    keepr = find(keep==1);
                    Uiter = Uiter(keepr,:);
                    dynamics = dynamics(keepr);
                    nR = nR+1;
                end
            end % reg
        end % dyn
    end % U
    [Pred,keep] = reduceunion(Piter);
    Pred = Piter;
    Ured = Uiter;
    dynred = dynamics;
    
    R = mldivide(Piter, Pmax, struct('simplecheck',1));
    if ~isfulldim(R),
        expands=0;
    else
        expands=1;
    end
        
    Pmax = [Pmax Pred];
    Pstore{iter} = Pred;
    Ustore{iter} = Ured;
    DynStore{iter} = dynred;
    if expands==0,
        break
    end
    fprintf('%d regions generated / %d regions in total\n\n',length(Piter),nR);
end %iter

if iscell(Pstore),
    len = length(Pstore);
else
    PP{1} = Pstore;
    Pstore = PP;
    len = 1;
end
Pn = [];
Pfinal = [];
Fi = {};
Gi = {};
Ai = {};
Bi = {};
Ci = {};
dynamics = [];
for ii=1:len,
    if ii==1 & userTset,
        continue
    end
    if ~isempty(DynStore{ii}),
        for reg=1:length(Pstore{ii}),
            Pn = [Pn Pstore{ii}(reg)];
            Fi{end+1} = zeros(nu,nx);
            Gi{end+1} = Ustore{ii}(reg,:)';
            Ai{end+1} = zeros(nx);
            Bi{end+1} = zeros(1,nx);
            Ci{end+1} = ii;
            dynamics = [dynamics DynStore{ii}(reg)];
        end
    end
end
ctrlStruct.sysStruct = sysStruct;
ctrlStruct.probStruct = probStruct;
ctrlStruct.Pn = Pn;
ctrlStruct.Pfinal = Pn;
ctrlStruct.Fi = Fi;
ctrlStruct.Gi = Gi;
ctrlStruct.Ai = Ai;
ctrlStruct.Bi = Bi;
ctrlStruct.Ci = Ci;
ctrlStruct.dynamics = dynamics;
ctrlStruct.overlaps = 1;
ctrlStruct.details.runTime = cputime-starttime;    


% ===========================================================================================
function [Pun,keep] = sub_mpt_expandlist(Pf,Pu,Options),
% given a polytope Pf and a polyarray Pu, using a subset check removes all fields of Pu
% which are covered by Pf

Options.reduce=0;        % for fast subset check
Options.constructset=0;  % for fast subset check
Options.elementwise=1;   % go elementwise, i.e. Pu(1)<=Pf, Pu(2)<=Pf, ... , Pu(n)<=Pf
expands = 1;
if ~isfulldim(Pu(1)),
    % if Pu is empty, returns Pf
    Pun=Pf;
    keep=1;
    return
end

PuExt=[Pu Pf];
lenPuExt = length(PuExt);
keep = [];
if lenPuExt>1,
    keep=(~le(Pu,Pf,Options))'; % returns indices of polytopes in Pu which are not a subset of Pf
    keep=[keep 1];
else
    keep=1;
end

Pun=PuExt(find(keep==1));    % the output consists of polytopes which are not a subset of Pf
return


% ===========================================================================================
function [Pn,dynamics]=sub_getInvSet(sysStruct, probStruct, Ustab, Options),


nPWA=length(sysStruct.A);
nx=length(sysStruct.A{1});
nu=size(sysStruct.B{1},2);

%obtain dynamics which contain origin
originindynamics=[];
if isfield(probStruct,'xref'),
    origin = probStruct.xref;
else
    origin = zeros(nx,1);
end
ctr=0;
for ii=1:nPWA
    if(all(sysStruct.f{ii}==0)) | isfield(probStruct,'xref') | isfield(probStruct,'uref')
        %.... othewise the origin would not be an equilibrium point
        if all(sysStruct.guardU{ii}==0) & max(sysStruct.guardX{ii}*origin - sysStruct.guardC{ii}) <= 0,
            originindynamics = [originindynamics ii];
        elseif(any(sysStruct.guardU{ii}~=0))
            tempP=polytope(sysStruct.guardU{ii},-sysStruct.guardX{ii}*origin + sysStruct.guardC{ii});
            if(isfulldim(tempP))
                originindynamics = [originindynamics ii];
            end
        end
    end
end

if(isempty(originindynamics))
    error('mpt_boolMinTime: No dynamic has the origin as an equilibrium point !! Aborting computation...');
end
if Options.verbose>=1,
    disp(['origin included in: ' num2str(originindynamics)]);
end

% here we compute the initial invariant target set for all dynamics which contain the origin in their interior
if Options.verbose>=1,
    disp('Computing target set');
end
userTset=0; % true if user provided the terminal set, false otherwise
%----------------------------------------------------------------------------
%Compute invariant target set for PWA system
%----------------------------------------------------------------------------

Fi = zeros(nu,nx);
P_CL = polytope;
for ii=1:length(originindynamics)
    dyn = originindynamics(ii);
    A_CL{ii}=sysStruct.A{dyn};
    f{ii}=sysStruct.f{dyn}+sysStruct.B{dyn}*Ustab;
    %%CONSTRUCT CONSTRAINT MATRICES FOR t=0
    Hx{ii}=[Fi;-Fi];                         %input constraints
    Kx{ii}=[sysStruct.umax;-sysStruct.umin];     %input constraints
    Hx{ii}=[Hx{ii}; sysStruct.C{dyn}*eye(nx);-sysStruct.C{dyn}*eye(nx)];
    Kx{ii}=[Kx{ii}; sysStruct.ymax;-sysStruct.ymin];
    Hx{ii}=[Hx{ii}; sysStruct.guardX{dyn}+sysStruct.guardU{dyn}*Fi];
    Kx{ii}=[Kx{ii}; sysStruct.guardC{dyn}];
    %P_CL=[P_CL polytope(Hx{ii},Kx{ii})];
    P_CL = [P_CL unitbox(2,2)];
end

%COMPUTE POSITIVE INVARIANT SUBSET OF PARTITION
[Pn,dynamics]=mpt_infsetPWA(P_CL,A_CL,f,sysStruct.noise,Options); %compute invariant target set
%Associate controllers to regions
if ~isfulldim(Pn),
    error('mpt_computePWATset: Invariant set is empty! Check your system definition.');
end
dynamics=originindynamics(dynamics);
