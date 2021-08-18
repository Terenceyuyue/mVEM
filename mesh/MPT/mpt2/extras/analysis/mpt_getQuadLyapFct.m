function [lyapP,decay,feasible,details]=mpt_getQuadLyapFct(ctrl,Options)
%MPT_GETQUADLYAPFCT Computes common Lyapunov function for PWA system
%
% [lyapP,decay,feasible]=mpt_getQuadLyapFct(ctrlStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% This function attempts to compute a quadratic Lyapunov function V(x)=x'lPx
% which guarantees exponential stability.
% PWQ(x) = x'LPx
% PWQ(x(k+1)) - PWQ(x(k)) <= rho * x^2 
%
% (i.e. rho must be negative to guarantee exponential stability)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl                   - Explicit controller (MPTCTRL objects)
% Options.decayEllipsoid      
%       - If set to 1, decay rate of PWQ function is constrained to be negative 
%         over ellipsoidal regions(reduces number of free variables in LMI)
%         If set to 0, decay rate is constrained to be negative 
%         over polytopic regions(more free variables in LMI)
% Options.posEllipsoid        
%       - If set to 1, positivity of PWQ function is constrained to be valid  
%         over ellipsoidal regions(reduces number of free variables in LMI)
%         If set to 0, positivity is constrained to be valid  
%         over polytopic regions(more free variables in LMI)
% Options.enforcePositivity   
%       - If set to zero, positivity constraints for PWQ function are not 
%         included in LMI (reduces constraints and computation time)
%         Post computation verification is performed to check if postivity holds
%         If not, positivity constraints are added and solution recomputed.
% Options.abs_tol       - Absolute tolerance
% Options.lpsolver      - Which LP solver to use (see help mpt_solveLP for details)
% Options.epsilon       - This is a tolerance factor which is introduced to turn
%                           LMI inequalities into strict inequalities.
% Options.debug_level   - If this is set to 1, the solution provided by the LMI
%                         solver will be double-checked manually. We strongly
%                         advise to set this to 1, since we've experienced
%                         numerous numerical issues with certain LMI solvers.
% Options.nicescale     - This will add additional constraints to obtain a nicely scaled
%                         Lyapunov function. There is no benefit in this apart from 
%                         getting nicer plots afterwards. Therefore this is switched 
%                         off by default
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% lyapP      - Quadratic Lyapunov function: Q(x)=x'*lyapP{r}*x
% decay      - the maximum Lyapunov decay rate over the partition
%              (is this is greater than zero, stability cannot be guaranteed)
%              The Lyapunov value decrease Delta V <= decay * ||x||^2
% feasible   - 1: stable or 0: no statement about stability possible
%
%
% see also MPT_GETPWQLYAPFCT, MPT_GETSTABFEEDBACK, MPT_GETCOMMONLYAPFCT

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
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


%This function returns a boolean on stability 
%the system dynamics are Giiven by x(k+1)=Ax(k)+Bu(k)
%over a polyhedron Hx<=K in which the feedback law u=Fx+Gi is optimal
%The function computes the maximum of the Lyapunov difference Giiven by:
% deltaP=x'Px-[(A+BF)x+BGi]'P[(A+BF)x+BGi]
% if the maximum is neGiative = > the system is stable
%
% Here A and B are cell arrays containing the vertices of the polytope...

error(nargchk(1,2,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if (nargin<2),
    Options=[];
end

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl)
        error('This function supports only explicit controllers!');
    end
    ctrlStruct = struct(ctrl);
else
    ctrlStruct = ctrl;
end

if ~isfield(Options,'debug_level'),
    Options.debug_level=mptOptions.debug_level;
end
verifySolution=Options.debug_level;
if ~isfield(Options,'epsilon'),
    Options.epsilon=mptOptions.abs_tol;   %epsilon is added to guarantee that a quadratic lower bound on PWQ exists
end
epsilon=Options.epsilon;
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'nicescale'),
    Options.nicescale=0;
end

if ~isfield(Options,'decayEllipsoid'),
    Options.decayEllipsoid=0;
end

if ~isfield(Options,'epsilon'),
    Options.epsilon=mptOptions.abs_tol;   %epsilon is added to guarantee that a quadratic lower bound on PWQ exists
end
epsilon=Options.epsilon;

if ~isfield(Options,'posEllipsoid'),
    Options.posEllipsoid=0;
end

if ~isfield(Options,'enforcePositivity'),
    Options.enforcePositivity=0;
end

start=cputime;


lpsolver=Options.lpsolver;

% if ~mpt_isValidCS(ctrlStruct)
%     error('mpt_getQuadLyapFct: First argument has to be a valid controller structure! See mpt_control for details.');
% end

sysStruct = ctrlStruct.sysStruct;
Pn = ctrlStruct.Pn;
Fi = ctrlStruct.Fi;
Gi = ctrlStruct.Gi;

verOptions.verbose=0;
if ~isfield(sysStruct,'verified')
    sysStruct = mpt_verifySysStruct(sysStruct, verOptions);
end

if isfield(ctrlStruct.probStruct, 'FBgain'),
    % handle pre-stabilization with feedback
    FBgain = ctrlStruct.probStruct.FBgain;
else
    FBgain = zeros(size(Fi{1}));
end

if (iscell(sysStruct.A))
    if(isfield(sysStruct,'Aunc') & ~isempty(sysStruct.Aunc))
        error('Cannot handle PWA systems with polytopic uncertainty')
    end
    pwasystem=1;
    Acell = sysStruct.A;
    Bcell = sysStruct.B;
    if length(sysStruct.A)~=length(Pn)
        disp('In PWA case, each dynamics has to be associated with exactly one region!');
        disp('Linking dynamics to regions...');
        nu=size(sysStruct.B{1},2);
        cc=0;
        Acell{1}=sysStruct.A{1}; %enter dummy values
        Bcell{1}=sysStruct.B{1}; %enter dummy values
        ABFcell = {};
        BGcell = {};
        for ii=1:length(Pn)
            [x,R] = chebyball(Pn(ii));            % compute center of the chebyshev's ball
            for jj=1:length(sysStruct.A)          % go through all dynamics description
                if max(sysStruct.guardX{jj}*x+sysStruct.guardU{jj}*((Fi{ii}(1:nu,:)+FBgain(1:nu,:))*x+Gi{ii}(1:nu,:))-sysStruct.guardC{jj})<Options.abs_tol,    % check which dynamics is active in the region
                    ABFcell{ii} = sysStruct.A{jj}+sysStruct.B{jj}*(Fi{ii}(1:nu,:) + FBgain(1:nu,:));
                    BGcell{ii} = sysStruct.B{jj}*Gi{ii}(1:nu,:) + sysStruct.f{jj};
                end
            end
            if(isempty(ABFcell{ii}))
                error('Faulty partition: Region could not be linked to any dynamic !!')
            end
        end
    else
        ABFcell = sysStruct.A;
        BGcell = sysStruct.f;
    end
        
elseif ~isfield(sysStruct,'Aunc') | ~isfield(sysStruct,'Bunc')
    %no polytopic uncertainty 
    Acell{1}=sysStruct.A;
    Bcell{1}=sysStruct.B;
    pwasystem=0;
else
    %polytopic uncertainty
    Acell=sysStruct.Aunc;
    Bcell=sysStruct.Bunc;
    pwasystem=0;
end

if isfield(sysStruct,'noise') & mpt_isnoise(sysStruct.noise)
    if isa(sysStruct.noise, 'polytope'),
        [Hnoise,Knoise]=double(sysStruct.noise);
        noiseVertex=extreme(sysStruct.noise);
    else
        % remember that V-represented noise has vertices stored column-wise,
        % therefore we need to convert them to a row-wise setup, since that's
        % what extreme() does.
        noiseVertex = sysStruct.noise';
    end
    havenoise = 1;
else
    havenoise = 0;    
    Hnoise=[];
    Knoise=[];    
    noiseVertex=zeros(1,size(Acell{1},2));
end

nu=size(Bcell{1},2);
nx=size(Acell{1},2);

try
    yalmip('clear') ;       %initialize yalmip
catch
    error('mpt_getQuadLyapFct: You need to download and install Yalmip for this function to work');
end

if havenoise,
	Hn={};
	Kn={};
	Fi={};
	Gi={};
    if iscell(sysStruct.C),
        error('PWA systems with additive noise are not supported.');
    end
    %extract RPI set from region 1
    X=polytope([sysStruct.C*eye(nx);-sysStruct.C*eye(nx)],[sysStruct.ymax;-sysStruct.ymin]);
    ctr=0; index=[];
    for j=1:length(Pn)
        for i=1:length(Acell)
            if(isinside(Pn(j),zeros(nx,1)))
                index=[index j];
                ctr=ctr+1;
                AA{ctr}=Acell{i}+Bcell{i}*(ctrlStruct.Fi{j}(1:nu,:) + FBgain(1:nu,:));
            else
                [HH,KK]=double(Pn(j));
                Hn{end+1}=HH;
                Kn{end+1}=KK;
                Fi{end+1}=ctrlStruct.Fi{j};
                Gi{end+1}=ctrlStruct.Gi{j};
            end
        end
    end
    if(length(index)>1)
        error('The case where multiple regions contain the origin is not treated by this function... yet.')
    end
    
    Oinf=mpt_infset(AA,X,[],sysStruct.noise);
    if(~isfulldim(Oinf))
        error('No robust invariant subset of region 1 around the origin exists')
    end
    Pnew=Pn(index)\Oinf;
    if(isfulldim(Pnew))
        for i=1:length(Pnew)
            [HH,KK]=double(Pnew(index));
            Hn{end+1}=HH;
            Kn{end+1}=KK;
            Fi{end+1}=ctrlStruct.Fi{index};
            Gi{end+1}=ctrlStruct.Gi{index};
        end
    end

else
    [Hn, Kn] = pelemfun(@double, Pn);
    Fi=ctrlStruct.Fi;
    Gi=ctrlStruct.Gi;
end

P = sdpvar(nx,nx);
rho=sdpvar(1,1);
ctr=0;

myprog=set([]);


for i=1:length(Hn)
    m = size(Hn{i},1); %facets
    HK = [-Hn{i} Kn{i}];
    F=Fi{i}(1:nu,:);
    G=Gi{i}(1:nu,:);
    fbgain = FBgain(1:nu,:);
    
    for j=1:length(Acell)
        for k=1:size(noiseVertex,1)
            noise=noiseVertex(k,:)';
            
            ctr=ctr+1;
            if pwasystem,
                ABF = ABFcell{i};
                BG = BGcell{i};
            else
                A=Acell{j};
                B=Bcell{j};
                ABF = A + B*(F + fbgain);
                BG = B*G;
            end
            
            if(all(Hn{i}*zeros(nx,1)<=Kn{i}) & all(noise==0))
                N{ctr} = 0;
                %%W = P-(A+B*F)'*P*(A+B*F)-rho*eye(nx);
                W = P-ABF'*P*ABF-rho*eye(nx);
                if(any(abs(G)>mptOptions.abs_tol))
                    error('If the origin is an equilibrium, the input for the region containing the origin may not have an affine term')
                end
            elseif(all(Hn{i}*zeros(nx,1)<=Kn{i}) & any(noise~=0))
                error('A PWA system containing the origin subject to noise CANNOT be asymptotically stable.')
            else
                N{ctr} = sdpvar(m,m);                
                N{ctr} = N{ctr}-diag(diag(N{ctr}));  % by Johan Loefberg
                %%HGam = [-((A+B*F)'*P*(A+B*F)-P)-rho*eye(nx)     -(G'*B'*P*(A+B*F))'-(noise'*P*(A+B*F))' ;...
                %%           -(G'*B'*P*(A+B*F))-noise'*P*(A+B*F)  -G'*B'*P*B*G-noise'*P*noise];
                temp = -(BG'*P*ABF)'-(noise'*P*ABF)'; % by Johan Loefberg
                HGam = [-(ABF'*P*ABF-P)-rho*eye(nx)    temp ;...
                             temp'                     -BG'*P*BG-noise'*P*noise];
                W = HGam - HK' * N{ctr} * HK ;
                indicies = find(triu(ones(m,m)-eye(m)));   % by Johan Loefberg
                myprog = myprog + set(N{ctr}(indicies)>0); % by Johan Loefberg
            end
            myprog = myprog + set(W>0);
        end
    end
end


myprog = myprog + set(rho>epsilon);
myprog = myprog + set(P>0);
myprog = myprog + set(P(:)<100);

details.setup_time=cputime-start;
start=cputime;

if ~isempty(mptOptions.sdpsettings),
    options=mptOptions.sdpsettings;
else
    options=sdpsettings('Verbose',double(Options.verbose>1));
end
solution = solvesdp(myprog,[],[],options);

details.solver_time=cputime-start;
start=cputime;

%analyze solver output
if(solution.problem==4) 
    res=checkset(myprog);
    if min(res)>0,
        infeasible=0;
        disp(['mpt_getQuadLyapFct: Numerical problems with LMI solver but solution is feasible.']);
    else
        infeasible=1;   %problem occured
        disp(['mpt_getQuadLyapFct: Numerical problems with LMI solver:  residual ' num2str(min(res)) ' (should be positive).'])
    end
elseif(solution.problem>0) 
    infeasible=1;   %problem occured
    res=checkset(myprog);
    disp(['mpt_getQuadLyapFct: LMI solver thinks problem is infeasible:  residual ' num2str(min(res)) ' (should be positive).'])
elseif(solution.problem<0)
    error('mpt_getQuadLyapFct: There is a problem with your LMI solver ! Go to the SeDuMi directory and type ''make''.')
else
    infeasible=0;
    disp('mpt_getQuadLyapFct: LMI solver found feasible solution.')
end


lyapunovP=double(P);
lyapP=lyapunovP/norm(lyapunovP);
decay=-double(rho);


if(infeasible==0 & decay<=0)
    feasible=1;
    if Options.verbose > 0
        disp(['The decay rate is:  ' num2str(decay) ' (should be negative)'])
    end
else 
    feasible=0;
end

details.postproc_time=cputime-start;

return
