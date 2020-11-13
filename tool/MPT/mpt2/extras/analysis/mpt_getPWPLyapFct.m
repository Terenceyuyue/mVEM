function [dQ,feasible,drho,runtime]=mpt_getPWPLyapFct(ctrl,ndeg,Options)
%MPT_GETPWPLYAPFCT Calculates Piecewise Polynomial (sum of squares) Lyapunov function
%
%  [dQ,feasible,drho]=mpt_getPWPLyapFct(ctrl,ndeg,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function attempts to compute a piecewise higher order Sum of Squares Lyapunov function
% using YALMIP which guarantees exponential stability. The following is satisfied
% alpha * x^2 <= V(x(k))<= beta * x^2  
% V(x(k+1)) - V(x(k))<= - rho * x^2      
% (alpha,beta,rho>=0)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------  
% ctrl                        - Explicit controller (EXPCTRL object)
% ndeg                        - Degree of Lyapunov Function(should be even)
%
% Options.debug_level  - If this is set to 1, the solution provided by the LMI
%                        solver will be double-checked manually via SDP. We
%                        strongly advise to set this to 1, since we've experienced
%                        numerous numerical issues with certain LMI solvers.
%                        If this value is set to 2, an additional check via
%                        state trajectories will be performed. 
% Options.enforcePositivity   - If set to zero, positivity constraints for PWQ
%                               function are not included in LMI (reduces
%                               constraints and computation time) Post
%                               computation verification is performed to check
%                               if postivity holds. If not, positivity
%                               constraints are added and solution recomputed. 
% Options.abs_tol      - Absolute tolerance
% Options.lpsolver     - Which LP solver to use (see help mpt_solveLP for details)
% Options.epsilon      - This is a tolerance factor which is introduced to turn
%                        LMI inequalities into strict inequalities.
% Options.nicescale    - This will add additional constraints to obtain a nicely
%                        scaled Lyapunov function. There is no benefit in this
%                        apart from getting nicer plots afterwards. Therefore
%                        this is switched off by default
% Options.useTmap      - If set to true (default), transition map will
%                        be computed to rule out certain transitions
% Options.sphratio     - Gives factor which governs maximum number of separating
%                        hyperplanes computed in transition maps. Number of
%                        separating  hyperplnaes computed at each step is given
%                        by length(Pn)^2 / Options.ratio
%                        Default value is 20.
%                        Set this option to 0 if you don't want to impose any
%                        limit on number of separating hyperplanes.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% dQ         - cell array of coefficients of Lyapunov polytope
%              In order to see the output polynomial in a "nice" way, proceed as follows:
%              x1=sdpvar(1,1), x2=sdpvar(1,1), ... , xn=sdpvar(1,1)
%              xmon=monolist([x1 x2 ... xn], ndeg);
%              sdisplay(dQ{i}'*xmon)
%
%    
% feasible   - 1: stable or 0: no statement about stability possible
% drho       - the maximum Lyapunov decay rate over the partition
%              (is this is greater than zero, stability cannot be guaranteed)
%              The Lyapunov value decrease Delta V <= drho * ||x||^2
% runtime    - structure that stores "setup_time", "solver_time" and "postproc_time"
%
% see also MPT_GETPWQLYAPFCT, MPT_GETPWALYAPFCT, MPT_GETCOMMONLYAPFCT, MPT_GETSTABFEEDBACK
%

% ---------------------------------------------------------------------------
%   LITERATURE:
% ---------------------------------------------------------------------------
% IFAC World Congress, Prague, Czech Republic, 2005
% "A Survey on Stability Analysis of Discrete-Time Piecewise Affine Systems"
% P. Biswas, P. Grieder, J. Loefberg, M. Morari
% Also available from: http://control.ee.ethz.ch

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003-2005 Pascal Grieder, Automatic Control Laboratory, ETH Zurich       
%               grieder@control.ee.ethz.ch
% (C) 2004 Pratik Biswas       
%          pbiswas@stanford.edu

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
if (nargin<2 | ~isreal(ndeg) | ndeg<2),
    error('You must specify a degree >1 for the polynomial function')
end
if (nargin<3),
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
if ~isfield(Options,'Sproc'),
    %apply S-procedure using constants instead of SOS functions
    Options.Sproc=1;
end
if ~isfield(Options, 'useTmap'),
    % if set to 1, transition map will be computed to rule out certain
    % transitions
    Options.useTmap = 1;
end
if ~isfield(Options, 'sphratio'),
    sphratio = 20;
else
    sphratio = Options.sphratio;
end
if sphratio == 0,
    % to avoid "division by zero" warnings
    sphratio = 1e-9;
end

if ~isfield(Options,'is_invariant'),
    %is the PWA partition guaranteed to be invariant?
    if isa(ctrl, 'mptctrl')
        Options.is_invariant = isinvariant(ctrl);
    else
        Options.is_invariant = 0;
    end
end

verifySolution=Options.debug_level;

if ~isfield(Options,'epsilon'),
    Options.epsilon=mptOptions.abs_tol;   %epsilon is added to guarantee that a quadratic lower bound on PWP exists
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
    Options.nicescale=1;
end

if ~isfield(Options,'decayEllipsoid'),
    Options.decayEllipsoid=0;
end

if ~isfield(Options,'posEllipsoid'),
    Options.posEllipsoid=0;
end

if ~isfield(Options,'enforcePositivity'),
    Options.enforcePositivity=0;
end


lpsolver=Options.lpsolver;


starttime=cputime; %measure runtime

sysStruct = ctrlStruct.sysStruct;
Pn = ctrlStruct.Pn;
Fi = ctrlStruct.Fi;
Gi = ctrlStruct.Gi;

verOptions.verbose=0;
if ~isfield(sysStruct,'verified')
    sysStruct = mpt_verifySysStruct(sysStruct, verOptions);
end

if mpt_isnoise(sysStruct.noise)
    error('Cannot compute PWP Lyapunov function for systems with additive disturbances.');
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
    Acell=sysStruct.A;
    Bcell=sysStruct.B;
    fcell=sysStruct.f;
    pwasystem=1;
    if length(Acell)~=length(Pn)
        disp('In PWA case, each dynamics has to be associated with exactly one region!');
        disp('Linking dynamics to regions...');
        noU=size(Bcell{1},2);
        cc=0;
        ABFcell=cell(1,length(Pn));
        BGcell=ABFcell;
        Acell=ABFcell;
        Bcell=ABFcell;
        nu=noU;
        for ii=1:length(Pn)
            [x,R] = chebyball(Pn(ii));            % compute center of the chebyshev's ball
            Acell{ii}=[];
            for jj=1:length(sysStruct.A)          % go through all dynamics description
                if max(sysStruct.guardX{jj}*x+sysStruct.guardU{jj}*((Fi{ii}(1:noU,:)+FBgain(1:noU,:))*x+Gi{ii}(1:noU,:))-sysStruct.guardC{jj})<Options.abs_tol,    % check which dynamics is active in the region
                    Acell{ii} = sysStruct.A{jj};
                    Bcell{ii} = sysStruct.B{jj};
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

noU=size(Bcell{1},2);

try
    yalmip('clear') ;       %initialize yalmip
catch
    error('mpt_getPWPLyapFct: You need to download and install Yalmip for this function to work');
end


%load parameters
nx = dimension(Pn(1));      %no of states
nu= size(Bcell{1},2);  %no of inputs
binaryOne=dec2bin(1);   %binary one: initialize for combinatorial search later
lookahead=1;            %only examine reachability for 1 step ahead


%---------------------------------------------------
isinOpt = Options;
isinOpt.fastbreak = 1;   % to allow quick escape from isinside as soon as one mathicng region is found

x=sdpvar(nx,1);                  %define vector x
xmon=monolist(x,ndeg);          %define vector of monomials of x
lyapmatsize=length(xmon);       


%initialize problem
F=set([]);                  %constraints
alpha=sdpvar(1,1);          %lower bound V > alpha*x'*x
rho=sdpvar(1,1);            %decay rate  V-V+ > rho*x'x
parametric_variables=[alpha ; rho]; 


for i=1:length(Pn)
    Q{i}=sdpvar(1,lyapmatsize);; %initialize PWP Lyapunov variables
    
    if(isinside(Pn(i),zeros(nx,1),isinOpt))
        containsOrigin(i)=1;%The Lyapunov function is quadratic around the origin...
        %...therefore the other elements are zero    
        parametric_variables=[parametric_variables ; Q{i}(nx+2:end)']; 
    else
        containsOrigin(i)=0;
        parametric_variables=[parametric_variables ; Q{i}(:)]; 
    end
end


F = F + set(alpha>epsilon);
F = F + set(rho>epsilon);    %set positivity constraints (stronger to prevent numerical problems)
%---------------------------------------------------


%---------------------------------------------------
%Compute bounding boxes for all regions 
%this serves to speed up the subsequent computation of transition sets
if ~Options.useTmap,
    bbOptions = Options;
    bbOptions.noPolyOutput = 1;
    % we don't need the bounding box as a polytope object, just it's extreme points
    
    for i=1:length(Pn)
        %compute bounding box for region
        [R,BoxMin{i},BoxMax{i}]=bounding_box(Pn(i),bbOptions);       %get the two most extreme points
        %now extract all other extreme points of the bounding box
        for j=1:2^nx
            index=dec2bin(j-1,nx);
            for k=1:nx
                if(index(k)==binaryOne)
                    boxPoint{i}{j}(k)=BoxMax{i}(k);
                else
                    boxPoint{i}{j}(k)=BoxMin{i}(k);
                end
            end%for k=1:nx
            if(size(boxPoint{i}{j},1)<size(boxPoint{i}{j},2))
                boxPoint{i}{j}= boxPoint{i}{j}';
            end
        end%for j=1:2^(nx)
    end%Pn
end
%---------------------------------------------------

mldivideOpt = Options;
mldivideOpt.simplecheck = 1;

%---------------------------------------------------
% Find PWP Lyapunov function such that the Lyapunov values decrease for each transition

transCtr=0; %counter for the number of feasible transitions
if(pwasystem)
    unc_loop=1;
else
    unc_loop=length(Acell);
end


lenP=length(Pn);
isfulldimP = zeros(1, lenP);
for i = 1:lenP
    isfulldimP(i) = isfulldim(Pn(i));
end

if ~Options.useTmap,
    disp(['Performing Reachability Analysis   0/' num2str(lenP)])
end
for dyn_ctr=1:unc_loop
    if(~pwasystem)
        %set current system matrices
        sys=dyn_ctr;
        A=Acell{sys};
        B=Bcell{sys};
    end
    
    if Options.useTmap,
        % compute transition map
        
        % prepare dynamics cells
        Acell = cell(1, lenP);
        Fcell = cell(1, lenP);
        if pwasystem,
            Acell = ABFcell;
            Fcell = BGcell;
        else
            for ii = 1:lenP,
                Acell{ii} = A+B*(Fi{ii}(1:nu,:) + FBgain(1:nu,:));
                Fcell{ii} = B*Gi{ii}(1:nu,:);
            end
        end
            
        % compute the transition map
        if Options.verbose > -1,
            fprintf('Computing transition map...\n');
        end
        tmapOptions.maxsph = ceil(lenP^2/sphratio);
        tmap = mpt_transmap(Pn, Acell, Fcell, tmapOptions);
        
        if Options.verbose > 1,
            fprintf('Transition map discarded %.2f%% of possible transitions.\n', 100*(1 - nnz(tmap)/numel(tmap)));
        end
    end
    
    for i=1:lenP
        if(mod(i,round(lenP/5))==0) & ~Options.useTmap
            disp(['Performing Reachability Analysis   ' num2str(i) '/' num2str(length(Pn))])
        end
        %dynamics for this start region
        if pwasystem
            ABF = ABFcell{i};
            BG = BGcell{i};
            sys = i;
        else
            ABF=A+B*(Fi{i}(1:nu,:) + FBgain(1:nu,:));
            BG =B*Gi{i}(1:nu,:);
        end
        
        if(~Options.is_invariant)
            %initialize invariance check polytope
            invP=polytope;
        end
        
        if ~Options.useTmap,
            %check if polytope is mapped onto a full dimensional polytope or onto
            %a facet (this has an impact on the subsequent reachability analysis)
            if(all(abs(eig(ABF))>Options.abs_tol))
                fullMap=1;
            else
                fullMap=0;
            end
            
            %compute bounding box for region at t+1
            upperCur=ones(nx,1)*-Inf;
            lowerCur=ones(nx,1)*Inf;
            for j=1:2^nx
                boxPoint_t1{i}{j}=ABF*boxPoint{i}{j}+BG;
                lowerCur=min([lowerCur boxPoint_t1{i}{j}]')';
                upperCur=max([upperCur boxPoint_t1{i}{j}]')';
            end
        end
        
        for j=1:lenP
            if Options.useTmap,
                possible_transition = tmap(i, j);
            else
                possible_transition = 1;
                if(fullMap)
                    if(~isfulldimP(j))
                        % no transition
                        continue;
                    else
                        for k=1:nx
                            %first check if bounding boxes intersect
                            if(upperCur(k)<BoxMin{j}(k))
                                % bounding boxes do not intersect => no transition
                                % exists
                                possible_transition = 0;
                                break
                            elseif (lowerCur(k)>BoxMax{j}(k))
                                % bounding boxes do not intersect => no transition
                                % exists
                                possible_transition = 0;
                                break
                            end
                        end %nx
                    end
                end
            end
            if(possible_transition)
                %possible target, i.e. bounding boxes intersect
                
                %extract subset from region i which enters region j in one time step.
                [tmpP,keptrows,feasible]=domain(Pn(j),ABF,BG,Pn(i),lookahead);
                
                if(feasible)  
                    %transition exists from 
                    start=i;    %subset of region i
                    target=j; %to region j
                    transCtr=transCtr+1; 
                    
                    if(~Options.is_invariant)
                        %store transition for verification later
                        invP=[invP tmpP];
                    end
                    
                    %store values for later verification
                    veriStore{transCtr}.start=start;
                    veriStore{transCtr}.target=target;
                    veriStore{transCtr}.sys=sys;
                    
                    % delta_P <0
                    m = nconstr(tmpP); %facets
                    [H ,K] = double(tmpP);             %Find Hx<=K representation of transition region 
                    
                    G=K-(H*x);
                    veriStore{transCtr}.G=G;
                    
                    xnextmon=monolist((ABF*x)+BG,ndeg);     %monomials after affine translation
                    
                    
                    if (containsOrigin(start)) 
                        %Lyap function has no linear or constant elements, if the associated region contains the origin
                        deltaV{transCtr}=(Q{start}(nx+2:end)*xmon(nx+2:end));   %   V(x(k))-V(x(k+1))
                    else
                        deltaV{transCtr}=(Q{start}*xmon);                       %   V(x(k))-V(x(k+1))
                    end
                    if (containsOrigin(target))
                        %Lyap function has no linear or constant elements, if the associated region contains the origin
                        deltaV{transCtr}=deltaV{transCtr} - (Q{target}(nx+2:end)*xnextmon(nx+2:end));          %   V(x(k))-V(x(k+1))
                    else
                        deltaV{transCtr}=deltaV{transCtr} - (Q{target}*xnextmon);          %   V(x(k))-V(x(k+1))
                    end
                    
                    %add constraint on decay rate
                    deltaV{transCtr}=deltaV{transCtr} - rho*(x'*x);
                    
                    if(containsOrigin(start)==0 | containsOrigin(target)==0)
                        %add S-procedure term
                        
                        if Options.Sproc==1
                            N{transCtr} = sdpvar(m,m);     
                            N{transCtr} = N{transCtr} - diag(diag(N{transCtr}));  
                            ix = find(triu(ones(m),1));    %get indices for the elements in the upper block diagonal
                            F = F + set(N{transCtr}(ix)>0);          % Non negative multipliers
                            parametric_variables=[parametric_variables ; N{transCtr}(ix)];
                            deltaV{transCtr}=deltaV{transCtr} - (G' * N{transCtr} *G);      %S procedure to ensure negative decay rate over polytope 
                        else 
                            %% S-procedure term consists is SOS
                            sumdecay=sdpvar(1,1);
                            parametric_variables=[parametric_variables ; sumdecay];
                            
                            smxmon=monolist(x,floor(ndeg/2));%use symmetric structure for efficient implementation
                            smlyap=length(smxmon);    
                            
                            temp = [];
                            for k=1:m
                                alphadecay{transCtr,k}=sdpvar(smlyap,smlyap);
                                alphadecayx{transCtr,k}=smxmon'*alphadecay{transCtr,k}*smxmon;
                                
                                % old
                                % sumdecay=sumdecay+alphadecayx{transCtr,k}*G(k);
                                % new
                                temp = [temp alphadecayx{transCtr,k}];
                                
                                parametric_variables=[parametric_variables; alphadecay{transCtr,k}(:)];
                                F= F + set(alphadecay{transCtr,k}>0);           
                            end;
                            % new
                            sumdecay = sumdecay + temp*G;
                            deltaV{transCtr}=deltaV{transCtr}  - sumdecay;
                        end;
                        
                    end
                    F= F + set(sos(deltaV{transCtr}));          %Decay rate negative through SOS constraint
                end%feasible
            end%strcmp
        end%for j
        
        if(~Options.is_invariant)
            %check that no state exits the feasible state space
            Paux = mldivide(Pn(i), invP, mldivideOpt);
            [xcenter,R]=chebyball(Paux);

            if(max(R)<Options.abs_tol*1e3);
                %everything is great
            else
                fprintf('\n\n')
                disp('mpt_getPWPLyapFct: Partition is not invariant and therefore it cannot be asymptotically stable !!')
                dQ=[];
                dL=[];
                dC=[];
                feasible=0;
                drho=[];
                runtime=[];
                return
            end
        end
        
    end%for i
end%for sys

if(transCtr==0)
    error('mpt_getPWPLyapFct: No transition between regions detected... aborting')
else
    disp(['mpt_getPWPLyapFct: Found ' num2str(transCtr) ' feasible transitions.'])
end
%---------------------------------------------------


%---------------------------------------------------

% Make sure the Lyapunov function is positive everywhere
for start=1:length(Pn)
    
    if(containsOrigin(start))
        V1{start}=(Q{start}(nx+2:end)*xmon(nx+2:end)) - (alpha*(x'*x)); 
        F= F + set(sos(V1{start}));  
        
    elseif(Options.enforcePositivity)
        V1{start}=(Q{start}*xmon) - (alpha*(x'*x)) ; 
        
        m=nconstr(Pn(start));    
        [H , K]=double(Pn(start));          %Find Hx<=K representation of region
        G=K-(H*x);
        
        if Options.Sproc==1 
            %use S-procedure with fixed variables
            M1{start} = sdpvar(m,m); %free variable 
            M1{start} = M1{start} - diag(diag(M1{start}));  
            ix = find(triu(ones(m),1));    %get indices for the elements in the upper block diagonal
            F = F + set(M1{start}(ix)>0);
            parametric_variables=[parametric_variables ; M1{start}(ix)];
            V1{start}=V1{start} - (G' * M1{start} *G);       % S-procedure to ensure positivity on polytope
        else 
            %use S-procedure with SOS variables
            sumpos=sdpvar(1,1);
            parametric_variables=[parametric_variables ; sumpos];
            
            smxmon=monolist(x,floor(ndeg/2));%use symmetric structure for efficient implementation
            smlyap=length(smxmon);    
            
            for k=1:m
                alphapos{start,k}=sdpvar(smlyap,smlyap);
                alphaposx{start,k}=smxmon'*alphapos{start,k}*smxmon;
                sumpos=sumpos+alphaposx{start,k}*G(k);
                
                parametric_variables=[parametric_variables ; alphapos{start,k}(:)];
                F= F + set(alphapos{start,k}>0);        %Higher order multiplier , currently disabled
            end;
            
            V1{start}=V1{start} - sumpos;
        end;
        F= F + set(sos(V1{start}));  
    end 
    
    
end


runtime.setup_time=cputime-starttime;  %measure runtime
%---------------------------------------------------


%---------------------------------------------------
%Solve LMI to find Lyapunov function

starttime=cputime;             %measure runtime
disp('Solving SOS Problem')
if ~isempty(mptOptions.sdpsettings),
    options=mptOptions.sdpsettings;
else
    options=sdpsettings('Verbose',0);
end
options.sos.congruence=0;
options.sos.newton=0;
options.sos.inconsistent=0;
options.sos.traceobj=0;
options.sos.clean=0;
options.sos.extlp = 0;
options.sedumi.eps=1e-12;

if(Options.nicescale)
    F = F + set(-100 < parametric_variables < 100); %bound variables 
end

%solution=solvesos(F , [] , options , parametric_variables);
[solution,outMonomials,Grammian,primalSlacks] = solvesos(F , [] , options , parametric_variables);

%analyze solver output
if(solution.problem==4) 
    infeasible=1;   %problem occured
    disp('mpt_getPWPLyapFct: SOS solver has numerical problems... analyzing solution...');
elseif(solution.problem>0) 
    infeasible=2;   %problem occured
    disp('mpt_getPWPLyapFct: SOS solver thinks problem is infeasible.');
elseif(solution.problem<0)
    error('mpt_getPWPLyapFct: There is a problem with your SOS solver ! Go to the SeDuMi directory and type ''make''.')
else
    infeasible=0;
    disp('mpt_getPWPLyapFct: SOS solver found feasible solution.')
end


%CHECKING IF SOS DECOMPOSITION IS CERTIFIABLE...
mineig = zeros(length(Grammian),1);
for i = 1:length(Grammian)
    mineig(i,1) = min(eig(Grammian{i}))/length(Grammian{i});
end             
if all(mineig>primalSlacks) & all(checkset(F(find(~is(F,'sos'))))>=0)
    %disp('CERTIFIABLE SOS DECOMPOSITION : NO POSTPROC NEEDED!');
    verifySolution=0;
    infeasible=0;
else
    verifySolution=1;
end


%---------------------------------------------------
% Extract data from sdpsolver
for i=1:length(Q)
    %extract Lyapunov Function
    dQ{i}=double(Q{i});
    if(containsOrigin(i))
        dQ{i}(1:nx+1)=0;
    end
end;
drho=double(rho);

%---------------------------------------------------



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(~Options.enforcePositivity)
    fprintf('mpt_getPWPLyap: Computing lower bound on Lyapunov function...');
    %compute lower bound now...
    
    alpha=sdpvar(1,1);                %lower bound 
    parametric_variables=[alpha];
    Fposver=set(alpha>0);            %initialize lower bound constraint
    
    for start=1:length(Pn)            %check if Lyap Function SOS over polytope
        V1{start}=(dQ{start}*xmon) - (alpha*(x'*x)); 
        
        if(~containsOrigin(start))
            %add S-procedure terms
            m=nconstr(Pn(start));    
            [H , K]=double(Pn(start));          %Find Hx<=K representation of region
            G=K-(H*x);
            if Options.Sproc==1 
                %use S-procedure with fixed variables
                M1{start} = sdpvar(m,m); %free variable 
                M1{start} = M1{start} - diag(diag(M1{start}));
                ix = find(triu(ones(m),1));    %get indices for the elements in the upper block diagonal
                
                Fposver = Fposver + set(M1{start}(ix)>0);
                parametric_variables=[parametric_variables ; M1{start}(ix)];
                V1{start}=V1{start} - (G' * M1{start} *G);       % S-procedure to ensure positivity on polytope
            else 
                %use S-procedure with SOS variables
                sumpos=sdpvar(1,1);
                parametric_variables=[parametric_variables ; sumpos];
                
                smxmon=monolist(x,floor(ndeg/2));%use symmetric structure for efficient implementation
                smlyap=length(smxmon);    
                
                for k=1:m
                    alphapos{start,k}=sdpvar(smlyap,smlyap);
                    alphaposx{start,k}=smxmon'*alphapos{start,k}*smxmon;
                    sumpos=sumpos+alphaposx{start,k}*G(k);
                    
                    parametric_variables=[parametric_variables ; alphapos{start,k}(:)];
                    Fposver= Fposver + set(alphapos{start,k}>0);        %Higher order multiplier , currently disabled
                end;
                V1{start}=V1{start} - sumpos;
            end;
        end
        
        Fposver= Fposver + set(sos(V1{start}));  
    end
    
    %setup lower bound computation and solve sos 
    Fposver=Fposver+ set(alpha<1);            %bound objective
    
    if(Options.nicescale)
        Fposver = Fposver + set(-100 < parametric_variables < 100); %bound variables 
    end
    
    %solution=solvesos(Fposver,-alpha,options,parametric_variables);
    [solution,outMonomials,Grammian,primalSlacks] = solvesos(Fposver,-alpha,options,parametric_variables);
    
    
    %analyze solver output
	if(solution.problem==4) 
        infeasible=1;   %problem occured
        fprintf('      numerical problems... analyzing solution... \n');
	elseif(solution.problem>0) 
        infeasible=2;   %problem occured
        fprintf('       SOS solver thinks problem is infeasible.\n');
	elseif(solution.problem<0)
        error('mpt_getPWPLyapFct: There is a problem with your SOS solver ! Go to the SeDuMi directory and type ''make''.')
	else
        infeasible=0;
        fprintf('           success!\n')
	end
end

dalpha = double(alpha); %extract lower bound parameter

runtime.solver_time=cputime-starttime;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(~verifySolution)
    %CHECKING IF SOS DECOMPOSITION IS CERTIFIABLE...
    %If previous certification failed, then manual verification will be run
    %regardless of the results here. Hence, the "if"-condition.

	mineig = zeros(length(Grammian),1);
	for i = 1:length(Grammian)
        mineig(i,1) = min(eig(Grammian{i}))/length(Grammian{i});
	end             
	if all(mineig>primalSlacks) & all(checkset(Fposver(find(~is(Fposver,'sos'))))>=0)
        %disp('CERTIFIABLE SOS DECOMPOSITION : NO POSTPROC NEEDED!');
        verifySolution=0;
        infeasible=0;
	else
        verifySolution=1;
	end
end


starttime=cputime;
if(verifySolution | infeasible==1)
    disp('mpt_getPWPLyapFct: SOS decomposition not certifiable. Running manual verification of results...');
    
    
    infeasible=0; %initialize
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fprintf('mpt_getPWPLyapFct: Verifying Decay rate...');
    
    tdec=sdpvar(1,1);
    Fdecayver=set(tdec>0);          %initialize decay constraints
    parametric_variables=[tdec];
    
    for i=1:transCtr            %check if decay rate is SOS over polytope 
        start=veriStore{i}.start;
        target=veriStore{i}.target;
        
        if pwasystem
            ABF = ABFcell{start};
            BG = BGcell{start};
        else
            ABF=A+B*(Fi{start}(1:nu,:) + FBgain(1:nu,:));
            BG =B*Gi{start}(1:nu,:);
        end
        
        G=veriStore{i}.G;
        
        xnextmon=monolist((ABF*x)+BG,ndeg);     %monomials after affine translation
        %   V(x(k+1))-V(x(k))
        if (containsOrigin(start) & containsOrigin(target)) 
            deltaVver{i}= (dQ{start}*xmon) - (dQ{target}*xnextmon);         
        elseif Options.Sproc==1
            %is N{transCtr} nonnegative
            %if not set to brute force zero
            dN=double(N{i});
            dN(find(dN<epsilon))=0;
            deltaVver{i}= (dQ{start}*xmon) - (dQ{target}*xnextmon) - (G' * dN *G);         
            
        else     
            sumdecay=sdpvar(1,1);
            parametric_variables=[parametric_variables ; sumdecay];
            
            deltaVver{i}= (dQ{start}*xmon) - (dQ{target}*xnextmon);   
            
            smxmon=monolist(x,floor(ndeg/2));%use symmetric structure for efficient implementation
            smlyap=length(smxmon);    
            %computing S-procedure terms...
            for k=1:m
                dAlpha  =   double(alphadecay{i,k});
                if(min(eig(dAlpha))<0)
                    disp(['mpt_getPWPLyapFct: Higher order S-procedure is not SOS!']);
                    infeasible=infeasible+1;
                    break
                end
                
                decayVerX{i,k}=smxmon'*dAlpha*smxmon;
                sumdecay=sumdecay+decayVerX{i,k}*G(k);
            end;
            deltaVver{i}=deltaVver{i}  - sumdecay;  
        end
        
        Fdecayver = Fdecayver + set(sos(deltaVver{i}-tdec*x'*x));
    end;    
    
    %setup decay versification and solve sos 
    Fdecayver = Fdecayver + set(tdec<1); %bound objective
    solution=solvesos(Fdecayver,-tdec,options,parametric_variables);
    
    if(solution.problem==4) 
        disp(['mpt_getPWPLyapFct: Decay Verification: Numerical Problems']);
        infeasible=infeasible+1;   %problem occured
    elseif(solution.problem>0)   
        infeasible=infeasible+1;   %problem occured
        disp(['mpt_getPWPLyapFct:Decay Verification: SOS solver has found infeasible transitions']);
    elseif(solution.problem<0)
        error('mpt_getPWPLyapFct: There is a problem with your SOS solver ! Go to the SeDuMi directory and type ''make''.')
    else
        fprintf('       success! \n');
    end
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    fprintf('mpt_getPWPLyapFct: Verifying Positivity...');
    
    
    %verify positivity results
    tpos=sdpvar(1,1);
    parametric_variables=[tpos];
    Fposver=set(tpos>epsilon);        %initialize positivity constraints
    smxmon=monolist(x,floor(ndeg/2));%use symmetric structure for efficient implementation
    smlyap=length(smxmon);    
            
    for start=1:length(Pn)            %check if decay rate is SOS over polytope
        Vver{start}= (dQ{start}*xmon);
        
        if(containsOrigin(start))
            Vver{start}= Vver{start} -tpos*x'*x;          
        elseif Options.Sproc==1
            %is M{start} nonnegative
            %if not set to brute force zero
            dM1=double(M1{start});
            dM1(find(dM1<epsilon))=0;
            
            m=nconstr(Pn(start));    
            [H , K]=double(Pn(start));          %Find Hx<=K representation of region
            G=K-(H*x);
            
            Vver{start}= Vver{start} - (G' * dM1 *G) -tpos*x'*x;          
        else
            sumpos=sdpvar(1,1);
            parametric_variables=[parametric_variables ; sumpos];
            
            %computing S-procedure terms...
            for k=1:m
                dAlpha  =   double(alphapos{start,k});
                if(min(eig(dAlpha))<0)
                    disp(['mpt_getPWPLyapFct: Higher order S-procedure is not SOS!']);
                    infeasible=infeasible+1;
                    break
                end
                
                posVerX{i,k}=smxmon'*dAlpha*smxmon;
                sumpos=sumpos+posVerX{i,k}*G(k);
            end;
             Vver{start} =  Vver{start}  - sumpos;  
        end
        
        Fposver=Fposver + set(sos(Vver{start}));
    end
    
    %setup positivity verification solve sos 
    Fposver=Fposver + set(tpos<1);  %bound objective
    solution=solvesos(Fposver,-tpos,options,parametric_variables);
    
    if(solution.problem==4) 
        disp(['mpt_getPWPLyapFct: Positivity Verification: Numerical Problems']);
    elseif(solution.problem>0) 
        infeasible=infeasible+1;   %problem occured
        disp(['mpt_getPWPLyapFct: Positivity Verification failed']);
    elseif(solution.problem<0)
        error('mpt_getPWPLyapFct: There is a problem with your SOS solver ! Go to the SeDuMi directory and type ''make''.')
    else
        fprintf('       success! \n');
    end
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
end;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%if you want to be really safe...
if(Options.debug_level>1 & ~infeasible)
    disp('mpt_getPWPLyapFct: Verification via gridding...')
    Options.gridpoints=15;
    [islyapfun] = mpt_checkLyapFct(ctrlStruct,dQ,Options);
    infeasible=~islyapfun;
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
runtime.postproc_time=cputime-starttime;



fprintf('\n\n');

if(dalpha>=epsilon & drho>=epsilon & ~infeasible)
    disp(['mpt_getPWPLyapFct: SUCCESS:  Found Piecewise Polynomial Lyapunov function of degree ' num2str(ndeg)])
    feasible=1;
else
    disp(['mpt_getPWPLyapFct: FAILURE:  Failed to find Piecewise Polynomial Lyapunov function of degree ' num2str(ndeg)])
    feasible=0;
end

return