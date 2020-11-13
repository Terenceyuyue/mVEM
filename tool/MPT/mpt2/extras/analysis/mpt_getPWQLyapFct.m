function [dQ,dL,dC,feasible,drho,runtime]=mpt_getPWQLyapFct(ctrl,Options)
%MPT_GETPWQLYAPFCT Calculates PWQ Lyapunov function
%
% [dQ,dL,dC,feasible,drho]=mpt_getPWQLyapFct(ctrlStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function attempts to compute a piecewise quadratic Lyapunov function
% PWQ(x) which guarantees exponential stability. The following is satisfied
% PWQ(x) = x'dQ{r}x + x'dL{r}+ dC{r}
% PWQ(x(k+1)) - PWQ(x(k))<= rho * x^2      
% (rho must be negative to guarantee exponential stability)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------  
% ctrl                        - Explicit controller (EXPCTRL object)
% Options.enforcePositivity   - If set to zero, positivity constraints for PWQ
%                               function are not included in LMI (reduces
%                               constraints and computation time) Post
%                               computation verification is performed to check
%                               if postivity holds. If not, positivity
%                               constraints are added and solution recomputed. 
% Options.abs_tol      - Absolute tolerance
% Options.lpsolver     - Which LP solver to use (see help mpt_solveLP for
%                        details) 
% Options.epsilon      - This is a tolerance factor which is introduced to
%                        turn LMI inequalities into strict inequalities.
% Options.debug_level  - If this is set to 1, the solution provided by the LMI
%                        solver will be double-checked manually. We strongly
%                        advise to set this to 1, since we've experienced
%                        numerous numerical issues with certain LMI solvers.
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
% dQ,dL,dC   - PWQ Lyapunov function:
%                PWQ(x)=x'dQ{r}x+x'*dL{r}+C{r}
%                iff x is in region r, i.e. Pn{r}.H*x<=Pn{r}.K
%    
% feasible   - 1: stable or 0: no statement about stability possible
% drho       - the maximum Lyapunov decay rate over the partition
%              (is this is greater than zero, stability cannot be guaranteed)
%              The Lyapunov value decrease Delta V <= drho * ||x||^2
%
% see also MPT_GETPWALYAPFCT, MPT_GETPWPLYAPFCT, MPT_GETCOMMONLYAPFCT, MPT_GETSTABFEEDBACK
%

% ---------------------------------------------------------------------------
%   LITERATURE:
% ---------------------------------------------------------------------------
% European Control Conference, Cambridge, UK, 2003
% "Stability & Feasibility of Constrained Receding Horizon Control"
% Grieder P., M. Luethi, P. Parrilo and M. Morari
%
% Automatica, Volume 38, issue 12, pp. 2139 - 2146 
% "Analysis of discrete-time piecewise affine and hybrid systems",
% Ferrari-Trecate G., F.A. Cuzzola, D. Mignone and M. Morari

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Pratik Biswas       
%          pbiswas@stanford.edu
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich       
%          grieder@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich
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
if ~isfield(Options,'enforcePositivity'),
    Options.enforcePositivity=0;
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

if ~isfield(Options, 'statusbar')
    Options.statusbar = 0;
end
if ~isfield(Options, 'status_min')
    Options.status_min = 0;
end
if ~isfield(Options, 'status_max')
    Options.status_max = 1;
end
if ~isfield(Options, 'closestatbar'),
    Options.closestatbar = 1;
end
statusbar = Options.statusbar;
closestatbar = Options.closestatbar;
Options.closestatbar = 0;
Options.statusbar = 0;
if statusbar,
    Options.verbose = -1;
end

dispwarning = 1;
if isa(ctrl, 'mptctrl')
    dispwarning = ~isinvariant(ctrl);
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
    error('Cannot compute PWQ Lyapunov function for systems with additive disturbances.');
end    

if isfield(ctrlStruct.probStruct, 'FBgain'),
    % handle pre-stabilization with feedback
    FBgain = ctrlStruct.probStruct.FBgain;
else
    FBgain = zeros(size(Fi{1}));
end

if (iscell(sysStruct.A))
    if(isfield(sysStruct,'Aunc') & ~isempty(sysStruct.Aunc))
        error('Cannot handLe PWA systems with polytopic uncertainty')
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
                    ABFcell{ii} = sysStruct.A{jj}+sysStruct.B{jj}*(Fi{ii}(1:nu,:) + FBgain(1:nu, :));
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
    error('mpt_getPWQLyapFct: You need to download and install Yalmip for this function to work');
end

progress = 0;
if statusbar,
    if ~isfield(Options, 'status_handle')
        Options.status_handle = mpt_statusbar('Computing invariant set...');
        %closestatbar = 1;
    end
end

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, progress, Options.status_min, Options.status_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

%load parameters
n = dimension(Pn(1));      %no of states
nu= size(Bcell{1},2);  %no of inputs
binaryOne=dec2bin(1);   %binary one: initialize for combinatorial search later
lookahead=1;            %only examine reachability for 1 step ahead


%---------------------------------------------------
%Initialize LMI Variables
%PWQ Lyapunov function: x'Q{i}x+x'L{i}+C{i}      if Hn{i}*x<=Kn{i}
isinOpt = Options;
isinOpt.fastbreak = 1;   % to allow quick escape from isinside as soon as one mathicng region is found
temp = repmat(1,1,length(Pn));
Q = sdpvar(n*temp,n*temp,'symmetric'); %initialize PWQ Lyapunov variables
L = sdpvar(n*temp,temp,'full');        %initialize PWQ Lyapunov variables
C = sdpvar(temp,temp,'full');          %initialize PWQ Lyapunov variables
if ~iscell(L),
    L = {L};
    C = {C};
    Q = {Q};
end
for i=1:length(Pn)
    if(isinside(Pn(i),zeros(n,1),isinOpt))
        containsOrigin(i)=1;%The Lyapunov function is quadratic around the origin...
        %...therefore the linear elements L{:}=0 and C{:}=0;
        L{i} = zeros(n,1);  
        C{i} = 0;  
    else
        containsOrigin(i)=0;
    end
end
rho=sdpvar(1,1);   %optimization variable for exponential stability
%---------------------------------------------------


%---------------------------------------------------
%Compute bounding boxes for all regions 
%this serves to speed up the subsequent computation of transition sets

bbOptions = Options;
bbOptions.noPolyOutput = 1;
% we don't need the bounding box as a polytope object, just it's extreme points

if ~Options.useTmap,
    for i=1:length(Pn)
        %compute bounding box for region
        [R,BoxMin{i},BoxMax{i}]=bounding_box(Pn(i),bbOptions);       %get the two most extreme points
        %now extract all other extreme points of the bounding box
        for j=1:2^n
            index=dec2bin(j-1,n);
            for k=1:n
                if(index(k)==binaryOne)
                    boxPoint{i}{j}(k)=BoxMax{i}(k);
                else
                    boxPoint{i}{j}(k)=BoxMin{i}(k);
                end
            end%for k=1:n
            if(size(boxPoint{i}{j},1)<size(boxPoint{i}{j},2))
                boxPoint{i}{j}= boxPoint{i}{j}';
            end
        end%for j=1:2^(n)
    end%Pn
end
%---------------------------------------------------


mldivideOpt = Options;
mldivideOpt.simplecheck = 1;

%----------------------------------------------------------------------------------------------------
%ADD CONSTRAINTS FOR LYAPUNOV DECAY

% Find PWQ Lyapunov function such that the Lyapunov values decrease for each transition

transCtr=0; %counter for the number of feasible transitions
myprog=set([]); %initialize LMI constraint to null
if(pwasystem)
    unc_loop=1;
else
    unc_loop=length(Acell);
end

isfulldimP = zeros(1, length(Pn));
lenP=length(Pn);
for i = 1:lenP
    isfulldimP(i) = isfulldim(Pn(i));
end

if Options.verbose > -1 & ~Options.useTmap,
    disp(['Performing Reachability Analysis   0/' num2str(lenP)])
end

for dyn_ctr=1:unc_loop
    if(~pwasystem)
        %set current system matrices
        sys=dyn_ctr;
        A=Acell{sys};
        B=Bcell{sys};
    end
    
    if statusbar,
        prog_min = Options.status_min;
        prog_max = prog_min + (Options.status_max - Options.status_min)/unc_loop;
    end

    progress = 0;
    if statusbar,
        if isempty(mpt_statusbar(Options.status_handle, progress, prog_min, prog_max)),
            mpt_statusbar;
            error('Break...');
        end     
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
        
        progress = (i-1)/lenP;
        if statusbar,
            if isempty(mpt_statusbar(Options.status_handle, progress, prog_min, prog_max)),
                mpt_statusbar;
                error('Break...');
            end     
        end
        
        if ~Options.useTmap,
            if(mod(i,round(lenP/5))==0) 
                if Options.verbose > -1,
                    disp(['Performing Reachability Analysis   ' num2str(i) '/' num2str(length(Pn))])
                end
            end
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
            upperCur=ones(n,1)*-Inf;
            lowerCur=ones(n,1)*Inf;
            for j=1:2^n
                boxPoint_t1{i}{j}=ABF*boxPoint{i}{j}+BG;
                lowerCur=min([lowerCur boxPoint_t1{i}{j}]')';
                upperCur=max([upperCur boxPoint_t1{i}{j}]')';
            end
        end
        
        for j=1:lenP
            
            if statusbar,
                if mod(j, 5)==0,
                    if isempty(mpt_statusbar(Options.status_handle, progress, prog_min, prog_max)),
                        mpt_statusbar;
                        error('Break...');
                    end     
                end
            end
        
            if Options.useTmap,
                % use information from the transition map
                possible_transition = tmap(i, j);
            else
                possible_transition = 1;
                if(fullMap)
                    if(~isfulldimP(j))
                        % no transition
                        continue;
                    else
                        for k=1:n
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
                        end %n
                    end
                end
            end
            
            if(possible_transition)
                %possible target, i.e. bounding boxes intersect
                
                %extract subset from region i which enters region j in one time step.
                [tmpP,keptrows,feasible]=domain(Pn(j),ABF,BG,Pn(i),lookahead);
                
                if(feasible)  
                    %transition exists from 
                    start  = i;    %subset of region i
                    target = j;    %to region j
                    transCtr=transCtr+1; 
                    
                    if(~Options.is_invariant)
                        %store transition for verification later
                        invP=[invP tmpP];
                    end
                    
                    %store values for later verification
                    veriStore{transCtr}.start=start;
                    veriStore{transCtr}.target=target;
                    veriStore{transCtr}.sys=sys;
                    veriStore{transCtr}.ABF=ABF;
                    veriStore{transCtr}.BG =BG;
                    if (containsOrigin(start) & containsOrigin(target)) 
                        quadTrans(transCtr) = 1;  %transition between quadratic functions
                        veriStore{transCtr}.HK = [];
                    else
                        quadTrans(transCtr)=0;  %not transition between quadratic functions
                        veriStore{transCtr}.HK = double(-tmpP);                                                
                    end
                    
%                     % delta_P <0
%                     if (containsOrigin(start) & containsOrigin(target)) 
%                         %In region 1 the Lyapunov-function is merely quadratic and not PWQ!
%                         delta_P{transCtr} = [ ABF'*Q{target}*ABF-Q{start}-rho*eye(n)];   
%                         W{transCtr} = -delta_P{transCtr};  
%                         quadTrans(transCtr) = 1;  %transition between quadratic functions
%                     else
%                         quadTrans(transCtr) = 0;  %not transition between quadratic functions
%                         
%                         %compute matrix for decay rate in Lyapunov function
%                         temp = ABF'*Q{target}*BG+(ABF'*L{target}-L{start})/2;
%                         delta_P{transCtr} = [ ABF'*Q{target}*ABF-Q{start}-rho*eye(n)            temp ;...
%                                 temp'    BG'*Q{target}*BG+C{target}+BG'*L{target}-C{start}];   
%                         
%                         
%                         m = nconstr(tmpP); %facets
%                         HK = double(-tmpP);             %Find Hx<=K representation of transition region 
%                         veriStore{transCtr}.HK=HK;
%                         N{transCtr} = sdpvar(m,m);     
%                         N{transCtr} = N{transCtr} - diag(diag(N{transCtr}));     
%                         ix = find(triu(ones(m),1));    %get indices for the elements in the upper block diagonal
%                         myprog = myprog + lmi('N{transCtr}(ix)>0');
%                         
%                         W{transCtr} = -delta_P{transCtr} - HK' * N{transCtr} * HK;       %Polytope based decay rate constraint 
%                         
%                     end
%                    
%                    %myprog = addLmi(myprog,'W{transCtr}>0');                 %eigenvalues must be greater equal zero
%                    
%                    myprog = myprog + lmi('W{transCtr}>0');                 %eigenvalues must be greater equal zero
                    
                end%feasible
            end%possible_transition
        end%for j
        
        if(~Options.is_invariant)
            %check that no state exits the feasible state space
            Paux = mldivide(Pn(i), invP, mldivideOpt);
            [xcenter,R]=chebyball(Paux);
            
            if(max(R)<Options.abs_tol*1e3);
                %everything is great
            else
                fprintf('\n\n')
                disp('mpt_getPWQLyapFct: Partition is not invariant and therefore it cannot be asymptotically stable !!')
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
    if closestatbar,
        mpt_statusbar;
    end
    error('mpt_getPWQLyapFct: No transition between regions detected... aborting')
else
    disp(['mpt_getPWQLyapFct: Found ' num2str(transCtr) ' feasible transitions.'])
end

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, progress, prog_min, prog_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

% Create sdpvar variable N in vectorized way
for i = 1:length(veriStore)
    dims(1,i) = size(veriStore{i}.HK,1);
end
dims(dims==0) = 1; % Dummy variables
N = sdpvar(dims,dims);
if all(size(dims)==1),
    N = {N};
end

% Setup LMIs for transition deacay
for transCtr = 1:length(veriStore)
    start = veriStore{transCtr}.start;
    target = veriStore{transCtr}.target;
    HK =  veriStore{transCtr}.HK;
    ABF = veriStore{transCtr}.ABF;
    BG = veriStore{transCtr}.BG;
    if quadTrans(transCtr)
        %In region 1 the Lyapunov-function is merely quadratic and not PWQ!
        delta_P{transCtr} = [ ABF'*Q{target}*ABF-Q{start}-rho*eye(n)];
        W{transCtr} = -delta_P{transCtr};
    else
        %compute matrix for decay rate in Lyapunov function
        temp = ABF'*Q{target}*BG+(ABF'*L{target}-L{start})/2;
        delta_P{transCtr} = [ ABF'*Q{target}*ABF-Q{start}-rho*eye(n)  temp ;...
            temp'    BG'*Q{target}*BG+C{target}+BG'*L{target}-C{start}];

        HK = veriStore{transCtr}.HK;
        m = size(HK,1);
        N{transCtr} = N{transCtr} - diag(diag(N{transCtr}));
        ix = find(triu(ones(m),1));    %get indices for the elements in the upper block diagonal
        myprog = myprog + set(N{transCtr}(ix)>0);
        W{transCtr} = -delta_P{transCtr} - HK' * N{transCtr} * HK;       %Polytope based decay rate constraint
    end   
    myprog = myprog + set(W{transCtr}>0);
end


%END OF ADDING CONSTRAINTS FOR LYAPUNOV DECAY
%----------------------------------------------------------------------------------------------------


%----------------------------------------------------------------------------------------------------
% Make sure the Lyapunov function is positive everywhere

for start=1:length(Pn)
    
    if statusbar,
        if mod(start, 5)==0,
            if isempty(mpt_statusbar(Options.status_handle, progress, prog_min, prog_max)),
                mpt_statusbar;
                error('Break...');
            end     
        end
    end
    
    %PWQ_P>=0 ,      P must be positive definite to be a Lyapunov function
    if(containsOrigin(start))
        lyap{start} = Q{start} - epsilon*eye(n);
        myprog = myprog +  set(lyap{start}>0);        %eigenvalues must be greater equal zero
        
    elseif(Options.enforcePositivity)
        
        m=nconstr(Pn(start));    
        HK=double(-Pn(start));          %Find Hx<=K representation of region
        
        M{start} = sdpvar(m,m); %free variable with nonnegative elements
        M{start} = M{start} - diag(diag(M{start}));  
        ix = find(triu(ones(m),1));    %get indices for the elements in the upper block diagonal
        %myprog = addLmi(myprog,'M{start}(ix)>0');      %elements must be greater equal zero
        myprog = myprog + set(M{start}(ix)>0);      %elements must be greater equal zero
        
        PWQ_P=[Q{start}-epsilon*eye(n)     0.5*L{start} ;  0.5*L{start}'   C{start}];
        lyap{start} = PWQ_P - HK' * M{start} * HK;            %Polytope based positivity constraint
        myprog = myprog +  set(lyap{start}>0);        %eigenvalues must be greater equal zero
    else
        %just for scaling
        lyap{start} = [Q{start}-epsilon*eye(n)     0.5*L{start} ;  0.5*L{start}'   C{start}];
    end
    
    %bound the PWQ P to get "nicely" scaled functions (not too big)  
    nn = size(lyap{start},1);
    myprog = myprog + set(lyap{start}<10*eye(nn));
    myprog = myprog + set(lyap{start}>-10*eye(nn));
    
    
    if(Options.nicescale)
        %add bunch of additional constraints to get nice plots
        maxVal=1e1;
        myprog = myprog + set(trace(Q{start})>1);
        for i=1:n
            myprog = myprog + set(norm(Q{start}(:,i))<maxVal);
        end
        if(~containsOrigin(start))
            myprog = myprog + set(norm(L{start})<maxVal);
            myprog = myprog + set(C{start}<maxVal);
            myprog = myprog + set(C{start}>-maxVal);
        end
    end
end

runtime.setup_time=cputime-starttime; %measure runtime

%END OF Make sure the Lyapunov function is positive everywhere
%----------------------------------------------------------------------------------------------------


if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, progress, prog_min, prog_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

%------------------------------------------------------------------------------------------------------
%SOLVE MAIN LMI

starttime=cputime;
disp('Solving LMI Problem')
myprog = myprog + lmi('rho<-epsilon');%objective to be minimized


%options=sdpsettings('Verbose',1,'solver','sdpt3');
if ~isempty(mptOptions.sdpsettings)
	options = mptOptions.sdpsettings;
	%options.solver = 'sedumi';
else
	%options=sdpsettings('Verbose',0,'solver','sedumi');
    options=sdpsettings('Verbose',0);
end
options.sedumi.maxiter=200;                                    %better to have quality than runtime
options.sedumi.stepdif=0;
options.sedumi.eps=1e-12;
solution = solvesdp(myprog,[],options);                        %find solution using LMI solver

%analyze solver output
if(solution.problem==4) 
    res=checkset(myprog);
    if min(res)>0,
        infeasible=1;   %double check result later...
        disp(['mpt_getPWQLyapFct: Numerical problems with LMI solver but solution is feasible.']);
    else
        infeasible=1;   %double check result later...
        disp(['mpt_getPWQLyapFct: Numerical problems with LMI solver:  residual ' num2str(min(res)) ' (should be positive).'])
    end
elseif(solution.problem>0) 
    infeasible=2;   %problem occured
    res=checkset(myprog);
    disp(['mpt_getPWQLyapFct: LMI solver thinks problem is infeasible:  residual ' num2str(min(res)) ' (should be positive).'])
elseif(solution.problem<0)
    if closestatbar,
        mpt_statusbar;
    end
    error('mpt_getPWQLyapFct: There is a problem with your LMI solver ! Go to the SeDuMi directory and type ''make''.')
else
    infeasible=0;
    disp('mpt_getPWQLyapFct: LMI solver found feasible solution.')
end

%extract solution from sedumi 
for i=1:length(Q)
    dQ{i}=double(Q{i});
    if(~containsOrigin(i))
        dL{i}=double(L{i});
        dC{i}=double(C{i});
    else
        dL{i}=zeros(n,1);
        dC{i}=0;
    end
end
for i=1:transCtr
    if(~quadTrans(i))
        dN{i}=double(N{i});
        if(~isreal(dN{i})) 
            if closestatbar,
                mpt_statusbar;
            end
            error('mpt_getPWQLyapFct: Error in LMI solver, variable should be real')
        end
        dN{i}(find(dN{i}<0))=0;
    end
end

drho=double(rho);

%END OF SOLVING MAIN LMI 
%------------------------------------------------------------------------------------------------------



%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%ENFORCING POSITIVITY

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, progress, prog_min, prog_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

if(~Options.enforcePositivity)
    disp('mpt_getPWQLyapFct: Computing lower bound on Lyapunov function...');
    %compute lower bound now...
    myprog=set([]);
    
    % Check if the Lyapunov function is positive everywhere
    for start=1:length(Pn)
        
        %PWQ_P>=0 ,      P must be positive definite to be a Lyapunov function
        if(containsOrigin(start))
            %these constraints are already enforced
        else
            
            m=nconstr(Pn(start));    
            HK=double(-Pn(start));          %Find Hx<=K representation of region
            
            M{start} = sdpvar(m,m); %free variable with nonnegative elements
            M{start} = M{start} - diag(diag(M{start}));  
            ix = find(triu(ones(m),1));    %get indices for the elements in the upper block diagonal
            myprog = myprog + set(M{start}(ix)>0);      %elements must be greater equal zero
            
            PWQ_P=[dQ{start}-epsilon*eye(n)     0.5*dL{start} ;  0.5*dL{start}'   dC{start}];
            lyap{start} = PWQ_P - HK' * M{start} * HK;            %Polytope based positivity constraint
            myprog = myprog +  set(lyap{start}>0);        %eigenvalues must be greater equal zero
        end
    end
    
    if all(containsOrigin),
        fprintf('mpt_getPWQLyapFct: Skipping verification because all regions contain the origin.\n');
    else
        solution = solvesdp(myprog,[],options);                        %find solution using LMI solver
    end
    
    %analyze solver output
    if(solution.problem==4) 
        res=checkset(myprog);
        if min(res)>0,
            infeasible=1; %double check result later...
            disp(['mpt_getPWQLyapFct: Numerical problems with LMI solver but solution is feasible.']);
        else
            infeasible=1;  %double check result later...
            disp(['mpt_getPWQLyapFct: Numerical problems with LMI solver:  residual ' num2str(min(res)) ' (should be positive).'])
        end
    elseif(solution.problem>0) 
        infeasible=2;   %problem occured
        res=checkset(myprog);
        disp(['mpt_getPWQLyapFct: LMI solver thinks problem is infeasible:  residual ' num2str(min(res)) ' (should be positive).'])
    elseif(solution.problem<0)
        if closestatbar,
            mpt_statusbar;
        end
        error('mpt_getPWQLyapFct: There is a problem with your LMI solver ! Go to the SeDuMi directory and type ''make''.')
    else
        infeasible=0;
        disp('mpt_getPWQLyapFct: LMI solver found feasible solution.')
    end
end


%extract solution from sedumi 
for i=1:length(Pn)
    if(~containsOrigin(i))
        dM{i}=double(M{i});
        if(~isreal(dM{i}))
            if closestatbar,
                mpt_statusbar;
            end
            error('mpt_getPWQLyapFct: Error in LMI solver, variable should be real')
        end
        dM{i}(find(dM{i}<0))=0;
    else
        dM{i}=[];
    end
end

runtime.solver_time=cputime-starttime;
%END OF ENFORCING POSITIVITY
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%-------------------------------------------------------------------------------------------
%VERIFICATIONS

starttime=cputime;

if(drho>0)
    disp(['mpt_getPWQLyapFct: Rho is not negative!!  (rho=' num2str(drho) ')'])
    infeasible=2;   %decay rate is not negative
end


if(infeasible<2)
    disp('mpt_getPWQLyapFct: verifying positivity of function...')
    %This is a manual double-check that the results given by SeDuMi really hold.
    %We've experienced some inconsistencies and therefore added this section
    
    %Also since the positivity constraints are sometimes not enforced, this also checks for positivity of Lyapunov function
    %If not satisfied, it adds the positivity constraints and recomputes a solution
    
    %check positivity...
    for start=1:length(Pn)
        if(containsOrigin(start))
            PWQ_P   = dQ{start};        %eigenvalues must be greater equal zero
        else
            HK  =   double(-Pn(start));          %Find Hx<=K representation of region
            
            PWQ_P   =   [dQ{start}  0.5*dL{start} ;  0.5*dL{start}'   dC{start}];
            PWQ_P   =   PWQ_P - HK' * dM{start} * HK;            %Polytope based positivity constraint
        end
        
        if(min(eig(PWQ_P))<0)
            disp(['mpt_getPWQLyapFct: Lyapunov function not positive in region ' num2str(start)])
            infeasible=2;
        end
    end
    
    if(infeasible<2)
        disp('mpt_getPWQLyapFct: verifying decay rate of function...')
        %no verification problems so far...
        
        ddrho=-eps*100; %set to very small dummy value for verification
        
        for kk=1:transCtr
            start=veriStore{kk}.start;   %moving from region i to region j
            target=veriStore{kk}.target;  %moving from region i to region j
            if(isfield(veriStore{kk},'sys'))
                sys=veriStore{kk}.sys;
            end
            if pwasystem
                ABF = ABFcell{start};
                BG  = BGcell{start};
            else
                ABF = A+B*(Fi{start}(1:noU,:) + FBgain(1:noU,:));
                BG  = B*Gi{start}(1:noU,:);
            end
            
            if (containsOrigin(start) & containsOrigin(target))
                delta_P = [ ABF'*dQ{target}*ABF-dQ{start}-ddrho*eye(n)];   
                W1 = -delta_P;  
            else
                delta_P = [ ABF'*dQ{target}*ABF-dQ{start}-ddrho*eye(n)           ABF'*dQ{target}*BG+(ABF'*dL{target}-dL{start})./2 ;...
                        BG'*dQ{target}*ABF+(ABF'*dL{target}-dL{start})'./2   BG'*dQ{target}*BG+dC{target}+BG'*dL{target}-dC{start}];  
                
                HK=veriStore{kk}.HK;
                W1 = -delta_P - HK' * dN{kk} * HK;  
            end
            
            if(min(eig(W1))<0)
                disp(['mpt_getPWQLyapFct: Decay Rate: Eigenvalue smaller zero when moving from region ' num2str(start) ' to ' num2str(target) '   (min. eig.=  ' num2str(min(eig(W1))) ')'])
                infeasible=2;
            end%min eig W1
        end
    end
    
    if(infeasible==1)
        %value is still at 1. This means no problems were detected
        disp('mpt_getPWQLyapFct: Verification showed that function is Lyapunov.')
        infeasible=0;
    end
end

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, 1, Options.status_min, Options.status_max)),
        mpt_statusbar;
        error('Break...');
    end
end

%END OF VERIFICATION
%-------------------------------------------------------------------------------------------


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%if you want to be really safe...
if(Options.debug_level>1 & ~infeasible)
    disp('mpt_getPWQLyapFct: Verification via gridding...')
    Options.gridpoints=15;
    [islyapfun] = mpt_checkLyapFct(ctrlStruct,dQ,dL,dC,Options);
    infeasible=~islyapfun;
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
runtime.postproc_time=cputime-starttime;


fprintf('\n\n')
if(~infeasible & drho<=0)
    disp(['mpt_getPWQLyapFct: SUCCESS:  Found Piecewise Quadratic Lyapunov function.'])
    feasible=1;
else 
    disp(['mpt_getPWQLyapFct: FAILURE:  No Piecewise Quadratic Lyapunov function could be found.'])
    feasible=0;
end

if closestatbar,
    mpt_statusbar;
end
return
