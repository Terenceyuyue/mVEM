function [lyapunovL,lyapunovC,feasible,drho]=mpt_getPWALyapFct(ctrl,Options)
%MPT_GETPWALYAPFCT Calculates PWA Lyapunov function
%
% [lyapunovL,lyapunovC,feasible,drho]=mpt_getPWALyapFct(ctrl,Options)
%
%-----------------------------------
%   DESCRIPTION:
%-----------------------------------
% This function attempts to compute a piecewise affine Lyapunov matrix
% which guarantees stability.
%
% PWQ = L*x+ C
% Delta PWA <= rho * ||x||_1
% (rho must be negative to guarantee exponential stability)
%
%-----------------------------------
%    INPUT:
%-----------------------------------   
% ctrl                - Explicit controller (MPTCTRL object)
% Options.abs_tol     - Absolute tolerance
% Options.lpsolver    - Which LP solver to use (see help mpt_solveLP for details)
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
% Note: If "Options" is missing or some of the fields are not defined, the default
%       values from mptOptions will be used.
%
%-----------------------------------
%   OUTPUT:
%-----------------------------------
% lyapunovL,lyapunovC       - PWA Lyapunov function:
%                             PWA(x)=lyapunovL{r}*x+lyapunovC{r}
%                             iff x is in region r, i.e. Pn{r}.H*x<=Pn{r}.K
%    
% feasible   - 1: asymptotically stable 0: no statement about stability possible
% drho       - the maximum Lyapunov decay rate over the partition
%              (is this is greater than zero, stability cannot be guaranteed)
%              The Lyapunov value decrease Delta V <= drho * ||x||_1
%              Note: If drho=0 then feasible=0, since the system is not exponentially stable 
%                    (it is merely stable)
%                    

%-----------------------------------
%   LITERATURE:
%-----------------------------------
% Based on results by Mikael Johanson in his book "Stability of PWA Systems"

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich
%          grieder@control.ee.ethz.ch
%

% ---------------------------------------------------------------------------
% Legal note:
%          This library is free software; you can redistribute it and/or
%          modify it under the terms of the GNU Lesser General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This library is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          Lesser General Public License for more details.
% 
%          You should have received a copy of the GNU Lesser General Public
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
verifySolution=Options.debug_level;
if ~isfield(Options,'abs_tol')
    Options.abs_tol = mptOptions.abs_tol;
end

if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
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

if ~mpt_isValidCS(ctrlStruct)
    error('mpt_getPWALyapFct: First argument has to be a valid controller structure! See mpt_control for details.');
end
% how many regions contain the origin?
dim = dimension(ctrlStruct.Pn(1));
[isin, contain_origin] = isinside(ctrlStruct.Pn, zeros(dim, 1));
if (ctrlStruct.probStruct.norm==2) && length(contain_origin)~=2^dim
    error('There is no chance of finding a PWA Lyapunov function for a controller computed in 2-norm. Aborting...')
end
epsilon=Options.abs_tol;

sysStruct = ctrlStruct.sysStruct;
Pn = ctrlStruct.Pn;
Fi = ctrlStruct.Fi;
Gi = ctrlStruct.Gi;

verOptions.verbose=0;
if ~isfield(sysStruct,'verified')
    sysStruct = mpt_verifySysStruct(sysStruct, verOptions);
end

if mpt_isnoise(sysStruct.noise)
    error('Cannot compute PWA Lyapunov function for systems with additive disturbances.');
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

if isfield(ctrlStruct.probStruct, 'FBgain'),
    % handle pre-stabilization with feedback
    FBgain = ctrlStruct.probStruct.FBgain;
else
    FBgain = zeros(size(Fi{1}));
end

if iscell(sysStruct.A),
    nu = size(sysStruct.B{1},2);
else
    nu = size(sysStruct.B,2);
end
if iscell(sysStruct.A)
    if(isfield(sysStruct,'Aunc') & ~isempty(sysStruct.Aunc))
        error('Cannot handle PWA systems with polytopic uncertainty')
    end
    Acell=sysStruct.A;
    Bcell=sysStruct.B;
    pwasystem=1;
    if length(Acell)~=length(Pn)
        disp('In PWA case, each dynamics has to be associated with exactly one region!');
        disp('Linking dynamics to regions...');
        noU=size(Bcell{1},2);
        cc=0;
        Acell={};
        Bcell={};
        ABFcell = cell(1,length(Pn));
        BGcell = ABFcell;
        Acell = ABFcell;
        Bcell = ABFcell;
        for ii=1:length(Pn)
            [x,R] = chebyball(Pn(ii));            % compute center of the chebyshev's ball
            for jj=1:length(sysStruct.A)          % go through all dynamics description
                if max(sysStruct.guardX{jj}*x+sysStruct.guardU{jj}*((Fi{ii}(1:noU,:)+FBgain(1:noU,:))*x+Gi{ii}(1:noU,:))-sysStruct.guardC{jj})<Options.abs_tol,    % check which dynamics is active in the region
                    Acell{ii} = sysStruct.A{jj};
                    Bcell{ii} = sysStruct.B{jj};
                    ABFcell{ii} = sysStruct.A{jj}+sysStruct.B{jj}*(Fi{ii}(1:noU,:) + FBgain(1:noU,:));
                    BGcell{ii} = sysStruct.B{jj}*Gi{ii}(1:noU,:) + sysStruct.f{jj};
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
%load parameters
n = dimension(Pn(1));      %no of states
noU= size(Bcell{1},2);  %no of inputs
binaryOne=dec2bin(1);   %binary one: initialize for combinatorial search later
lookahead=1;            %only examine reachability for 1 step ahead


isinOpt = Options;
isinOpt.fastbreak = 1;   % to allow quick escape from isinside as soon as one mathicng region is found
ctr=0;
containsOrigin = zeros(length(Pn),1);
varindex = cell(length(Pn),1);
if isfield(ctrlStruct.probStruct, 'xref')
    x_origin = ctrlStruct.probStruct.xref;
else
    x_origin = zeros(n, 1);
end
if isfield(ctrlStruct.probStruct, 'uref')
    uref = ctrlStruct.probStruct.uref;
else
    uref = zeros(nu, 1);
end
for i=1:length(Pn)
    if(isinside(Pn(i), x_origin ,isinOpt))
        containsOrigin(i)=1;    %The Lyapunov function is linear around the origin...             
        varindex{i}=ctr+1:(ctr+n);
        ctr=ctr+n;
    else
        containsOrigin(i)=0;
        varindex{i}=ctr+1:(ctr+n+1); %...and affine everywhere else.
        ctr=ctr+n+1;
    end
end

totvar=ctr+1; %total number of variables

%---------------------------------------------------
%Compute bounding boxes for all regions 
%this serves to speed up the subsequent computation of transition sets
if ~Options.useTmap,
    bbOptions = Options;
    bbOptions.noPolyOutput = 1;
    BoxMin = cell(1,length(Pn));
    BoxMax = cell(1,length(Pn));
    for i=1:length(Pn)
        %compute bounding box for region
        [R,BoxMin{i},BoxMax{i},lv,uv]=bounding_box(Pn(i),bbOptions);       %get the two most extreme points
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

progress = 0;
if statusbar,
    if ~isfield(Options, 'status_handle')
        Options.status_handle = mpt_statusbar('Computing invariant set...');
        closestatbar = 1;
    end
end
if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, progress, Options.status_min, Options.status_max)),
        mpt_statusbar;
        error('Break...');
    end     
end


%---------------------------------------------------
% Find PWQ Lyapunov function such that the Lyapunov values decrease for each transition
lenP=length(Pn);
isfulldimP = zeros(1, lenP);
for i = 1:lenP
    isfulldimP(i) = isfulldim(Pn(i));
end

Aconstr=[];
constr_origin=[];
transCtr=0; %counter for the number of feasible transitions
if(pwasystem)
    unc_loop=1;
else
    unc_loop=length(Acell);
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
        
        if(mod(i,round(lenP/10))==0) & ~Options.useTmap,
            if Options.verbose>-1,
                disp(['Building constraint matrices... ' num2str(i) '/' num2str(lenP)])
            end
        end
        if(pwasystem)
            %set current system matrices
            sys=i;
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

        if ~Options.useTmap,
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
            %check all possible target regions
            
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
                    start=i;    %subset of region i
                    target=j; %to region j
                    transCtr=transCtr+1; 
                    [vertices] = extreme(tmpP);
                    
                    % delta_P <0
                    for vctr=1:size(vertices,1)
                        Atmp=zeros(1,totvar);
                        %V(x(k+1))-V(x(k))-rho*??????x(k)??????_1<=0
                        Atmp(totvar)=-norm(vertices(vctr,:),1);       %-rho*??????x(k)??????_1
                        if(~containsOrigin(i))
                            Atmp(varindex{i})=[-vertices(vctr,:) -1];%-V(x(k))
                        else
                            Atmp(varindex{i})=[-vertices(vctr,:)];%-V(x(k))
                        end
                        %%Atmp(((j-1)*(n+1)+1):(j*(n+1)))=[((A+B*Fi{i})*vertices(vctr,:)')'+Gi{i}'*B' +1];%V(x(k+1))
                        FF = Fi{i}; FF=FF(1:nu,:);
                        GG = Gi{i}; GG=GG(1:nu) + uref;
                        
                        Atmp2=zeros(1,totvar);
                        if(~containsOrigin(j))
                            %Atmp2(varindex{j})=[((A+B*FF)*vertices(vctr,:)')'+GG'*B' +1];%+V(x(k+1))
                            Atmp2(varindex{j})=[(ABF*vertices(vctr,:)')'+BG' +1];%+V(x(k+1))
                        else
                            %Atmp2(varindex{j})=[((A+B*FF)*vertices(vctr,:)')'+GG'*B'];%+V(x(k+1))
                            Atmp2(varindex{j})=[(ABF*vertices(vctr,:)')'+BG'];%+V(x(k+1))
                        end
                        if(~all(abs(Atmp)<mptOptions.abs_tol))
                            constr_origin{size(Aconstr,1)+1}=[i j];
                            Aconstr=[Aconstr;Atmp+Atmp2];
                        elseif(~all(abs(Atmp2)<mptOptions.abs_tol))
                            if closestatbar,
                                mpt_statusbar;
                            end
                            error(['origin is not an equilibrium point for region ' num2str(i)]);
                        end
                    end
                    
                end%feasible
            end%strcmp
        end%for j
    end%for i
end%for sys
Bconstr=zeros(size(Aconstr,1),1);

if(transCtr==0)
    if closestatbar,
        mpt_statusbar;
    end
    error('mpt_getPWALyapFct: No transition between regions detected... aborting')
else
    if Options.verbose>-1,
        disp(['mpt_getPWALyapFct: Found ' num2str(transCtr) ' feasible transitions.'])
    end
end

%---------------------------------------------------
%Add positivity constraints
for i=1:lenP
    [vertices] = extreme(Pn(i));
    for vctr=1:size(vertices,1)
        
        Atmp=zeros(1,totvar);
        if(~containsOrigin(i))
            Atmp(varindex{i})=[-vertices(vctr,:) -1];
        else
            Atmp(varindex{i})=[-vertices(vctr,:)];
        end
        
        if(~all(abs(Atmp) < mptOptions.abs_tol) & ...
                ~all(abs(vertices(vctr,:)) < mptOptions.abs_tol))
            constr_origin{size(Aconstr,1)+1}=[i];
            Aconstr=[Aconstr;Atmp];
            Bconstr=[Bconstr;-norm(vertices(vctr,:),1)*epsilon];  %linear lower bound on Lyap fct.
            
            Aconstr=[Aconstr;-Atmp];
            Bconstr=[Bconstr;norm(vertices(vctr,:),1)*100];  %linear upper bound on Lyap fct.
        end
    end
end



f=zeros(totvar,1);
f(end)=1;
disp('Normalizing constraint polytope ...')

%try CPLEX first because CPLEX is MUCH more reliable for large LPs than the other solvers
defaultLP = Options.lpsolver;
prefered_solvers = [2 8 15 0 9 10 11 12 13 14 15 7 4 1 3 5];
for ii = 1:length(prefered_solvers),
    if any(ismember(mptOptions.solvers.lp_all, prefered_solvers(ii)))
        lpsolver = prefered_solvers(ii);
        if lpsolver~=defaultLP,
            fprintf('Switching to solver %s (faster computation)...\n', mpt_solverInfo('lp', lpsolver));
        end
        mpt_options('lpsolver', lpsolver);
        break
    end
end

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, progress, prog_min, prog_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

H = [Aconstr;-f']; K = [Bconstr;100];
% P=polytope([Aconstr;-f'],[Bconstr;100],0,1); %normalize polytope
P = unitbox(2);
if(~isfulldim(P))
    disp('Positivity / Decay polytope is not full dimensional.')
    how='infeasible';
    drho=[];
    lyapunovL=[];
    lyapunovC=[];
    mpt_options('lpsolver', defaultLP);
    feasible = 0;
    return
else
    disp('Solving LP...')
%     [H,K]=double(P);
    try,
        [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f(:)',H,K,[],[],[],lpsolver);
        %[xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f(:)',H,K,[],[],chebyball(P),lpsolver);
    catch
        mpt_options('lpsolver', defaultLP);
    end
    drho=xopt(end);
end

% set the global LP solver back to default value
mpt_options('lpsolver', defaultLP);

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, 0.9, Options.status_min, Options.status_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

fprintf('\n')
feasible=0;
while(abs(feasible)<=1)
	if(~strcmp(how,'ok'))
        if Options.verbose>-1,
            disp(['mpt_getPWALyapFct:   Problem with LP - ' how])
        end
        feasible=feasible-1;
	elseif(drho>-Options.abs_tol)
        if Options.verbose>-1,
            disp(['PWA Lyapunov function not found; minimal decay bound is (almost) positive: rho=' num2str(drho)]);
        end
        [v,p]=max(H*xopt-K);
        if Options.verbose>-1,
            disp(['Problem in region/transition ' num2str(constr_origin{p})])
        end
        feasible=feasible-1;
	elseif(max(H*xopt-K)>Options.abs_tol)
        if Options.verbose>-1,
            disp(['LP constraints are violated by ' num2str(max(H*xopt-K))])
        end
        feasible=feasible-1;
    else
        %everything ok
        feasible=2;
	end
    
    if(abs(feasible)<=1)
        %try again
        if Options.verbose>-1,
            disp('Numerical Problems: tightening constraints and solving LP again...')
        end
        %LPi%[xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,H,K-Options.abs_tol,[],[],[],Options.lpsolver); 
        [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f',H,K-Options.abs_tol,[],[],[],lpsolver); 
        drho=xopt(end);
    end
end

if statusbar,
    if isempty(mpt_statusbar(Options.status_handle, 1, Options.status_min, Options.status_max)),
        mpt_statusbar;
        error('Break...');
    end     
end

lyapunovL=[];
lyapunovC=[];
if(feasible==2)
    feasible=1;
    if Options.verbose>-1,
        disp(['SUCCESS: System is asymptotically stable! Minimal decay bound ' num2str(drho)]);
    end
    %extract solution
	for i=1:lenP
    lyapunovL{i}=xopt(varindex{i}(1:n))';
    if(~containsOrigin(i))
        lyapunovC{i}=xopt(varindex{i}(n+1));
    else
        lyapunovC{i}=0;
    end
	end
else
    feasible=0;
    if Options.verbose>-1,
        disp('NO SUCCESS: System may NOT be asymptotically stable')
    end
end

if closestatbar,
    mpt_statusbar;
end
