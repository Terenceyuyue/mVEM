function [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver1(Matrices,Options)
%MPT_MPLP Explicitly solves the given linear program (LP)
%
% [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp(Matrices,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Multiparametric linear programming
%
% Solves the problem
%   V(x) = min H U
%           U
%   s.t.   G U <= W + E x
%          bndA*x <= bndb
%
% As a solution we get 'nR' regions
%   Pn(i)={x : H x <= K}
% 
% with the optimal control law
%   U = Fi{i} x + Gi{i}
%
% and the corresponding cost function expression
%   V(x) = Bi{i} x + Ci{i}
%  
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Matrices - a struct with all the parameters which are needed.
%            See description above for explanation.
%   Matrices.H=H;
%   Matrices.G=G;   
%   Matrices.W=W;
%   Matrices.E=E;
%   Matrices.bndA=bndA;   Limits on exploration space, i.e. bndA*x<=bndb
%   Matrices.bndb=bndb;
%
% Options.verbose      - level of verbosity
% Options.lpsolver     - which LP solver to use (help mpt_solveLP)
% Options.max_iter     - maximum number of iterations of the algorithm
% Options.step_size    - length of step over a facet
% Options.f_perturb    - Perturbation of the optimization direction
% Options.nu           - How many elements to extract from the optimizer (to
%                        deal with slacks) 
% Options.debug_level  
%           Due to numerical problems tiny regions are sometimes difficult to
%           calculate, i.e. are not identified at all. This may create "gaps"
%           in the computed control law. For the exploration, these will be
%           jumped over and the exploration in the state space will continue.
%           "debug_level" can have three values:
%       
%           0: No debug done
%           1: A tolerance is given to find gap in the region partition,
%              small empty regions inside the region partition will be discarded.
%              Note that this is generally not a problem, since the feedback law 
%              is continuous and can therefore be interpolated easily.
%              Correction to the calculation of the outer hull.
%           2: Zero tolerance to find gap in the region partition, empty regions
%              if they exist, will be detected, i.e. the user will be notified.
%              Correction to the calculation of the outer hull.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Pn,Fi,Gi           - for region Pn(i).H*x <= Pn(i).K computed input is
%                      U=Fi{i}*x+Gi{i}   
% activeConstraints  - Cell Array which stores the active constraints 
%                      of the optimizer in each region.
% Phard              - Defines the feasible state space partition (i.e. union of
%                      all regions) as Phard.H*x<=Phard.K
% details            - a structure with the following fields:
%     nR      number of regions
%     Pn      polyhedral partition
%     Fi      control law
%     Gi      control law
%     BC      connection list
%     Bi      value function
%     Ci      value function
%     nHard   number of hard constraints
%     Phard   polytope given by hard constraints
%     nb      number of constraints for each region
%     LISTa   list of active constraints
%
% see also MPT_CONSTRUCTMATRICES, MPT_MPQP, MPT_OPTCONTROL, MPT_OPTCONTROLPWA

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%     baotic@control.ee.ethz.ch
% (C) 2002 Francesco Borrelli, Automatic Control Laboratory, ETH Zurich

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
if nargin<2,
    Options=[];
end

if ~isfield(Options,'debug_level'),
    Options.debug_level=mptOptions.debug_level;
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'rel_tol'),
    Options.rel_tol=mptOptions.rel_tol;
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'max_iter'),
    Options.max_iter=20;       % maximum number of iterations of the algorithm
end
if ~isfield(Options,'f_perturb'),
    Options.f_perturb=0.5;     % Perturbation of the optimization direction up to 50% 
end
if ~isfield(Options,'step_size'),
    Options.step_size=mptOptions.step_size;
end
if ~isfield(Options,'ispwa'),
    Options.ispwa=0;
end
if ~isfield(Options,'nu'),
    % how many elements to extract from the optimizer (to deal with slacks)
    Options.nu=size(Matrices.G,2);
end

lpsolver=Options.lpsolver;
EMPTY_ROW_TOL=Options.abs_tol;  % Tolerance for declaring that row is empty
CONSTR_TOL = Options.rel_tol;   % Tolerance for declaring that constr. is redundant
ALPHA = max(Options.step_size,1e-7);
ALPHAit  = 50;       % number of iterations
ALPHAmax = Options.step_size;     % maximum step
ALPHAinc = ALPHAmax/(ALPHA*ALPHAit);        

DEBUG=(Options.debug_level==2); % DEBUG is true if debug_level==2

emptypoly=polytope;

L1=Matrices.H;
G=Matrices.G;
S=Matrices.E;
W=Matrices.W;
bndA=Matrices.bndA;
bndb=Matrices.bndb;



nx=size(S,2);
nu=size(G,2);
nC=size(G,1);

if(~isfield(Options,'center'))
    center = zeros(nx,1);
else
    center = Options.center;
end


nR = 0;     % number of regions

nHard = 0;  % number of hard constraints
hardA=[];   % hard constraints
hardb=[];   % hard constraints
nb=[];      % number of constraints for each region

LISTa={};   % list of active constraints

Pn=emptypoly;
Hn={};      % normalized polyhedron description
Kn={};      % normalized polyhedron description
Fi={};      % control law
Gi={};      % control law
BC={};      % connection list
Bi={};      % value function
Ci={};      % value function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                   %%
%%      FIND THE STARTING POINT      %%
%% (IF ORIGINAL PROBLEM IS FEASIBLE) %%
%%                                   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(bndA),
    crA=[G -S];
    crb=W;
else
    crA=[G -S; zeros(size(bndA,1),size(G,2)+size(S,2)-size(bndA,2)) bndA];
    crb=[W; bndb];
end

Pcr=polytope(crA, crb);
opt_poly.lpsolver=lpsolver;
[aux,R]=chebyball(Pcr,Options);
indeq=[];


if R<=EMPTY_ROW_TOL,
    Ri={};
    nR=0;
    nHard=0;
    hardA=[];
    hardb=[];
    disp('mpt_mplp: No feasible starting point found');
    
    Ri.LISTa=LISTa;
    Ri.nb=nb;
    Ri.Pn=polytope;
    Ri.Fi=Fi;
    Ri.Gi=Gi;
    Ri.Bi=Bi;
    Ri.Ci=Ci;
    Ri.BC=BC;
    Ri.nR=nR;
    Ri.nHard=nHard;
    Ri.Phard = polytope;
    activeConstraints=Ri.LISTa;
    Phard=Ri.Phard;
    for i=1:Ri.nR
        Ri.Bi{i}=L1(:)'*Ri.Fi{i};
        Ri.Ci{i}=L1(:)'*Ri.Gi{i};
    end
    details=Ri;
    Pn=Ri.Pn;
    Fi=Ri.Fi;
    Gi=Ri.Gi;
    return;
end
xFeasible=aux(nu+1:nu+nx);

%================================================%
% Solve optimization problem for point xFeasible %
%================================================%

if ~isfield(Options,'max_iter'),
    Options.max_iter=7;       % maximum number of iterations of the algorithm
end
if ~isfield(Options,'f_perturb'),
    Options.f_perturb=0.5;     % Perturbation of the optimization direction up to 50% 
end
[found,ii,indexmaxdiff]=findas(L1,G,W,S,[],xFeasible,[],Options);

if found<1
    % no active constraints were found for the given initial point xFeasible
    % we try to perturb it
    x0=xFeasible;
    indexmaxdiff=[];
    Pub=unitbox(nx,1);
    Vub=extreme(Pub);
    iter=0;
    while iter<=Options.max_iter   % try at most Options.max_iter times
        iter=iter+1;
        for jj=1:2^nx,
            [found,ii,indexmaxdiff]=findas(L1,G,W,S,indeq,x0,indexmaxdiff,Options);
            %found=1=> region not flat and optimal vertex found.
            %found=0=> region flat
            %found=-1=> problem infeasible
            if found
                break;
            end
            if isinf(R)
                R = 1e-1;
            end
            % ad-hoc strategy to perturb a state vector
            pert = (4+Options.f_perturb/(jj+1))/(Options.max_iter+2-iter);
            F=Vub(jj,:);
            F(1)=sign(xFeasible(1))*F(1)*Options.f_perturb/jj;
            if nx>1,
                F(2)=F(2)*Options.f_perturb/2*jj;
            end
            addvec = F'*pert;
            x0=xFeasible+R*addvec; 
        end
        if found,
            break
        end
    end
end

if found<1
    % ad-hoc perturbation failed, feasible starting point was still not found
    % we use random perturbations
    x0=xFeasible;
    disp('mpt_mplp: using random perturbation to find feasible starting point. Results may differ!');
    iter=0;
    while iter<=Options.max_iter
        iter=iter+1;
        [found,ii,indexmaxdiff]=findas(L1,G,W,S,indeq,x0,indexmaxdiff,Options);
        if found
            break;
        end
        if isinf(R)
            R = 1e-1;
        end
        x0=xFeasible+R*(2*abs(rand(nx,1))-1); 
    end
end

if found<=0,
    Ri={};
    nR=1;
    nHard=nC;
    hardA=bndA;
    hardb=bndb;
    disp('No feasible starting point in mpt_mplp');
    
    %This is if you what to store the initial region where or the problem is infeasible or there are only flat CR
    Ri.LISTa{nR}=LISTa;
    Ri.nb=nb;
    Ri.Pn=polytope;
    Ri.Fi{nR}=inf(nu,nx);
    Ri.Gi{nR}=inf(nu,1);
    Ri.Bi{nR}=inf(nx,1);
    Ri.Ci{nR}=Inf;
    Ri.BC{nR}=BC;
    Ri.nR=0;
    Ri.nHard=nHard;
    Ri.Phard = polytope(hardA, hardb);
    activeConstraints=Ri.LISTa;
    Phard=Ri.Phard;
    for i=1:Ri.nR
        Ri.Bi{i}=L1(:)'*Ri.Fi{i};
        Ri.Ci{i}=L1(:)'*Ri.Gi{i};
    end
    details=Ri;
    Pn=Ri.Pn;
    Fi=Ri.Fi;
    Gi=Ri.Gi;
    return;
end

nii=length(ii);

nR=nR+1; % new region
LISTa{nR}=ii;

if nii==0,
    error('mpt_mplp: Critical error detected! Check your constraint and problem matrices.');
else
    
    % STEP 3: Compute matrices the solution X and the polyhedra crA crb in this way:
    %         apply Gauss to the following system 
    %         Gtilde   -St     |z  =  Wt  
    %                          |x  
    %           to obtain
    %        
    %         Iz       0    L     |z    =  |Wt 
    %         0        Ix   N     |x      
    %
    %
    %
    Gt=G(ii,:);
    Gcap=G;
    Gcap(ii,:)=[];
    Wt=W(ii);
    Wcap=W;
    Wcap(ii,:)=[];
    St=S(ii,:);
    Scap=S;
    Scap(ii,:)=[];
    
    BigA0=[Gt -St];
    BigB0=[Wt];
    sA0=size(BigA0,2);
    [Alo,Qlo]=qr(BigA0);
    BigB0=Alo\BigB0;
    T=Qlo(1:nu,1:nu);
    Af=inv(T)*Qlo(1:nu,:);
    bf=inv(T)*BigB0(1:nu,:);
    
    
    % Linear term of solution (lambda does not depend on x)
    Dx = -Af(1:nu,nu+1:nu+nx);
    % Constant term
    cx = bf;
    tsA = Gcap*Dx- Scap; 
    tsb = Wcap - Gcap*cx;   
    Fi{nR}=Dx;      % control law
    Gi{nR}=cx;      % control law
    
    crA=[bndA;  tsA];   % region definition
    crb=[bndb;  tsb];   % region definition
end

Pcr=polytope(crA, crb);
Pn=[Pn Pcr];
isemptypoly=~isfulldim(Pcr);

if isemptypoly,
    if Options.verbose>0,
        disp('mpt_mplp: Empty polyhedron after removal');
    end
end

Pquadrant = {};
Pquadrant = cell(1,2^nx); % number of quadrant is 2^nx
for ii=1:length(Pquadrant),
    Pquadrant{ii} = emptypoly;
end
q = sub_whichquadrant(Pcr,center);
for qqq=1:length(q),
    Pquadrant{q(qqq)} = [Pquadrant{q(qqq)} Pcr];
end

nb(nR)=nconstr(Pn(nR));
BC{nR}=zeros(nb(nR),1);


%%%%%%%%%%%%%%%%%%%%%%%%
%%                    %%
%% DO THE EXPLORATION %%
%%                    %%
%%%%%%%%%%%%%%%%%%%%%%%%
region=1;
while region<=nR,
    if Options.verbose>1,
        fprintf('Region: %d/%d, %d borders, %d hard        \r', region,nR,nb(region),nHard);
    elseif Options.verbose==1,
        if mod(region,20)==0,
            fprintf('Region: %d/%d        \r', region,nR);
        end
    end
    
    for border=1:nb(region),
        
        % Create a point close to the border
        [xBorder, RBorder] = facetcircle(Pn(region), border, Options);
        [Hnr, Knr] = double(Pn(region));
        xBeyond=xBorder+ALPHA*Hnr(border,:)';
        
        % Check if the point is inside of any existing polyhedra
        nMatch=0;
        isinOpt.abs_tol = 0;
        isinOpt.fastbreak = 1;
        vi = (xBeyond - center) > 0;
        q = 1;
        for i=0:nx-1
            q = q + (vi(i+1))*2^i;    
        end
        %[dummy,neighbor] = isinside(Pn,xBeyond,isinOpt);
        [dummy,neighbor] = isinside(Pquadrant{q},xBeyond,isinOpt);
        
        nMatch = length(neighbor);
        if nMatch>0
            neighbor = neighbor(end);
        end
        
        if nMatch>=1,
            if nMatch>1,
                if Options.verbose>0 & ~Options.ispwa,
                    disp('Polyhedra should not overlap');
                end
            end
            BC{region}(border)=neighbor;
            
        elseif nHard==0 | all(hardA*xBeyond <= hardb),
            
            %==============================================%
            % Solve optimization problem for point xBeyond %
            %==============================================%
            
            niter=1;
            found=0;
            while found==0 & niter<=ALPHAit,
                [found,ii,indexmaxdiff]=findas(L1,G,W,S,indeq,xBeyond,indexmaxdiff,Options);
                %if found=0 then 1. alpha is too small or 
                %                2.x beyound is on a facet of a different critical region
                
                xBeyond=xBorder+ALPHA*ALPHAinc*niter*Hnr(border,:)'; 
                niter=niter+1;
            end
            nii=length(ii);
            if found==0,
                disp('Storing a flat region');
            end
            
            
            if found==-1 | isequal(ii,LISTa{region}),
                nHard=nHard+1;
                hardA(nHard,:)=Hnr(border,:);
                hardb(nHard,:)=Knr(border,:);
            else
                nR=nR+1; % new region
                BC{region}(border)=nR;
                LISTa{nR}=ii;
                
                
                
                if nii==0,
                    error('mpt_mplp: Critical error detected! Check your constraint and problem matrices.');
                else
                    
                    % STEP 3: Compute matrices the solution X and the polyhedra crA crb in this way:
                    %         apply Gauss to the following system 
                    %         Gtilde   -St     |z  =  Wt  
                    %                          |x  
                    %           to obtain
                    %        
                    %         Iz       0    L     |z    =  |Wt 
                    %         0        Ix   N     |x      
                    %
                    %
                    %
                    Gt=G(ii,:);
                    Gcap=G;
                    Gcap(ii,:)=[];
                    Wt=W(ii);
                    Wcap=W;
                    Wcap(ii,:)=[];
                    St=S(ii,:);
                    Scap=S;
                    Scap(ii,:)=[];
                    
                    BigA0=[Gt -St];
                    BigB0=[Wt];
                    sA0=size(BigA0,2);
                    [Alo,Qlo]=qr(BigA0);
                    BigB0=Alo\BigB0;
                    T=Qlo(1:nu,1:nu);
                    Af=inv(T)*Qlo(1:nu,:);
                    bf=inv(T)*BigB0(1:nu,:);
                    
                    
                    % Linear term of solution (lambda does not depend on x)
                    Dx = -Af(1:nu,nu+1:nu+nx);
                    
                    % Constant term
                    cx = bf;
                    
                    tsA = Gcap*Dx- Scap; 
                    tsb = Wcap - Gcap*cx;   
                    
                    
                    Fi{nR}=Dx;      % control law
                    Gi{nR}=cx;      % control law
                    
                    crA=[bndA;  tsA];   % region definition
                    crb=[bndb;  tsb];   % region definition
                end            
                
                %Add the border to handle overlapping regions with dual degeneracy
                
                crA=[crA; -Hnr(border,:)];
                crb=[crb; -Knr(border)];
                Ptemp = polytope(crA, crb);
                isemptypoly=~isfulldim(Ptemp);
                
                if isemptypoly,
                    if Options.verbose>0,
                        disp('mpt_mplp: Empty polyhedron after removal. Removing this polyhedron from the candidate list');
                    end
                    
                    %added by Mato Baotic, 05.11.2002
                    BC{region}(border)=0;
                    LISTa{nR}=[];
                    nR=nR-1;
                else
                    Pn = [Pn Ptemp];
                    q = sub_whichquadrant(Ptemp,center);
                    for qqq=1:length(q),
                        Pquadrant{q(qqq)} = [Pquadrant{q(qqq)} Ptemp];
                    end
                    nb(nR)=nconstr(Ptemp);
                    BC{nR}=zeros(nb(nR),1);
                end
                
            end
        end % nMatch==0
        
    end % border
    region=region+1;
    
end % END EXPLORATION




Ri.LISTa=LISTa;
Ri.nb=nb;
Ri.Pn=Pn;
Ri.Fi=Fi;
Ri.Gi=Gi;
Ri.Bi=Bi;
Ri.Ci=Ci;
Ri.BC=BC;
Ri.nR=nR;
hardA=[hardA; bndA];
hardb=[hardb; bndb];
Ri.Phard = polytope(hardA, hardb);
Ri.nHard = nconstr(Ri.Phard);
Ri.activeConstraints = Ri.LISTa;
activeConstraints=Ri.LISTa;
Phard=Ri.Phard;
for i=1:Ri.nR
    Ri.Bi{i}=L1(:)'*Ri.Fi{i};
    Ri.Ci{i}=L1(:)'*Ri.Gi{i};
    Ri.Fi{i}=Ri.Fi{i}(1:Options.nu,:);
    Ri.Gi{i}=Ri.Gi{i}(1:Options.nu,:);
end
Fi=Ri.Fi;
Gi=Ri.Gi;
details=Ri;
if Options.verbose>0,
    disp(sprintf('mpt_mplp: %d regions', region-1));
end


return;



function [out,ii,indexmaxdiff]=findas(f,A,b,F,indeq,tq,indexmaxdiff,Options)
% Find set of active constraints.

% out=0 => tq belong to a flat region
% out=1 => tq belong to a no-flat region and ii is list of active constraints defining a base
% out=-1 => problem infeasible
% f - cost
% A z <= b+F tq

xN=size(F,2);
Aeq=A(indeq,:);
Aineq=A;
Aineq(indeq,:)=[];
bigb2 = b + F*tq; 
% solve a smaller system just involving X variables
ii=[];      %Index of active constraints
out=0;         %out=1 if the algorithms finds a basis for for theta=tq.

if ~isempty(indeq),
    bigb2eq=bigb2(indeq);
    bigb2ineq=bigb2;
    bigb2ineq(indeq)=[];
else
    bigb2ineq=bigb2;
    bigb2eq=[];
end

[xq,fval,lambda,exitflag,how]=mpt_solveLPi(f(:)',Aineq,bigb2ineq,Aeq,bigb2eq,[],Options.lpsolver);

if ~strcmp(how,'ok'), 
    if Options.verbose>1,
        disp('mpt_mplp: Infeasible LP subproblem #2')
    end
    out=-1;
    return 
else 
    %the problem is feasible so find the the set of active constraints:
    ii=find(A*xq-bigb2>=-Options.abs_tol);% Index of active constraints

    if rank(A(ii,:))<size(A,2)  %Possible dual degeneracy.
        if Options.verbose>=2,
            disp('mpt_mplp: Possible Dual Degeneracy: trying to find a basis');
        end
        
        ii=find(A*xq-bigb2>=-Options.abs_tol);% Index of active constraints
        if rank(A(ii,:))<size(A,2)  %Possible dual degeneracy.
            if Options.verbose>2,
                disp('mpt_mplp: Possible Dual Degeneracy: trying new method to find a basis');
            end
            newAeq=Aeq;
            newAineq=Aineq;
            newbigb2eq=bigb2eq;
            newbigb2ineq=bigb2ineq;
            in=1;
            jcou=1;
            xqold=xq;
            while jcou<10 & in;
                ii=find(newAineq*xq-newbigb2ineq>=-Options.abs_tol);
                if isempty(ii),
                    in=0;
                end
                newAeq=[newAeq;newAineq(ii,:)];
                newbigb2eq=[newbigb2eq;newbigb2ineq(ii,:)];
                newAineq(ii,:)=[];
                newbigb2ineq(ii,:)=[];
                
                [xq,fval,lambda,exitflag,how]=mpt_solveLPi(f(:)',newAineq,newbigb2ineq,newAeq,newbigb2eq,[],Options.lpsolver);
                
                mi=find(abs(xq-xqold)>100*Options.abs_tol);
                for j=1:length(mi)
                    if ~any(mi(j)==indexmaxdiff),
                        %if ~ismember(mi,indexmaxdiff),
                        %%indexmaxdiff=[indexmaxdiff;mi];
                        indexmaxdiff=[indexmaxdiff;mi(j)];
                    end
                end
                ii=find(A*xq-bigb2>=-Options.abs_tol);% Index of active constraints
                if rank(A(ii,:))==size(A,2),
                    in=0;
                end
                jcou=jcou+1;
                xqold=xq;
            end
            if rank(A(ii,:))<size(A,2)  %Possible dual degeneracy.
                if Options.verbose>0,
                    disp('Possible Dual Degeneracy: trying to modify f in the direction of degeneracy');
                end
                for qq=1:10,
                    degenerate = 1;
                    for k=1:length(indexmaxdiff);
                        f(indexmaxdiff(k))=rand(1);
                    end 
                    %test: perturb it all!
                    fpert=rand(length(f),1);
                    f=f(:)+fpert(:);
                    
                    [xq,fval,lambda,exitflag,how]=mpt_solveLPi(f(:)',newAineq,newbigb2ineq,newAeq,newbigb2eq,[],Options.lpsolver);
                    
                    ii=find(A*xq-bigb2>=-100*Options.abs_tol);% Index of active constraints
   
                    if rank(A(ii,:))<size(A,2) %Possible dual degeneracy.
                        continue
                        error('Not possible to find a basic solution. Try to modify Options.abs_tol')
                    end
                    degenerate = 0;
                end
                if degenerate,
                    error('Not possible to find a basic solution. Try to modify Options.abs_tol')
                end
            end
        end
    end
end

P=A(ii,:);
L=F(ii,:);
if rank([P L])>rank(P),
    %I'm on the boundary of CR
    if Options.verbose>1,
        disp('mpt_mplp: solution is on the boundary of 2 or more CR')
    end
    return
end
out=1;


function q = sub_whichquadrant(P,center)

Options.noPolyOutput = 1;
[R,l,u] = bounding_box(P,Options);

nx = dimension(P);
vert = ones(nx,2^nx);
binOne=dec2bin(1);
for i=1:2^nx
    decVec=dec2bin(i-1,nx);
    for j=1:nx
        if(decVec(j)==binOne)
            vert(j,i)=l(j);
        else
            vert(j,i)=u(j);
        end
    end
end

q1 = zeros(1,2^nx);
for i=1:2^nx
    q2 = sub_whichquadrant_point(vert(:,i),center);
    q1(q2) = 1;
end
q = find(q1==1);

return



%New function to find to which quadrant a point belongs

function q = sub_whichquadrant_point(xBeyond,center)

q = [];
q1 = [];

nx = length(xBeyond);
%center = zeros(nx,1);       % center per default
twoQuadrant = find(xBeyond==0);

if isempty(twoQuadrant)   
    %One = eye(nx);
    vi = (xBeyond - center) > 0;
    q1 = 1;
    for i=0:nx-1
        q1 = q1 + (vi(i+1))*2^i;    
    end   
else    
    for j=1:length(twoQuadrant)
        xBeyond(twoQuadrant(j))=1;
        %One = eye(nx);
        vi = (xBeyond - center) > 0;
        q = 1;
        for i=0:nx-1
            q = q + (vi(i+1))*2^i;    
        end 
        
        q1 = [q1 q];
        xBeyond(twoQuadrant(j))=-1;
        %One = eye(nx);
        vi = (xBeyond - center) > 0;
        q = 1;
        for i=0:nx-1
            q = q + (vi(i+1))*2^i;    
        end 
        xBeyond(twoQuadrant(j))=0;
        q1 = [q1 q]; 
    end      
end
q = q1;

% return
