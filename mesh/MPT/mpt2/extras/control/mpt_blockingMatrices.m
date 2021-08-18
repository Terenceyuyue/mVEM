function [Matrices]=mpt_blockingMatrices(Matrices,sysStruct,probStruct,Options)
%MPT_BLOCKINGMATRICES Constructs matrices for the CFTOC problem for move blocking strategies
%
% [Matrices]=mpt_blockingMatrices(Matrices,sysStruct,probStruct,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Constructs cost and constraint matrices for the finite time constrained 
% optimal control problem on linear time-invariant systems with quadratic
% cost (2-norm) in the case where inputs or their differences are fixed to be
% constant over a certain number of steps during the prediction horizon N.
% Matrices received from function mpt_constructMatrices are compressed
% depending on probStruct.inputblocking and probStruct.deltablocking. 
%
% Degrees of freedom are reduced from full degrees of freeedom m * N (m=number 
% of inputs, N=prediction horizon) in the non blocking case to m * M (M < N,
% M independent decision variables in blocked input sequence) in the case where
% inputs are their differences are fixed to be constant during the prediction
% horizon.
%
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Matrices          NORM=2
% --------
%                   G,W,E,H,F,Y,Cf,Cx,Cc  - matrices of the problem, i.e. 
%
% sysStruct        - System structure in the sysStruct format
%
% probStruct       - Problem structure in the probStruct format
%                  
%                   - NORM=2
%                                   
%                  - horizon          
%		             Prediction horizon N; How many time steps are considered.
%
%                  - inputblocking
%                    Fix inputs u_k to be constant over a certain number of steps
%                    in the prediction horizon.
%                     
%                    inputblocking=[n1 n2 ... nk], \sum_{n1}^{nk} = N.
%                    Entries define how many consecutive inputs are fixed to be 
%                    constant. Sum of all entries has to be equal the prediction
%                    horizon N.
%                     
%                    Example:      N = 5, inputblocking = [1 4]
%                                  First predicted input is independent and the
%                                  next 4 predicted inputs are fix to be constant,
%                                  i.e. u1, u2=u3=u4=u5.
%
%                  - deltablocking
%                    Fix difference u_k - u_{k+1} to be constant over a certain number 
%                    of steps in the prediction horizon.
%                       
%                    deltablocking = [1 ... nk ... N], 1 < .. < .. nk < .. < N
%                    Entries defines wich inputs are independent (=decision variables),
%                    inputs in between are interpolated (=constant differences).
%                    
%                    Example:      N = 5, deltablocking = [1 5]
%                                  First input u1 and last input u5 are independent,
%                                  inputs in between are interpolated, i.e.
%                                  u1-u2 = u2-u3 = u3-u4 = u4-u5.
%
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% Matrices              NORM=2
% --------
%                       Compressed matrices for fixed inputs or fixed differences in
%                       prediction horizon. 
%
%                       G,W,E,H,F,Y,Cf,Cx,Cc  - matrices of the problem, i.e. 
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% "Move Blocking Strategies in Receding Horizon Control"
% R. Cagienard, P. Grieder, E. C. Kerrigan and M. Morari, 2004, submitted
% check http://control.ee.ethz.ch for latest info
%
% see also MPT_CONSTRUCTMATRICES
%
% Copyright is with the following author(s):
%
% (C) 2004 Raphael Cagienard, Automatic Control Laboratory, ETH Zurich,
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          cagienard@control.ee.ethz.ch

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

error(nargchk(3,4,nargin));

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end
if nargin<4,
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
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'includeLQRset')
    Options.includeLQRset=1;
end
if ~isfield(sysStruct,'verified'),
    sysStruct=mpt_verifySysStruct(sysStruct);
end
if ~isfield(probStruct,'verified'),
    probStruct=mpt_verifyProbStruct(probStruct);
end
if ~isfield(Options,'pwa_index'),
    Options.pwa_index=1;  % which dynamics to extract
end
if ~isfield(Options,'noConstraintReduction')
    Options.noConstraintReduction=0;
end
if ~isfield(Options,'noReduce'),
    Options.noReduce=0;
end

if(sum(sysStruct.type~='LTI')~=0)
    fprintf('\n\n')
    error('mpt_blockingMatrices: no LTI-system')
end
if(probStruct.norm~=2)
    fprintf('\n\n')
    error('mpt_blockingMatrices: move blocking is implemented only for 2-norm')
end
if(probStruct.subopt_lev~=0)
    fprintf('\n\n')
    error('mpt_blockingMatrices: probStruct.subopt_lev not equals 0')
end


if (isfield(probStruct,'inputblocking') & isfield(probStruct,'deltablocking')),
    if(isempty(probStruct.inputblocking) & isempty(probStruct.deltablocking))
        disp('mpt_blockingMatrices: inputblocking and deltablocking vectors empty, returned without blocking');
        return
    end
    if(~isempty(probStruct.inputblocking) & ~isempty(probStruct.deltablocking))
        fprintf('\n\n')
        error('mpt_blockingMatrices: inputblocking and deltablocking vectors defined non empty')
    end
end
if (isfield(probStruct,'inputblocking') & ~isfield(probStruct,'deltablocking')),
    if(isempty(probStruct.inputblocking))
        disp('mpt_blockingMatrices: inputblocking empty and deltablocking not defined, returned without blocking');
        return
    end
end
if (~isfield(probStruct,'inputblocking') & isfield(probStruct,'deltablocking')),
    if(isempty(probStruct.deltablocking))
        disp('mpt_blockingMatrices: inputblocking not defined and deltablocking empty, returned without blocking');
        return
    end
end

G = Matrices.G;
W = Matrices.W;
E = Matrices.E;
H = Matrices.H;
F = Matrices.F;
Y = Matrices.Y;
Cf = Matrices.Cf;
Cx = Matrices.Cx;
Cc = Matrices.Cc;
symmetric = Matrices.symmetric;
bndA = Matrices.bndA;
bndb = Matrices.bndb;
Pinvset = Matrices.Pinvset;

% fix inputs to be constant
if (isfield(probStruct,'inputblocking') & ~isempty(probStruct.inputblocking)),
    if(sum(floor(probStruct.inputblocking)~=probStruct.inputblocking)~=0)
        fprintf('\n\n')
        error('mpt_blockingMatrices: inputblocking vector entries not integers')
    end
    if(sum(probStruct.inputblocking<=0))
        fprintf('\n\n')
        error('mpt_blockingMatrices: inputblocking vector entries zero or minus')
    end
    if(sum(probStruct.inputblocking)~=probStruct.N)
        fprintf('\n\n')
        error('mpt_blockingMatrices: sum of inputblocking vector entries not equal to prediction horizon N')
    end
    
    %output matlab command window
    if isfield(probStruct,'feedback')
        if(probStruct.feedback==0)
            vec=[' inputblocking: '];
            blocked_variable=['u'];
        elseif(probStruct.feedback==1)
            vec=[' offsetblocking: '];
            blocked_variable=['c'];
        end
    else
        vec=[' inputblocking: '];
        blocked_variable=['u'];
    end
    k=1;
    for i=1:length(probStruct.inputblocking)
        for j=1:probStruct.inputblocking(i)
            vec=[vec blocked_variable num2str(k) ];
            if((probStruct.inputblocking(i)-j)==0)
                vec=[vec ' '];
                k=k+1;
            else
                vec=[vec '='];    
                k=k+1;
            end
        end
        
    end
    disp(['mpt_blockingMatrices:' num2str(vec)])
    
    %compress matrices inputblocking
    nu=size(sysStruct.B,2);
    n_input=nu; 
    inputblocking=probStruct.inputblocking;
    for k=1:length(inputblocking)
        for i=1:(inputblocking(k)-1)
            
            kk=[(n_input*k - n_input + 1) : n_input*k];
            kkp= kk + ones(1,n_input)*n_input;
            
            H(:,kk)=H(:,kk)+H(:,kkp);
            H(kk,:)=H(kk,:)+H(kkp,:);
            H(:,kkp)=[];
            H(kkp,:)=[];
            F(:,kk)=F(:,kk)+F(:,kkp);
            F(:,kkp)=[];
            Cf(:,kk)=Cf(:,kk)+Cf(:,kkp);
            Cf(:,kkp)=[];
            G(:,kk)=G(:,kk)+G(:,kkp);
            G(:,kkp)=[];
        end
    end    
end

% fix differences to be constant
if (isfield(probStruct,'deltablocking') & ~isempty(probStruct.deltablocking)),
    if(isempty(probStruct.deltablocking))
        fprintf('\n\n')
        error('mpt_blockingMatrices: deltablocking vector is empty')
    end
    %error check on deltablockin definition vector
    if(sum(abs((floor(probStruct.deltablocking)-(probStruct.deltablocking)))~=0))
        fprintf('\n\n')
        error('mpt_blockingMatrices: deltablocking vector not integers')
    end
    if(min(probStruct.deltablocking)<1)
        fprintf('\n\n')
        error('mpt_blockingMatrices: deltablocking vector entry smaller than 1')
    end
    if(max(probStruct.deltablocking)>probStruct.N)
        fprintf('\n\n')
        error('mpt_blockingMatrices: deltablocking vector entry larger than prediction horizon N')
    end
    if(probStruct.deltablocking(1)~=1)
        fprintf('\n\n')
        error('mpt_blockingMatrices: deltablocking vector first entry not equal to 1')
    end
    if(probStruct.deltablocking(end)~=probStruct.N)
        fprintf('\n\n')
        error('mpt_blockingMatrices: deltablocking vector last entry not equal to prediction horizon N')
    end
    if(length(probStruct.deltablocking)~=length(unique(probStruct.deltablocking)))
        fprintf('\n\n')
        error('mpt_blockingMatrices: deltablocking vector contains several times same entry')
    end
    if(sum(sort(probStruct.deltablocking)~=probStruct.deltablocking)~=0)
        fprintf('\n\n')
        error('mpt_blockingMatrices: deltablocking vector entries not increasing')
    end
    
    %output matlab command window
    if isfield(probStruct,'feedback')
        if(probStruct.feedback==0)
            vec=[' deltainputblocking: '];
            blocked_variable=['u'];
        elseif(probStruct.feedback==1)
            vec=[' deltaoffsetblocking: '];
            blocked_variable=['c'];
        end
    else
        vec=[' deltainputblocking: '];
        blocked_variable=['u'];
    end
    for i=1:length(probStruct.deltablocking)
        vec=[vec blocked_variable num2str(probStruct.deltablocking(i)) ' '];
    end
    disp(['mpt_blockingMatrices:' num2str(vec) 'decision variables'])
     
    %compress matrices deltablocking
    if((max(probStruct.deltablocking)<=probStruct.N)|(length(probStruct.deltablocking)>1))
        
        nu=size(sysStruct.B,2);
        n_input=nu; 
        deltablocking=probStruct.deltablocking;
        for k=2:length(deltablocking)
            delta_N=deltablocking(k)-deltablocking(k-1);
            for i=(deltablocking(k-1)+1):(deltablocking(k)-1)
                               
                scalecoeff1=(deltablocking(k)-i)/delta_N;
                scalecoeff2=(i-deltablocking(k-1))/delta_N;
                                
                kk=[(n_input*deltablocking(k-1) - n_input + 1) : n_input*deltablocking(k-1)];
                kkp=[(n_input*deltablocking(k) - n_input + 1) : n_input*deltablocking(k)];
                kki=[(n_input*i - n_input + 1) : n_input*i];
                
                H(:,kk)=H(:,kk)+H(:,kki)*scalecoeff1;
                H(:,kkp)=H(:,kkp)+H(:,kki)*scalecoeff2;
                
                H(kk,:)=H(kk,:)+H(kki,:)*scalecoeff1;
                H(kkp,:)=H(kkp,:)+H(kki,:)*scalecoeff2;
                
                F(:,kk)=F(:,kk)+F(:,kki)*scalecoeff1;
                F(:,kkp)=F(:,kkp)+F(:,kki)*scalecoeff2;
                
                Cf(:,kk)=Cf(:,kk)+Cf(:,kki)*scalecoeff1;
                Cf(:,kkp)=Cf(:,kkp)+Cf(:,kki)*scalecoeff2;
                
                G(:,kk)=G(:,kk)+G(:,kki)*scalecoeff1;
                G(:,kkp)=G(:,kkp)+G(:,kki)*scalecoeff2;
                        
            end
        end
        
        kkt=[];
        for k=1:length(deltablocking)
            kkt=[kkt [(n_input*deltablocking(k) - n_input + 1) : n_input*deltablocking(k)]];
        end
        
        H=H(:,kkt);
        H=H(kkt,:);
        F=F(:,kkt);
        Cf=Cf(:,kkt);
        G=G(:,kkt);
      
    end
end

%this automatically removes redundant constraints in the formulation
%this needs to be switched off for the infinite horizon algorithm for LTI systems
if (sysStruct.type==0 | strcmp(sysStruct.type,'LTI')) & probStruct.N==inf & ~isfield(sysStruct,'guardC'),
    reduce_constraints=0;   
else
    reduce_constraints=1;
end
if Options.noConstraintReduction,
    reduce_constraints=0;
end

horizon=probStruct.N;    
nu=size(sysStruct.B,2);
nx=length(sysStruct.A);
ny=size(sysStruct.D,1);
nuH=nu*horizon;         %degrees of freedom = horizon  * number of inputs  

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%   Post-Processing of the Constraints 
%   (e.g. removal of redundant rows, reordering of constraints, etc.)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%Check if problem is symmetric / reorder constraints
%If the mp-solvers want to take advantage of symmetry, the constraints need to be ordered in a 
%predefined manner, i.e. they must alternate symmetrically. This is enforced by the lines below. 
GWE=[G W E];
tG=[];
tW=[];
tE=[];

foundmatch=0;
symmetric=probStruct.useSymmetry;
while(~isempty(GWE) & probStruct.useSymmetry)
    tG(end+1,:)   =GWE(1,1:nuH);
    tW(end+1,:)   =GWE(1,nuH+1);
    tE(end+1,:)   =GWE(1,(nuH+2):end);
    
    foundmatch=0;
    for i=2:size(GWE,1)
        normVec1=GWE(i,:);
        normVec1(find(isinf(normVec1)))=[];
        normVec2=[-tG(end,:) tW(end) -tE(end,:)];
        normVec2(find(isinf(normVec2)))=[];
        
        if(norm(normVec1-normVec2)<Options.abs_tol*size(GWE,2))
            foundmatch=1;
            tG(end+1,:)=GWE(i,1:nuH);
            tW(end+1,:)=GWE(i,nuH+1);
            tE(end+1,:)=GWE(i,(nuH+2):end);
            GWE([1 i]',:)=[];
            break
        end
    end
    if(foundmatch==0)
        if Options.verbose>1,
            symmetric=0;
            disp('mpt_constructMatrices: Constraints are not symmetric')
        end
        break
    end     
end
if(probStruct.useSymmetry)
    G=tG;
    W=tW;
    E=tE;
end
if Options.verbose>1,
    if(foundmatch==1)
        disp('mpt_constructMatrices: Constraints are symmetric')
    end
end



% Take out limits which are Inf
aux=find(isinf(W));
G(aux,:)=[];
W(aux,:)=[];
E(aux,:)=[];

% Take out possible constraints where the matrix [G K] has null rows
aux=find(all(([G E]==zeros(size([G E])))')'); % Rows which are all 0
G(aux,:)=[];
W(aux,:)=[];
E(aux,:)=[];

if(reduce_constraints)
    %remove all redundant constraints
    if Options.noReduce,
        GEW=polytope([G -E],W,0,2);
    else
        GEW=polytope([G -E],W);
    end
    if(~isfulldim(GEW)) 
        if Options.verbose>0,
            disp('mpt_constructMatrices: Problem is infeasible.')
        end
        G=zeros(1,nuH);
        E=zeros(1,nx);
        W=-Inf;
        S=E;
        H=-1;
        F=0;
        Hinv=1;
        Pinvset = polytope;
        symmetric = 0;
        if nargout==1,
            Matrices.G = G;
            Matrices.W = W;
            Matrices.E = E;
            Matrices.H = H;
            Matrices.F = F;
            Matrices.Y = Y;
            Matrices.Cf = Cf;
            Matrices.Cx = Cx;
            Matrices.Cc = Cc;
            Matrices.symmetric = symmetric;
            Matrices.bndA = bndA;
            Matrices.bndb = bndb;
            Matrices.Pinvset = Pinvset;
            Matrices.GEW = GEW;
            G = Matrices;
        end
        return
    end
    [GE,W]=double(GEW);
    G=GE(:,1:(end-nx));
    E=-GE(:,(end-nx+1):end);
end


if(probStruct.norm==2)
    H=(H+H')/2;      %make sure hessian is symmetric
    Hinv=inv(H);
    S=E+G*Hinv*F';
    
    if(cond(H)>1e8)
        disp('mpt_constructMatrices:')
        disp(' +++++++++++++    WARNING     WARNING     WARNING     WARNING     +++++++++++++')
        disp('                       THE CONDITION NUMBER OF THE H MATRIX IS VERY LARGE              ')
        disp('                   THIS MAY LEAD TO NUMERICAL ERRORS AND INCONSISTENCIES        ')
        disp(' CHECK THE VALIDITY OF THE OBTAINED RESULTS AFTER THE ALGORITHM HAS COMPLETED')
    end
    
else
    %set to empty for 1 / Inf norm
    Hinv=[];
    S=[];
end


if nargout==1,
    Matrices.G = G;
    Matrices.W = W;
    Matrices.E = E;
    Matrices.H = H;
    Matrices.F = F;
    Matrices.Y = Y;
    Matrices.Cf = Cf;
    Matrices.Cx = Cx;
    Matrices.Cc = Cc;
    Matrices.symmetric = symmetric;
    Matrices.bndA = bndA;
    Matrices.bndb = bndb;
    Matrices.Pinvset = Pinvset;
    G = Matrices;
end


return
