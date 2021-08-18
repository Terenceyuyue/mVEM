function [res]=regiondiffXU(Pbnd,Pn,G,W,E,Options)
%REGIONDIFFXU Computes region difference in lifted XU space
%
% R=regiondiffXU(Pbnd,Pn,G,W,E,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Function computes the difference between polytope Pbnd
% and union of polytopes Pn:  res=Pbnd\(U Pn) in the lifted XU space.
% In general output "res" is a polytope array containing
% non-overlapping polytopes whose union is equal to
% the difference Pbnd\(U Pn).
%
% Note:
% Algorithm works even for a non-minimal representation of Pn,
% However, if Pn are minimal, then code runs faster...
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pbnd                  - polytope
% Pn                    - array of polytopes
% G,W,E                 - Matrices of the mpQP optimization problem
% Options.lpsolver      - Solver for LPs (see help mpt_solveLP)
% Options.verbose       - level of verbosity (see help mpt_init)
% Options.intersect_tol - tolerance to detect intersection
% Options.abs_tol       - absolute tolerance
% Options.tolerance     - nonegative tolerance to decide if constraint is redundant
% Options.reduce        - if region R is not given in minimal representation
%                         this Option should be set to 1, default=1 (yes)
% Options.reduce_output - should the description of regions Pj={x | P.Hn{j} x <= P.Kn{j}}
%                         j=1,...,P.nR, be reduced; default=1 (yes)
% Options.constructset  - if difference is not convex, should set P be
%                           constructed; default=1 (yes)
%
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% res   - polytope array containing the region difference Pbnd\(U Pn)
%         if the difference is empty (i.e. Pbnd is fully covered by (U Pn))
%         returns an empty polytope
%
% see also REGIONDIFF
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch
%
% History
% --------
% 11.10.2002    Initial version
% 30.07.2003    Options.reduce_output added
% 02.08.2003    Unnecessary cuts are avoided. Huge speed up of the code.

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

%emptypoly=polytope;
emptypoly = mptOptions.emptypoly;
if nargin<4,
      error('mpt_regiondiffXU: Not enough input arguments in mpt_regiondiffXU!');
end
if(nargin<5)
    Pbnd=emptypoly;
end
if nargin<6,
    Options = [];
end
if ~isfield(Options,'lpsolver'),
   Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'intersect_tol'),
    Options.intersect_tol=sqrt(eps);
end
if ~isfield(Options,'reduce'),
    Options.reduce=0;         % 1 = reduce the description of regions R
end
if ~isfield(Options,'constructset'),
    Options.constructset=0;   % should we construct the difference, or just see if U Ri covers R
end
if ~isfield(Options,'presolved'),
    Options.presolved=0;      % have "active constraints" of Ri already been isolated
end
if ~isfield(Options,'checkoverlap'),
    Options.checkoverlap=1;      % remove elements of Hn and Kn which do not intersect polytope given by H and K 
end
if ~isfield(Options,'reduce_output'),
    Options.reduce_output=0;      % remove elements of Hn and Kn which do not intersect polytope given by H and K 
end

abs_tol = Options.abs_tol;

%%Translate Pn into lifted (X,U)-space
nu=size(G,2);
if isempty(Pn.Array),
    Pn.Array{1}=Pn;
end
for i=1:length(Pn.Array)
    Pn.Array{i}.H=[Pn.Array{i}.H zeros(size(Pn.Array{i}.H,1),nu)];
end
if(isfulldim(Pbnd))
    [bndA,bndb]=double(Pbnd);
    bndA=[bndA zeros(size(bndA,1),nu)];
    Pbnd=set(Pbnd,'H',bndA);
else
    bndA=[];
    bndb=[];
end
%Obtain key lifted  polytope H*[x U]<=K
%%Phk = polytope([-E G; bndA],[W;bndb]);
Hk = [-E G; bndA];
Kk = [W; bndb];
[xHK, rHK] = chebyball_f(Hk,Kk,Options);
if rHK<Options.abs_tol,
    res=emptypoly;
    return
end

% initial value for the output
%-----------------------------
res = emptypoly;

% check which Ri intersect R
%---------------------------
N=length(Pn);   % initial number of regions
Rc=zeros(N,1);  % radius of Chebyshev ball
Options.reduce_intersection=0;
for ii=1:N
    %% [Hk,Kk]=double(Phk);
    Hn=Pn.Array{ii}.H;
    Kn=Pn.Array{ii}.K;
    [xc,Rc(ii)]=chebyball_f([Hk; Hn], [Kk; Kn], Options);
end
[val,ind]=sort(-Rc);            % sort Ri in descending order


N=sum(Rc>=Options.intersect_tol);   % consider only Ri that do intersect R

% no Ri intersect R => return R
%------------------------------
if N==0
    res = polytope(Hk,Kk);
    return;
end    

Pdummy = Pbnd;

A=[-E G; bndA];
B=[W;bndb];
H=A;
K=B;
m=size(A,1);      % number of rows in R
mi=zeros(N,1);    % beginning of the storage of the region Ri constraints



% Find "active constraints" of polyhedra Ri and put them all in one big matrix
% "active constraints" are those that do not show up in the description of R
%-----------------------------------------------------------------------------
if Options.presolved
    % "active constraints" separated in advance (description of Ri comprises only such constraints)
    for ii=1:N
        i=ind(ii);
        mi(ii)=nconstr(Pn(i));
        Phk=Phk & Pn(i);
    end
else
    % find "active constraints"
    for ii=1:N
        i=ind(ii); 
        %%[Hni,Kni]=double(Pn(i));
        Hni=Pn.Array{i}.H;
        Kni=Pn.Array{i}.K;
        for j=1:size(Hni,1)
            % find unique facets (different from those defining R)
            HKnij = [Hni(j,:) Kni(j)];
            if all( sum( abs([H,K]-HKnij(ones(m,1),:)), 2) >= abs_tol )
                % was before:
                % if all( sum( abs([H,K]-repmat([Hni(j,:) Kni(j)],m,1)), 2) >= Options.abs_tol )
                mi(ii)=mi(ii)+1;
                A=[A; Hni(j,:)];
                B=[B; Kni(j)];
            end
        end
    end
end

% if some Ri has no "active constraints" we are done (Ri is covering R)
%----------------------------------------------------------------------
if any(mi==0)
    res=emptypoly;
    return;
end


M=sum(mi);                          % total number of constraints
beg_mi=m+cumsum([0; mi(1:end-1)]);  % beginning position of constraints description for all Ri

A=[A; -A(m+1:m+M,:)];               % constraints negative expressions (for faster computation afterwards)
B=[B; -B(m+1:m+M,:)];               % constraints negative expressions (for faster computation afterwards)


counter=zeros(1,N);                 % constraints counters
INDICES=1:m;                        % indices of currently valid constraints

level=1;                            % level in the tree (= index of the Ri currently considered)

%===================
% Explore the tree
%===================
while level~=0
    
    
    
    if counter(level)==0
        % if we are for the first time on this level
        % first check if Ri{level} is intersecting particular node (described through INDICES)
        % if Ri{level} is not intersecting this node increase the level (go to the next Ri)
        % until we find one Ri that does intersect the node
        %=====================================================================================
        for j=level:N
            auxINDICES=[INDICES beg_mi(j)+1:beg_mi(j)+mi(j)];
            [xopt,R]=chebyball_f(A(auxINDICES,:),B(auxINDICES),Options);
            if R>=abs_tol
                % we do intersect the node on a level j: start exploring from there
                %------------------------------------------------------------------
                level=j;
                counter(level)=1;
                INDICES=[INDICES beg_mi(level)+1+M];
                break;
            end
        end
        % if no Ri intersect the node, go one level up
        % use node description (from that higher level)
        % and store it as a leaf node => part of R\(U Ri)
        %------------------------------------------------
        if R<abs_tol
            level=level-1;
            % consider this to be a leaf node. store result for the output
            %--------------------------------------------------------------
            if Options.reduce_output
                res = [res polytope(A(INDICES,:),B(INDICES))];
            else
                %%res = [res polytope(A(INDICES,:),B(INDICES),0,1)];
                res = Pdummy; 
            end
            if ~Options.constructset
                return;                 % WE ARE NOT LOOKING FOR THE FULL SOLUTION -> END SEARCH
            end                
            
            % inrease counters
            %-----------------
            nzcount=find(counter>0); % find nonzero counters
            for jj=length(nzcount):-1:1
                level=nzcount(jj);
                counter(level)=counter(level)+1;
                if counter(level)<=mi(level)
                    INDICES(end)=INDICES(end)-M;                        % revert last constraint to "+"
                    INDICES=[INDICES beg_mi(level)+counter(level)+M];   % add "-" constraint
                    break;
                else
                    counter(level)=0;
                    INDICES(m+sum(counter)+1:end)=[];
                    level=level-1;
                    % if we reach the root of the tree we are done
                    if level==0
                        return;
                    end
                end
            end
        end
        
    else
        % counter(level)>0
        % if we are not for the first time on this level
        %===============================================
        
        % increase counters
        %-----------------
        nzcount=find(counter>0); % find nonzero counters
        for jj=length(nzcount):-1:1
            level=nzcount(jj);
            counter(level)=counter(level)+1;
            if counter(level)<=mi(level)
                INDICES(end)=INDICES(end)-M;                        % revert last constraint to "+"
                INDICES=[INDICES beg_mi(level)+counter(level)+M];   % add "-" constraint
                break;
            else
                counter(level)=0;
                INDICES(m+sum(counter)+1:end)=[];
                level=level-1;
                % if we reach the root of the tree we are done
                if level==0
                    return;
                end
            end
        end
    end
    
    % is the problem feasible?
    %-------------------------
    [xopt,R]=chebyball_f(A(INDICES,:),B(INDICES),Options);
    if R>=abs_tol
        % if we are at the leaf level store the node
        if level==N
            if Options.reduce_output
                res = [res polytope(A(INDICES,:),B(INDICES))];
            else
                %%res = [res polytope(A(INDICES,:),B(INDICES),0,1)];
                res = Pdummy;
            end
            if ~Options.constructset
                return;                 % WE ARE NOT LOOKING FOR THE FULL SOLUTION -> END SEARCH
            end
        else
            level=level+1;              % if we are not at the leaf level, go further down the tree
        end
    end
end

return;
