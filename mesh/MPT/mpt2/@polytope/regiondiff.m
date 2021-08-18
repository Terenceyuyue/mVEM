function [res]=regiondiff(P,Pn,Options)
%REGIONDIFF Region difference
%
% R=refiondiff(P,Pn,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Function computes the difference between polytope P
% and union of polytopes Pn:  R=P\(U Pn).
% In general output P is a polytope array containing
% non-overlapping polytopes whose union is equal to
% the difference P\(U Pn).
%
% USAGE:
%   R=regiondiff(P,Pn)
%   R=regiondiff(P,Pn,Options)
%
% Note:
%  Function is actually designed to answer the question
%   is P subset of (U Pn) ?
%  if answer is YES then P=empty polytope.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                     - polytope
% Pn                    - polytope array
% Options.abs_tol       - absolute tolerance
% Options.intersect_tol - tolerance to detect intersection
% Options.lpsolver      - Which LP solver to use (see help mpt_solveLP)
% Options.constructset  - if 1, construct the whole difference
% Options.reduce        - if 0, the polyarray Pn will not be reduced (for fast
%                         coverage detection)
% Options.reduce_output - if 0, the set P\Pn will not be reduced, i.e.
%                         redundant constraints will not be removed (for fast
%                         coverage detection)
% Options.presolved     - set to 1 if elements of Pn consist just of "active
%                         constraints". 0 otherwise
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R   - polytope array containing non-overlapping polytopes whose union is equal to P\Pn
%
% see also MLDIVIDE, REGIONDIFFXU
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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
    Options=[];
    if nargin<2,
        error('REGIONDIFF: Not enough input arguments!');
    end
end

if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;      % Use NAG Toolbox
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;   % Nonnegative tolerance to decide if a constraint is redundant
end
if ~isfield(Options,'intersect_tol'),
    Options.intersect_tol=sqrt(eps);
end
if ~isfield(Options,'reduce'),
    Options.reduce=1;         % 1 = reduces the description of input region R (NOT Ri!!!)
end
if ~isfield(Options,'reduce_output'),
    Options.reduce_output=1;  % 1 = reduce the description of output regions Pj={x | P.Hn{j} x <= P.Kn{j}}
end
if ~isfield(Options,'constructset'),
    Options.constructset=1;   % should we construct the difference, or just see if U Ri covers R
end
if ~isfield(Options,'presolved'),
    Options.presolved=0;      % have "active constraints" of Ri already been isolated
end
if ~isfield(Options,'infbox'),
    Options.infbox=mptOptions.infbox;
else
    Options.infbox=min(mptOptions.infbox,Options.infbox); % internal infbox cannot exceed MPT infbox
end

abs_tol = Options.abs_tol;
reduce_output = Options.reduce_output;
constructset = Options.constructset;

if ~isa(P,'polytope') | ~isa(Pn,'polytope'),
    error('REGIONDIFF: Input argument MUST be a polytope object!');
end

if ~isfulldim(P),
    res = mptOptions.emptypoly;
    return
end
if ~isfulldim(Pn),
    res = P;
    return
end

[cx,cr]=chebyball(P);
if isinf(cr) & all(cx==0),
    % P is R^n, convert it to infinity box polytope
    dimbox=dimension(Pn(1));
    P=polytope([eye(dimbox); -eye(dimbox)],ones(dimbox*2,1)*Options.infbox);
end


% initial value for the output
%-----------------------------
%emptypoly=polytope;
emptypoly = mptOptions.emptypoly;
res=emptypoly;

% find minimal representation for the region R
%---------------------------------------------
if Options.reduce,
    if ~isminrep(P),
        P=reduce(P);
    end
    if ~isfulldim(P)
        res=emptypoly;
        return
    end
    if ~isfulldim(Pn(1)),
        res=P;
        return
    end
end
    
Pdummy = P;

N=length(Pn);   % initial number of regions
if isempty(Pn.Array),
    Pn.Array{1}=Pn;
end

% check which Ri intersect R
%---------------------------
Rc=zeros(N,1);  % radius of Chebyshev ball

PH = P.H;
PK = P.K;
for ii=1:N
    PnAii = Pn.Array{ii};
    if ~isfulldim(PnAii),
        H=[];
        K=[];
    else
        H=PnAii.H;
        K=PnAii.K;
    end
    [xc,Rc(ii)]=chebyball_f([PH; H],[PK; K],Options);
end

[val,ind]=sort(-Rc);            % sort Ri in descending order

N=sum(Rc>=Options.intersect_tol);   % consider only Ri that do intersect R

% no Ri intersect R => return R
%------------------------------
if N==0
    res=P;
    return;
end

A=P.H;
B=P.K;
H=A;
K=B;              % storing all constraints in one matrix
m=size(A,1);      % number of rows in R
mi=zeros(N,1);    % beginning of the storage of the region Ri constraints

% Find "active constraints" of polyhedra Ri and put them all in one big matrix
% "active constraints" are those that do not show up in the description of R
%-----------------------------------------------------------------------------
if Options.presolved
    % "active constraints" separated in advance (description of Ri comprises only such constraints)
    for ii=1:N
        i=ind(ii);
        mi(ii)=size(Pn.Array{i}.H,1);
        A=[A; Pn.Array{i}.H];
        B=[B; Pn.Array{i}.K];
    end
else
    % find "active constraints"
    for ii=1:N
        i=ind(ii); 
        Hni=Pn.Array{i}.H;
        Kni=Pn.Array{i}.K;
        if ~isfulldim(Pn.Array{i}),
            Hni=[];
            Kni=[];
        end
        for j=1:size(Hni,1)
            % find unique facets (different from those defining R)
            HKnij = [Hni(j,:) Kni(j)];
            if all( sum( abs([H,K]-HKnij(ones(m,1),:)), 2) >= abs_tol )
                % replaces if all( sum( abs([H,K]-repmat([Hni(j,:) Kni(j)],m,1)), 2) >= abs_tol )
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


counter = zeros(1,N);               % constraints counters
INDICES = 1:m;                      % indices of currently valid constraints

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
            if reduce_output
                res = [res polytope(A(INDICES,:),B(INDICES))];
            else
                %%res = [res polytope(A(INDICES,:),B(INDICES),0,1)];
                res = Pdummy;
            end
            if ~constructset
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
    
    % is the problem feasible?
    %-------------------------
    [xopt,R]=chebyball_f(A(INDICES,:),B(INDICES),Options);
    if R>=abs_tol
        % if we are at the leaf level store the node
        if level==N
            if reduce_output
                res = [res polytope(A(INDICES,:),B(INDICES))];
            else
                %%res = [res polytope(A(INDICES,:),B(INDICES),0,1)];
                res = Pdummy;
            end
            if ~constructset
                return;                 % WE ARE NOT LOOKING FOR THE FULL SOLUTION -> END SEARCH
            end
        else
            level=level+1;              % if we are not at the leaf level, go further down the tree
        end
    end

end

return;
