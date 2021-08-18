function [V,R,P,adjV,adjf] = extreme(P,Options)
%EXTREME Calculates extreme points of a given polytope
%
% [V,R,P,adjV] = extreme(P,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes extreme points (vertices) of the given polytope
%
% NOTE:
%   Since computing extreme points is very expensive, there is a possibility to
%   store the extreme points in the internal POLYTOPE structure. If the vertices
%   are already stored in the polytope object, no computation will be performed!
%
%   That's why it is always recommended to use the function as follows:
%    [V,R,P]=extreme(P)
%
% USAGE:
%   V=extreme(P)             - returns extreme points
%   [V,R]=extreme(P)         - returns extreme points and rays
%   [V,R,P]=extreme(P)       - returns extreme points, rays, and the update
%                                    POLYTOPE structure
%   [V,R,P,adjV]=extreme(P)  - returns list of adjacent vertices
%
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                       - Polytope
% Options.extreme_solver  - Which method to use for vertex enumeration
%                           (0 - analytical enumeration, 3 - CDD)
%                           (see help mpt_init)
% Options.debug_level     - Sets the level of error checking 
%                           0: no additional checks performed
%                           1: check if the computed points are really extreme
%                           2: checks if hull of the computed vertices is
%                              identical to the initial polytope P
% Options.abs_tol         - absolute tolerance
% Options.roundat         - if CDD is used, it usually helps to round the
%                           H-representation of a polytope to certain number of
%                           decimal places. This option defines at which decimal
%                           point the representation should be rounded.
%                           (Default is Options.roundat=15, which means that the
%                           representation will be rounded to 15 decimal points)
%                           NOTE! rounding is only used if extreme_solver=3
%                           NOTE! set Options.roundat=Inf to disable rounding
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% V    - extreme points of the polytope
% R    - rays (if R is non-empty, the polytope is unbounded)
% P    - updated polytope object in which the vertices are now stored for faster
%          future computation
% adjV - indices of adjacent vertices
%
% see also HULL
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
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

% if ~isa(P,'polytope')
%     error('EXTREME: input argument MUST be a polytope');
% end

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<2,
    Options=[];
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end
if ~isfield(Options,'debug_level')
    Options.debug_level=mptOptions.debug_level;
end
if ~isfield(Options,'extreme_solver')
    Options.extreme_solver=mptOptions.extreme_solver;
end
if ~isfield(Options,'verbose')
    Options.verbose = mptOptions.verbose;
end
if ~isfield(Options, 'lpsolver'),
    Options.lpsolver = mptOptions.lpsolver;
end
if ~isfield(Options, 'roundat'),
    % to how many decimal points should the input points be rounded for CDD    
    Options.roundat = 15;
end

if length(P.Array)>0,
    error('EXTREME: this function does not work for arrays of polytopes!');
end

if ~P.minrep,
    P=reduce(P);   % remove redundant constraints
end

[nc,nx]=size(P.H);

V=[];
R=[];
adjV={};
adjf=[];
if ~isfulldim(P),
    return
end

Vstored = ~isempty(P.vertices);
Porig = P;

if Options.extreme_solver==3 || Options.extreme_solver==4,
    % cddmex
    if nargout>3,
        [V,R,P,adjV,adjf] = cddmex_extreme(P,Options);
    else
        [V,R,P] = cddmex_extreme(P,Options);
    end
    
elseif Options.extreme_solver==0,
    % analytical solution
    if nargout>3,
        [V,R,P,adjV,adjf] = matlab_extreme(P,Options);
    else
        [V,R,P] = matlab_extreme(P,Options);
    end
    
elseif Options.extreme_solver==1,
    if nargout>3,
        [V,R,P,adjV,adjf] = lrs_extreme(P,Options);
    else
        [V,R,P] = lrs_extreme(P,Options);
    end
    
elseif Options.extreme_solver==2,
    % analytic solution - alternative method
    if nargout>3,
        [V,R,P,adjV,adjf] = matlab_extreme_alt(P,Options);
    else
        [V,R,P] = matlab_extreme_alt(P,Options);
    end
    
else
    error('extreme: unknown solver!');
end


if(Options.debug_level>0) %& ~Vstored
    [result,i]=checkextreme(P,V,nx,Options);
    if result~=0,
        % something went wrong, try another solver
        cursol = mptOptions.extreme_solver;
        solvers = mptOptions.solvers.extreme;

        for isol = 1:length(solvers),
            nextsolver = solvers(isol);
            if nextsolver==cursol,
                % don't try current solver
                continue
            end
            cur_str = mpt_solverInfo('extreme', cursol);
            next_str = mpt_solverInfo('extreme', nextsolver);
            if Options.verbose>1,
                fprintf('EXTREME: "%s" failed, trying "%s"...\n', cur_str, next_str);
            end

            switch nextsolver
                case 3
                    % cddmex
                    if nargout>3,
                        [V,R,P,adjV,adjf] = cddmex_extreme(Porig,Options);
                    else
                        [V,R,P] = cddmex_extreme(Porig,Options);
                    end
                    
                case 1
                    % LRS - solve the problem twice to guarantee success
                    if nargout>3,
                        [V,R,P,adjV,adjf] = lrs_extreme(Porig,Options);
                    else
                        [V,R,P] = lrs_extreme(Porig,Options);
                    end
                    
                case 0
                    % matlab
                    if nargout>3,
                        [V,R,P,adjV,adjf] = matlab_extreme(Porig,Options);
                    else
                        [V,R,P] = matlab_extreme(Porig,Options);
                    end

                case 2
                    % matlab, alternative method
                    if nargout>3,
                        [V,R,P,adjV,adjf] = matlab_extreme_alt(Porig,Options);
                    else
                        [V,R,P] = matlab_extreme_alt(Porig,Options);
                    end
                    
                case 4
                    % do nothing, this solver is only utilized by hull.m
                otherwise
                    error(sprintf('EXTREME: unknown solver %d', nextsolver));
                    
            end
            cursol = nextsolver;
            [result,i]=checkextreme(P,V,nx,Options);
            if result==0,
                break
            end
        end
    end

    if result ~= 0,
        % something went wrong, try analytical computation
        V = sub_extreme_through_hull(P,Options);
        if isempty(V),
            P = polytope;
            nx = 0;
            adjV = {};
            adjf = [];
            R = [];
        else
            nx = size(V,2);
            P.vertices = V;
            adjf = [];
            adjV = {};
            R = [];
            if nargout > 3,
                [adjV,R] = sub_adjacent_vertices(V,P,Options);
            elseif nargout > 1,
                if ~isbounded(P),
                    R = 1;
                end
            end
        end
        
        [result,i]=checkextreme(P,V,nx,Options);
        
        if result==1,
            if(Options.extreme_solver~=0) && 0
                % this code does nothing, therefore we skip it...
                Options.extreme_solver=0;
                disp(['EXTREME: Point ' num2str(i) ' is not a vertex! Recomputing with brute force method...'])       
            else
                error(['EXTREME: Point ' num2str(i) ' is not a vertex!'])       
            end
            return
        elseif result==2,
            %is point inside polytope, i.e. on facet            
            if(Options.extreme_solver~=0)
                Options.extreme_solver=0;
                disp(['EXTREME: Point ' num2str(i) ' is not inside hull! Recomputing with brute force method...'])     
                if nargout<4,
                    [V,R,P] = matlab_extreme(P,Options);
                else
                    [V,R,P,adjV,adjf] = matlab_extreme(P,Options);
                end
            else
                error(['EXTREME: Point ' num2str(i) ' is not inside hull!'])     
            end
            return
        end
    end
end

if(Options.debug_level>1)
    %check if hull is equal
    if(hull(V,Options)~=P)
        if(Options.extreme_solver~=0)
            Options.extreme_solver=0;
            disp(['EXTREME: Hull of extreme vertices is not identical to original Polytope; Recomputing with brute force method...'])  
            [V,R,P,adjV,adjf] = matlab_extreme(P,Options);
        else
            error(['EXTREME: Hull of extreme vertices is not identical to original Polytope!'])  
        end
    end
end

% assign polytope with updated vertices in caller's workspace
if ~isempty(inputname(1)) && nargout==0,
    assignin('caller',inputname(1),P);
end



%-------------------------------------------------------------------------
function [V,R,P,adjV,adjf] = lrs_extreme(P,Options)

V = []; R = []; adjV={}; adjf=[];

[V,R] = mpt_lrs('extreme', P);
P.vertices = V;
adjf = [];
adjV = {};
if nargout > 3,
    [adjV] = sub_adjacent_vertices(V,P,Options);
end




%-------------------------------------------------------------------------
function [V,R,P,adjV,adjf] = cddmex_extreme(P,Options)

V = []; R = []; adjV={}; adjf=[];

if isempty(P.vertices) || nargout>=4,
    % if no vertices are stored in the internal structure, compute them
    % also perform the computation if adjancency information is requested
    
    % change almost-zero elements to zero
    P.H(abs(P.H)<1e-12) = 0;
    P.K(abs(P.K)<1e-12) = 0;
    
    if ~isinf(Options.roundat),
        % round H-representation to certain number of decimal places
        roundfactor = 10^Options.roundat;
        P.H = round(P.H*roundfactor) / roundfactor;
        P.K = round(P.K*roundfactor) / roundfactor;
    end
    
    tempH.A=P.H;
    tempH.B=P.K;
    failed=0; tempV=[]; adjV={};
    try
        if nargout>=4,
            % compute also adjancency information
            [tempV,adjV] = cddmex('adj_extreme',tempH);
        else
            tempV=cddmex('extreme',tempH);      % compute extreme points only
            adjV={};
        end
        V=tempV.V;
        R=tempV.R;
        adjf=[];
        P.vertices=V;
    catch
        failed=1;
    end
    if size(R,1)>0 || failed==1,
        % CDD sometimes fails to compute extreme points correctly
        % this happens often when slopes of two hyperplanes are too close
        % that's why we use a fixed-point arithmetics
        for ii=1:size(P.H,1),
            for jj=1:size(P.H,2),
                P.H(ii,jj)=fix(P.H(ii,jj)*(1/1e-5))*1e-5;
            end
            P.K(ii)=fix(P.K(ii)*(1/1e-5))*1e-5;
        end
        P.minrep=0;
        P=reduce(P);
        tempH.A=P.H;
        tempH.B=P.K;
        failed=0; tempV=[]; adjV={};
        try
            if nargout>=4,
                % compute also adjacency information
                [tempV,adjV] = cddmex('adj_extreme',tempH);
            else
                tempV=cddmex('extreme',tempH);
                adjV={};
            end
            V=tempV.V;
            R=tempV.R;
            P.vertices=V;
        end
    end
else
    R=[];
    V=P.vertices;     % vertices are already present, just return them
end



%-------------------------------------------------------------------------
function [V,R,P,adjV,adjf] = matlab_extreme(P,Options)

V = []; R = []; adjV={}; adjf=[];

wstatus = warning;
[nc,nx]=size(P.H);
if nx==1,
    % polytope is a line segment
    for ii=1:size(P.H,1),
        V=[V; P.K(ii)/P.H(ii)];
    end
    if size(P.H)==1,
        R=1;
    end
    adjV={};
    adjf=[];
    return

elseif nx==2,
    % 2D polytope
    if isempty(P.vertices)    % if no vertices are stored in the internal structure, calculate them
        H=P.H;
        K=P.K;
        alf=angle(H(:,1)+j*H(:,2));
        [Y,I]=sort(alf);                      % order hyperplanes in cyclic way
        H=[H; H(1,:)];
        K=[K; K(1)];
        I=[I; I(1)];
        X=[];
        for ii=1:nc
            HH=[H(I(ii),:); H(I(ii+1),:)];
            KK=[K(I(ii)); K(I(ii+1))];
            if isinf(cond(HH)),               % if matrix is singular, polytope is unbounded
                R=[R 1];
            else
                x=HH\KK;
                X=[X;x'];
            end
        end
    else
        X=P.vertices;
    end
    if nargout<4,
        V=X;
        adjf=[];
        adjV={};
    else
        H=P.H;
        K=P.K;
        for ii=1:size(X,1),
            vertex=X(ii,:)';                                % take ii'th extreme point
            vind=find(abs(H*vertex-K)<=Options.abs_tol);    % get indices of active constraints
            adjvertices=[];
            for jj=1:size(X,1),                             % go through all extreme points except of ii
                if ii==jj, continue, end
                ind=find(abs(H*vertex-H*X(jj,:)')<=Options.abs_tol);  % get indices of active constraints
                if length(intersect(vind,ind))>=1,          % if the two points have at least one common active constraint, they are adjacent
                    adjvertices=[adjvertices jj];
                end
            end
            if length(adjvertices)~=nx,                     % each point has to have 2 adjacent vertices, otherwise the polytope is not bounded!
                disp('warning: each vertex should have have 2 adjacent vertices!');
            end
            adjV{ii}=adjvertices;
        end


    end
elseif nx==3 && nchoosek(nc, nx)<=161700,
    % 3D polytope with less than 100 facets
    % if it has more than 100 facets, the number of all possible combinations
    % we have to explore is just so high that it is more efficient to use the
    % hull-approach via sub_extreme_through_hull()
    if isempty(P.vertices),
        H=P.H;
        K=P.K;
        X=[];
        adjf=[];

        warning off
        nadjv=0;
        adjW={};
        nc = size(P.H,1);
        nx = size(P.H,2);
        allcombs = nchoosek(1:nc,nx); %allcombs = combnk(1:nc,nx);
        for ii=1:length(allcombs),
            onecomb = allcombs(ii,:);
            HH = H(onecomb,:);
            KK = K(onecomb,:);
            x=HH\KK;                               % compute H representation
            hxk=H*x-K;
            if all(H*x-K<Options.abs_tol), % does the point lie in the original polytope?
                if sum(abs(hxk)<=Options.abs_tol)>=nx,
                    if size(X,1)==0,
                        X=[X; x'];
                    else
                        isthere=0;
                        for kk=1:size(X,1),
                            if all(abs(X(kk,:)-x')<=Options.abs_tol),
                                % was this point already stored?
                                isthere=1;
                                break
                            end
                        end
                        if ~isthere,
                            X=[X; x'];
                            adjf=[adjf; onecomb];
                            nadjv=nadjv+1;
                            adjW{nadjv}=find(abs(hxk)<=Options.abs_tol);
                        end
                    end
                end
            end
        end
    else
        X=P.vertices;
    end
    warning(wstatus);
    V = X;
    adjf = {};
    if nargout>3,
        [adjV, R] = sub_adjacent_vertices(X,P,Options);
    elseif nargout>1,
        if ~isbounded(P),
            R = 1;
        end
    end
        
else
    % general nD case
    X = sub_extreme_through_hull(P,Options);
    P.vertices = X;
    adjf = [];
    adjV = {};
    R = [];
    if isempty(X),
        V = [];
        P = polytope;
        return
    end
    if nargout > 3,
        [adjV,R] = sub_adjacent_vertices(X,P,Options);
    elseif nargout > 1,
        if ~isbounded(P),
            R = 1;
        end
    end
end


V=X;



%-------------------------------------------------------------------------
function [V,R,P,adjV,adjf] = matlab_extreme_alt(P,Options)
% alternative method for extreme points enumeration

V = []; R = []; adjV={}; adjf=[];

if ~isbounded(P),
    % this method requires a bounded polytope
    return
end

H = P.H;
K = P.K;
[xcheb, rcheb] = chebyball_f(H, K, Options);

% shift polytope such that it contains the origin
K = K - H*xcheb;
D = H ./ repmat(K,[1 size(H,2)]);
k = mpt_convhulln(D);
G  = zeros(size(k,1),size(D,2));
for ii = 1:size(k,1)
    F = D(k(ii, :), :);
    % create inequalities
    G(ii,:) = (F\ones(size(F,1),1))';
end
% shift polytope back to original location
V = G + repmat(xcheb', [size(G,1),1]);

if nargout>3,
    [adjV, R] = sub_adjacent_vertices(V,P,Options);
end


%-------------------------------------------------------------------------
function [result,i]=checkextreme(P,V,nx,Options)
% result = 0 - all ok
% result = 1 - point 'i' is not a vertex
% result = 2 - point 'i' is not inside of convex hull

result = 0;
i = 0;
if isempty(V),
    result = 1;
    return
end
for i=1:size(V,1)
    d = P.H*V(i,:)'-P.K;
    tmp=find(abs(d)<=Options.abs_tol); % find intersections with hyperplanes
    if length(tmp)<nx
        % check if each point is really a vertex
        result = 1;
        return
    elseif any(d>Options.abs_tol),
        % is point inside polytope, i.e. on facet
        result = 2;
        return
    end
end


%------------------------------------------------------------------------
function V = sub_extreme_through_hull(P,Options)

xmid = chebyball(P);
H = P.H;
K = P.K;
Ai = zeros(size(H));
for ii=1:size(H,1),
    Ai(ii,:) = H(ii,:)/(K(ii)-H(ii,:)*xmid);
end

Q = hull(Ai,Options);

if ~isfulldim(Q),
    V = [];
    return
end

H = Q.H;
K = Q.K;

V = H;

nx = size(V,2);
Vnew = zeros(size(V));
for iv = 1:size(V,1),
    for ix = 1:nx,
        Vnew(iv,ix) = H(iv,ix)/K(iv) + xmid(ix);
    end
end
V = Vnew;


%-----------------------------------------------------------------------
function [adjV,R] = sub_adjacent_vertices(X,P,Options),

H = P.H;
K = P.K;
abs_tol = Options.abs_tol;

adjV=[];
R = [];
nx = size(X,2);
for ii=1:size(X,1),
    vertex=X(ii,:)';
    vind=find(abs(H*vertex-K)<=abs_tol);
    adjvertices=[];
    for jj=1:size(X,1),
        if ii==jj, continue, end
        ind=find(abs(H*X(jj,:)'-K)<=abs_tol);
        if length(intersect(vind,ind))>=nx-1,
            adjvertices=[adjvertices jj];
        end
    end
    if length(adjvertices)<nx,
        %means that the polytope is unbounded
        R=[R 1];
    end
    adjV{ii}=adjvertices;
end
