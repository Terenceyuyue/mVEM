function [P,Vconv]=hull(V,Options)
%HULL Converts vertices/polytope array into an H-representation polytope
%
% [P,Vconv]=hull(V,Options)
% [P,Vconv]=hull(PA,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Creates convex hull of vertices or of a polytope array and returns the 
% H-representation of the hull
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% V                      - matrix containing vertices of the polytope
% PA                     - polytope array
% Options.extreme_solver - which method to use for convex hull computation
%                          (see help mpt_init for details)
% Options.abs_tol        - absolute tolerance
% Options.roundat        - if CDD is used, it usually helps to round the
%                          V-representation of a polytope to certain number of
%                          decimal places. This option defines at which decimal
%                          point the representation should be rounded.
%                          (Default is Options.roundat=14, which means that the
%                          representation will be rounded to 14 decimal points)
%                          NOTE! rounding is only used if extreme_solver=3
%                          NOTE! set Options.roundat=Inf to disable rounding
%
% Note: Initial values of the Options variable are given by mptOptions 
%       (see help mpt_init)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% P      - an H-representation polytope P={x | H x <= K}
% Vconv  - extreme points (i.e. points forming the convex hull)
%
% see also POLYTOPE/HULL, EXTREME, UNION, ENVELOPE

% Copyright is with the following author(s):
%
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich
%     fjc@control.ee.ethz.ch
% (C) 2005 Mario Vasak, FER, Zagreb
%     mario.vasak@fer.hr
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich
%     kvasnica@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich
%     kvasnica@control.ee.ethz.ch

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

% if no options given, use default values
if ~isfield(Options,'extreme_solver')
    Options.extreme_solver=mptOptions.extreme_solver;
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'infbox')
    Options.infbox=mptOptions.infbox;
end
if ~isfield(Options,'verbose')
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'lpsolver')
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'debug_level')
    Options.debug_level=mptOptions.debug_level;
end
if ~isfield(Options,'noReduce')
    Options.noReduce = 0;
end
if ~isfield(Options,'novertred')
    % if set to 1, vertices will not be reduced to kick out non-extreme points
    % before calling cddmex('hull', V) - this may help in certain tricky cases
    Options.novertred = 0;
end
if ~isfield(Options, 'roundat'),
    % to how many decimal points should the input points be rounded for CDD
    Options.roundat = 14;
end

P = mptOptions.emptypoly;
Options.emptypoly = P;

if isempty(V)
    Vconv=[];
    return
end

Vorig = V;
nx = size(V, 2);

% initiliaze variables ([issue91])
Vconv = []; H = []; K = [];

if Options.extreme_solver==3,
    % cddmex
    [P, Vconv, H, K, lowdim] = hull_cddmex(V, Options);
    if lowdim,
        return
    end
    
elseif Options.extreme_solver==1,
    % LRS (matlab implementation)
    [P, Vconv, H, K] = hull_lrs(V, Options);

elseif Options.extreme_solver==2,
    % analytic solution, alternative method
    [P, Vconv, H, K] = hull_matlab_alt(V, Options);
    
else
    % analytical solution
    [P, Vconv, H, K] = hull_matlab(V, Options);
end

if Options.debug_level>0
    if size(H, 2) == nx,
        % do not check hull if H is empty, otherwise an error about
        % incompatible dimensions
        % [issue91]
        [result,i] = checkhull(P,H,K,Vconv,Options);
    else
        result = 1;
        i = 1;
    end
    if result~=0,
        % something went wrong, try another solver
        cursol = mptOptions.extreme_solver;
        
        % introduce new solver (4) - hull_cddmex called with Options.novertred=1
        % we only include it if extreme_solver=3 is available
        solvers = mptOptions.solvers.extreme;
        solvers(find(solvers==Options.extreme_solver)) = [];
        
        for isol = 1:length(solvers),
            nextsolver = solvers(isol);
            if nextsolver==cursol,
                % don't try current solver
                continue
            end
            cur_str = mpt_solverInfo('extreme', cursol);
            next_str = mpt_solverInfo('extreme', nextsolver);
            if Options.verbose>1,
                fprintf('HULL: "%s" failed, trying "%s"...\n', cur_str, next_str);
            end
            
            switch nextsolver
                case 4
                    % cddmex - hull based on full vertices
                    cddopt = Options;
                    cddopt.novertred = 1;
                    [P, Vconv, H, K] = hull_cddmex(Vorig, cddopt);
                case 3
                    % cddmex
                    cddopt = Options;
                    cddopt.novertred = 0;
                    [P, Vconv, H, K] = hull_cddmex(Vorig, cddopt);
                case 1
                    % LRS
                    [P, Vconv, H, K] = hull_lrs(Vorig, Options);
                case 0
                    % matlab
                    [P, Vconv, H, K] = hull_matlab(Vorig, Options);
                case 2
                    % matlab alternative
                    [P, Vconv, H, K] = hull_matlab_alt(Vorig, Options);
                otherwise
                    error(sprintf('HULL: unknown solver %d', nextsolver));
            end
            cursol = nextsolver;
            if ~isempty(H),
                [result,i] = checkhull(P,H,K,Vconv,Options);
            else
                result = 1;
                i = 1;
            end
            if result==0,
                break
            end
        end
    end
    
    if result~=0,
        % still something wrong, try another method...
        
        % use approach based on extreme points enumeration
        % as described in:
        %
        % A Problem in Enumerating Extreme Points
        % Katta G. Murty
        % Dep. of Industrial and Operations Engineering
        % University of Michigan
        % http://www-personal.engin.umich.edu/~murty/segments-10.pdf

        disp('HULL: problem detected, using alternative method...');
        V = unique(Vorig, 'rows');     % first kick out identical points
        Vconv=mpt_convhulln(V);        % identify points forming the convex hull
        Vconv=V(unique(Vconv),:);      % store extreme points
        [nvert,ndim] = size(Vconv);
        xmid = (sum(Vconv)/nvert)';    % compute interior point
        Yt = zeros(nvert,ndim);
        for iv = 1:nvert,              % make transformation y = x - xmid
            for id = 1:ndim,
                Yt(iv,id) = Vconv(iv,id) - xmid(id);
            end
        end
        PP = polytope(Yt,ones(nvert,1),1); % form polytope in q-space {q : q y <= 1}
        V = extreme(PP,Options);           % enumerate it's extreme points
        nvert = size(V,1);
        K = zeros(nvert,1);
        H = V;
        for iv = 1:nvert,
            K(iv) = 1 + V(iv,:)*xmid;
        end
        P = polytope(H,K,0,2);               % convex hull is given by {x : q_i x <= 1 + q_i xmid}
        %[H,K] = double(P);
        % Reason why the line above is commented out (by Mario Vasak):
        %
        % the command P=polytope(H,K,0,2) just a few rows back could return an empty
        % polytope if cheby radius of P falls below the abs_tol, and that is 
        % correct. But taking H & K in row 490 then from P brings empty H and K in 
        % the superimposed function, and this ends in a hull-error (result=1), but 
        % nothing is really wrong: H and K are hopefully correctly computed, but 
        % the resulting hull is too small, so H and K should be checked in 
        % checkhull.m, while hull.m should return an empty polytope if the check 
        % is all right.
        
        [result,i] = checkhull(P,H,K,Vconv,Options);
        if result~=0,
            error(['HULL: Point ' num2str(i) ' is not a vertex!'])
            return
        end
    end
end

if Options.noReduce==0,
    % reduce the polytope representation if required
    if ~isminrep(P),
        P = reduce(P);
    end
end

%--------------------------------------------------------------------------
function [P,Vconv,H,K,lowdim]=hull_cddmex(V, Options)

% change almost-zero elements to zero
V(abs(V)<1e-12) = 0;

if ~isinf(Options.roundat),
    % round V-representation to certain number of decimal places
    roundfactor = 10^Options.roundat;
    V = round(V*roundfactor) / roundfactor;
end

Vorig = V;
vert.V=V;

% initiliaze variables ([issue91])
P=[]; Vconv = []; H = []; K = []; lowdim = 0;

if Options.novertred==0,
    try
        vert=cddmex('reduce_v',vert);  % kick out points which are not extreme
        %vert=cddmex('v_hull_extreme',vert);  % kick out points which are not extreme
    catch
        % return on error, other solver will be tried ([issue91])
        return
    end
    if isempty(vert.V),
        Vconv=[];
        P = Options.emptypoly;
        return
    end
    if size(vert.R,1)>0,
        % rays detected, hull of points is not bounded
        % in that case we make an intersection with the infinity box to guarantee boundedness
        if Options.verbose>0,
            disp('HULL(V): polytope is unbounded! making an intersection with infbox...')
        end
        Vinf=extreme(unitbox(size(V,2),Options.infbox),Options);   % extreme points of the unit hypercube
        vert_inf.V=[V; Vinf];                                      % make the intersection
        
        try
            vert=cddmex('v_hull_extreme',vert_inf);                    % kick out points which are not extreme
        catch
            % return on error, other solver will be tried ([issue91])
            P=[]; Vconv = []; H = []; K = [];
            return
        end
        if size(vert.R,1)>0,
            error('HULL(V): polytope still unbounded after correction!');
        end
    end
    Vconv = vert.V;
end

wstatus = warning;
warning off
lowdim = 0;

if Options.novertred,
    vert2.V = Vorig;
else
    vert2.V = Vconv;
end
try
    HK=cddmex('hull',vert2);        % compute H-representation convex hull
    if ~isempty(HK.lin),
        lowdim = 1;
        HK = [];
    end
catch
    HK = [];
end
warning(wstatus);

if isempty(HK),
    % return on error, other solver will be tried ([issue91])
    P=Options.emptypoly; Vconv = []; H = []; K = [];
    return
end

H=HK.A; %needed for comparison/debug later
K=HK.B; %needed for comparison/debug later
if isempty(H),
    P = Options.emptypoly;
    Vconv = [];
else
    % do not reduce the polytope yet, it will be done at the end
    P=polytope(H,K,0,2);
    P = set(P,'vertices',Vconv);
    %[H,K] = double(P);
    % Reason why the line above is commented out (by Mario Vasak):
    %
    % the command P=polytope(H,K,0,2) just a few rows back could return an empty
    % polytope if cheby radius of P falls below the abs_tol, and that is 
    % correct. But taking H & K in row 490 then from P brings empty H and K in 
    % the superimposed function, and this ends in a hull-error (result=1), but 
    % nothing is really wrong: H and K are hopefully correctly computed, but 
    % the resulting hull is too small, so H and K should be checked in 
    % checkhull.m, while hull.m should return an empty polytope if the check 
    % is all right.
end
return


%--------------------------------------------------------------------------
function [P,Vconv,H,K]=hull_matlab(V, Options)

V = unique(V, 'rows');             % first kick out identical points
if size(V,2)==3,                   % points in 3D space
    Vconv=mpt_convhulln(V);        % identify points forming the convex hull
    H=[];                          % initialize the matrices
    K=[];
    for ii=1:size(Vconv,1),
        r=Vconv(ii,:);             % pick up points which define extreme hyperplane
        p1=V(r(1),:)';
        p2=V(r(2),:)';
        p3=V(r(3),:)';
        h=-cross(p1-p2,p1-p3)';    % do the cross-product to get hyperplane description
        k=h*p1;
        H=[H; h];
        K=[K; k];
    end
    if isempty(H),
        P = Options.emptypoly;
    else
        % do not reduce the polytope yet, it will be done at the end
        P=polytope(H,K,0,2);
        Vconv=V(unique(Vconv),:);      % store extreme points
        P=set(P,'vertices',Vconv);
        %[H,K] = double(P);
        % Reason why the line above is commented out (by Mario Vasak):
        %
        % the command P=polytope(H,K,0,2) just a few rows back could return an empty
        % polytope if cheby radius of P falls below the abs_tol, and that is 
        % correct. But taking H & K in row 490 then from P brings empty H and K in 
        % the superimposed function, and this ends in a hull-error (result=1), but 
        % nothing is really wrong: H and K are hopefully correctly computed, but 
        % the resulting hull is too small, so H and K should be checked in 
        % checkhull.m, while hull.m should return an empty polytope if the check 
        % is all right.
    end
    return

elseif size(V,2)==2,               % points in 2D space
    Vconv=mpt_convhulln(V);        % identify points forming the convex hull
    puni=V(unique(Vconv),:);       % pick up extreme points

    % create a point in the middle by taking a convex combination of all extreme points
    pp=[1/length(puni)*sum(puni(:,1)); 1/length(puni)*sum(puni(:,2))];
    H=[];
    K=[];
    for ii=1:size(Vconv,1),
        r=Vconv(ii,:);
        p1=V(r(1),:);              % pick up 1st point
        p2=V(r(2),:);              % pick up 2nd point
        if (p1(1)==p2(1)),         % special case, first coordinates are identical
            h=[1 0];
            k=p1(1);
        elseif (p1(2)==p2(2)),     % special case, second coordinates are identical
            h=[0 1];
            k=p1(2);
        else
            m=(p2(2)-p1(2))/(p2(1)-p1(1));  % general case, compute the slope
            b=p1(2)-m*p1(1);                % and an offset
            h=[m -1];
            k=-b;
        end
        if h*pp-k>Options.abs_tol,   % identify proper sign of the hyperplane
            h=-h;
            k=-k;
        end
        H=[H; h];
        K=[K; k];
    end
    if isempty(H),
        P = Options.emptypoly;
    else
        % do not reduce the polytope yet, it will be done at the end
        P = polytope(H,K,0,2);
        Vconv=puni;                     % store extreme points
        P=set(P,'vertices',Vconv);
        %[H,K] = double(P);
        % Reason why the line above is commented out (by Mario Vasak):
        %
        % the command P=polytope(H,K,0,2) just a few rows back could return an empty
        % polytope if cheby radius of P falls below the abs_tol, and that is 
        % correct. But taking H & K in row 490 then from P brings empty H and K in 
        % the superimposed function, and this ends in a hull-error (result=1), but 
        % nothing is really wrong: H and K are hopefully correctly computed, but 
        % the resulting hull is too small, so H and K should be checked in 
        % checkhull.m, while hull.m should return an empty polytope if the check 
        % is all right.
    end
    return
else
    % general n-D case, version by Mario Vasak, FER, 20th July, 2005
    
    V_all=V;                
    
    Vconv = mpt_convhulln(V);         
    V=V(unique(Vconv),:);       % pick up only extreme points
    Vconv=mpt_convhulln(V);     % take the combinations of extreme points only
    nx = size(Vconv,2);           % dimension
    nv = size(V,1);               % number of extreme points
    I=1:nv;                       % I is the set of all extreme points
    H=[];
    K=[];
    H_my=[];
    K_my=[];
    combinations_taken=zeros(1,size(Vconv,1));      %used to find which points are really extreme
    not_encircled_points=zeros(1,size(V_all,1));    %for debugging, not yet used, intended to check whether ALL the points are on the one side of the facet
    badly_computed_q={};                            %for debugging, not yet used, intended to show which of the taken facets do not contain all the points on one side
    for k=1:size(V_all,1)
        badly_computed_q{k}=[];        
    end
    max_violation=0;                                %to find the greatest violation of constraints among all of the points
    for q=1:size(Vconv,1), 
        J=Vconv(q,:);           % J is a set of points lying on the same extreme hyperplane (facet)
        ImJ = setdiff(I,J);     % get indices of the set I-J
        AAA=[V(J,:) -ones(length(J),1)];
        nnn=null(AAA);          %find the null-space of AAA
        
        if size(nnn,2)==1       %if the null-space is one-dimensional (non-degenerated facet)....
            V_ImJ=V(ImJ,:);
        
            row_norm=norm(nnn(1:end-1));
            test=(V_ImJ*nnn(1:end-1)-nnn(end))/row_norm;
            if all(test<=Options.abs_tol)                %...and all other extreme points lie strictly on one side of the hyperplane...
                H_my=[H_my; (nnn(1:end-1))'/row_norm];      %...take that facet
                K_my=[K_my; nnn(end)/row_norm];
                combinations_taken(q)=1;
                test_2=H_my(end,:)*V_all'-K_my(end);
                affected_indices=find(test_2>=Options.abs_tol);
                if ~isempty(affected_indices)
                    max_violation=max(max_violation,max(test_2));
                end
                not_encircled_points(affected_indices)=1;
                for k=1:length(affected_indices)
                    badly_computed_q{affected_indices(k)}(end+1)=q;
                end
            elseif all(test>=-Options.abs_tol)           %...and all other extreme points lie strictly on one side of the hyperplane...
                H_my=[H_my; (-nnn(1:end-1))'/row_norm];     %...take that facet
                K_my=[K_my; -nnn(end)/row_norm];
                combinations_taken(q)=1;
                test_2=H_my(end,:)*V_all'-K_my(end);
                affected_indices=find(test_2>=Options.abs_tol);
                if ~isempty(affected_indices)
                    max_violation=max(max_violation,max(test_2));
                end
                not_encircled_points(affected_indices)=1;
                for k=1:length(affected_indices)
                    badly_computed_q{affected_indices(k)}(end+1)=q;
                end
            end
        end
    end

    if Options.verbose > 0,
        if ~isempty(find(not_encircled_points)) %display warning if some of the initially given points is left outside - 3rd change-for-Michal
            disp(['Warning: Maximal occurring violation in H-rep of convex hull is ',num2str(max_violation/Options.abs_tol),' times the absolute tolerance.']);
        end
    end
    
    H=H_my;
    K=K_my;
    
    if(~isempty(H))
        % do not reduce the polytope yet, it will be done at the end
        P = polytope(H,K,0,2);
        Vconv=V(unique(Vconv(find(combinations_taken),:)),:);  %so, extreme points are only those combinations from Vconv recognized as facets
        P=set(P,'vertices',Vconv);
    else
        P=Options.emptypoly;
        Vconv=[];
    end
    %[H,K] = double(P);
    % Reason why the line above is commented out (by Mario Vasak):
    %
    % the command P=polytope(H,K,0,2) just a few rows back could return an empty
    % polytope if cheby radius of P falls below the abs_tol, and that is 
    % correct. But taking H & K in row 490 then from P brings empty H and K in 
    % the superimposed function, and this ends in a hull-error (result=1), but 
    % nothing is really wrong: H and K are hopefully correctly computed, but 
    % the resulting hull is too small, so H and K should be checked in 
    % checkhull.m, while hull.m should return an empty polytope if the check 
    % is all right.
end
return


%--------------------------------------------------------------------------
function [P,Vconv,H,K]=hull_matlab_alt(V, Options)

w = warning;
warning off
V = unique(V, 'rows');
k = mpt_convhulln(V);
Vconv = V(unique(k),:);  % store extreme points
c = mean(Vconv);
V = V - repmat(c,[size(V,1) 1]);
H = zeros(size(k, 1), size(V, 2));
ind_remove = [];
for ix = 1:size(k,1)
    F = V(k(ix,:),:);
    aa = (F\ones(size(F,1),1))';
    if any(isinf(aa)),
        ind_remove = [ind_remove ix];
        continue
    end
    H(ix,:) = aa;
end
ind_keep = setdiff(1:size(k, 1), ind_remove);
H = H(ind_keep, :);
K = ones(size(H, 1), 1);
K = K + H*c';
if isempty(H),
    P = Options.emptypoly;
else
    % do not reduce the polytope yet, it will be done at the end
    P = polytope(H, K, 0, 2);
    %[H,K] = double(P);
    % Reason why the line above is commented out (by Mario Vasak):
    %
    % the command P=polytope(H,K,0,2) just a few rows back could return an empty
    % polytope if cheby radius of P falls below the abs_tol, and that is 
    % correct. But taking H & K in row 490 then from P brings empty H and K in 
    % the superimposed function, and this ends in a hull-error (result=1), but 
    % nothing is really wrong: H and K are hopefully correctly computed, but 
    % the resulting hull is too small, so H and K should be checked in 
    % checkhull.m, while hull.m should return an empty polytope if the check 
    % is all right.
end
warning(w);



%--------------------------------------------------------------------------
function [P,Vconv,H,K]=hull_lrs(V, Options)

V = unique(V, 'rows');     % kick out identical points
Vconv=mpt_convhulln(V);    % identify points forming the convex hull
Vconv=V(unique(Vconv),:);  % store extreme points
[P,dummy,H,K] = mpt_lrs('hull',Vconv,Options);


%--------------------------------------------------------------------------
function [result,i]=checkhull(P,H,K,Vconv,Options)
% result = 0 - all ok
% result = 1 - point 'i' is not a vertex
% result = 2 - point 'i' is not inside of convex hull

abs_tol = Options.abs_tol;
result = 0;
i = 0;
if isempty(Vconv) & ~isempty(H),
    % something is definitelly wrong. we don't have any vertices, 
    % even though H-representation is not empty
    result = 1;
    return
end

nx = size(Vconv,2);
for i=1:size(Vconv,1)
    HVK = H*Vconv(i,:)'-K;
    tmp=find(abs(HVK) <= abs_tol); %find intersections with hyperplanes
    if length(tmp) < nx
        % point is not a vertex
        result = 1;
        return
    elseif (any(HVK > abs_tol))
        % point is NOT inside of the polytope
        result = 2;
        return
    end
end
