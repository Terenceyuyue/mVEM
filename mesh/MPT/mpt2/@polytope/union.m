function [Pu,how] = union(Pn, Options)
%UNION convex union computation
%
% [Pu,how]=union(Pn,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function tries to compute convex union of polytopes in the given polytope
% array Pn 
%
% USAGE:
%   U = union([P1 P2 P3 ...])
%   U = union([P1 P2 P3 ...], Options)
%   [U,how] = union([P1 P2 P3 ...], Options)
%
% NOTE! if Pn has more than 2 elements, the "convex union" problem is in fact
% solved with polytope/isconvex.m, unless Options.useisconvex=0 is specified.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn                  - Polytope array
% Options.useisconvex - if enabled, uses isconvex() to compute convex union.
%                     if number of regions is smaller than 3 and this flag is
%                     not defined, the default value will be 0. otherwise, it
%                     will be set to 1.
%                     
% Options.lpsolver  - solver used for LP problems
% Options.verbose   - level of verbosity
% Options.abs_tol   - absolute tolerance
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% Pu       - polytope describing convex union in case such one exists, return
%            array of input arguments otherwise
% how = 1    if union is convex
%     = 0    union is not convex (in this case Pu=Pn)
%
% see also ENVELOPE, HULL, ISCONVEX
%

% For the simple case of finding convex union of 2 polytopes, the following
% algorithm is used:
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
% Algorithm used is described in:
% Alberto Bemporad, Komei Fukuda, Fabio D. Torrisi,
% "Convexity recognition of the union of polyhedra",
% Computational Geometry, 18, pp. 141-154, 2001.
%
% If the input polyarray has more than two polytopes, we use the approach
% described in isconvex.m

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch
%


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

if ~isa(Pn,'polytope')
    error('UNION: you MUST pass polytopes to this function!');
end

emptypoly = mptOptions.emptypoly;

if ~isfulldim(Pn),  %check if polytope is empty
    Pu = emptypoly;
    how = 0;
    return
end

if isempty(Pn.Array),
    Pu=Pn;
    how=1;
    return
end

if ~isstruct(mptOptions),
    mpt_error;
    return
end

if nargin<2,
    Options=[];
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end

if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'reduce'),
    Options.reduce=0;           % Reduce regions to minimal representation
end
if ~isfield(Options,'verbose'),
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'useisconvex'),
    % if enabled, uses alternative approach to compute convex unions (faster on
    % higher dimensions)
    if length(Pn)<=2,
        % if we have only 2 regions, we switch to "standard" approach as
        % described in the corresponding literature (see above for citation)
        Options.useisconvex = 0;
    else
        Options.useisconvex = 1;
    end
end

nR=length(Pn.Array);

if Options.reduce,
    Pn=reduce(Pn);
end

if Options.useisconvex,
    % use isconvex() which is much faster on higher dimensions and/or with more
    % than 3 regions
    [how, Pu] = isconvex(Pn, Options);
    if ~how,
        Pu = Pn;
    end
    return
end

nAi=zeros(1,nR);
nBi=zeros(1,nR);
nxi=zeros(1,nR);


for i=1:nR,
    nBi(i)=size(Pn.Array{i}.K,1);
    [nAi(i),nxi(i)]=size(Pn.Array{i}.H);
end
if any(nBi-nAi),
    error('UNION: polytopes in the array must have the same number of constraints!');
end
if any(nxi-nxi(1)),
    error('UNION: polytopes in the array must have the same number of columns!');
end

% construct envelope by 
% removing non-valid constraints

nonvalid=cell(1,nR);
nNV=zeros(1,nR);

Hu=[];
Ku=[];
for region1=1:nR,

    i1=ones(1,nAi(region1));
    
    for region2=1:nR,
        if region2==region1,
            continue;
        end

        for i=1:nAi(region1),
            if ~i1(i),  % I have already discarded this constraint
                continue;
            end
            % solve an LP
            h = [Pn.Array{region2}.H; -Pn.Array{region1}.H(i,:)];
            k = [Pn.Array{region2}.K; -Pn.Array{region1}.K(i)];
            [x,R]=chebyball_f(h,k,Options);
            if R>Options.abs_tol,
                i1(i)=0;
            end
        end
    end

    % positions of non-valid constraints
    nonvalid{region1}=find(~i1);
    nNV(region1)=length(nonvalid{region1});
    
    Hu=[Hu; Pn.Array{region1}.H(find(i1),:)];
    Ku=[Ku; Pn.Array{region1}.K(find(i1),1)];
end
if(isempty(Hu))    %less constraints than simplex => unbounded
    how=0;
    Pu = Pn;
    return
end
Pu = polytope(Hu,Ku);   %envelope
[Hu,Ku] = double(Pu);
nHu=size(Hu,1);
if(nHu<nxi(1)+1)    %less constraints than simplex => unbounded
    how=0;
    %% Pu=polytope;
    Pu = Pn;
    return
end


how=1;

pos=ones(1,nR);
total=prod(nNV);

f=[zeros(1,nxi(1)) 1];

% check through all possible combinations of non-valid constraints
for i=1:total,
    Hmat=[];
    Kmat=[];

    for j=1:nR,
        Hmat=[Hmat; -Pn.Array{j}.H(nonvalid{j}(pos(j)),:) -1];
        Kmat=[Kmat; -Pn.Array{j}.K(nonvalid{j}(pos(j)),1)];
    end
    Hmat=[Hmat; Hu zeros(nHu,1)];
    Kmat=[Kmat; Ku];

    [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,Hmat,Kmat,[],[],[],Options.lpsolver);
    if strcmp(how,'unbounded') | (strcmp(how,'ok') & xopt(nxi(1)+1) < -Options.abs_tol),
        how=0;
        if Options.verbose>=2,
            disp('UNION warning: union of polytopes is not convex!!!');
        end
        Pu=Pn;
        return;
    end

    % adjust pos for next call
    for j=nR:-1:1,
        pos(j)=pos(j)+1;
        if pos(j)<=nNV(j),
            break;
        end
        pos(j)=1;
    end
end
how=1;
