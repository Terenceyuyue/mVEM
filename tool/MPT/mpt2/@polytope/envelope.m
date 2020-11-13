function R=envelope(P,Options)
%ENVELOPE Computes envelope of n polytopes
%
% R=envelope(P,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Compute envelope of polytopes
%
% USAGE:
%   R=envelope([P1 P2 P3]) - computes envelope of 3 polytope P1, P2 and P3
%   R=envelope(PA)         - if PA is a polytope array, computes envelope of the
%                            polytopes contained in that array.
%
% NOTE: the envelope may not exist! In that case the function returns R^n
%       (read the input section on how R^n is implemented)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                 - Polytope array
% Options.lpsolver  - LP solver to use (see help mpt_solveLP)
% Options.abs_tol   - absolute tolerance
% Options.infbox    - Internally, R^n is converted to a large box. This argument
%                       specifies bounds of this box.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R   - Polytope describing the envelope
%
% see also BOUNDING_BOX, UNION, HULL
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

% if ~isa(P, 'polytope')
%   error('ENVELOPE: Argument MUST be a polytope object');
% end

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<2
    Options=[];
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'lpsolver')
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'infbox'),
    Options.infbox=mptOptions.infbox;
end

nR=length(P.Array);
if nR==0,
    % if P is a polytope, just reduce the polytope
    if isfulldim(P,Options) & ~P.minrep,
        P=reduce(P);
    end
    R=P;
    return
end

nonvalid=cell(1,nR);
nNV=zeros(1,nR);
Ae=[];
Be=[];
nAi=zeros(1,nR);
nBi=zeros(1,nR);
nxi=zeros(1,nR);

for i=1:nR,
    if ~P.Array{i}.minrep,
        P.Array{i}=reduce(P.Array{i});   % reduce the polytope if it is not in minimal representation
    end
    nBi(i)=size(P.Array{i}.K,1);
    [nAi(i),nxi(i)]=size(P.Array{i}.H);
end

for region1=1:nR,    % go through all regions in P
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
            [x,R]=chebyball_f([P.Array{region2}.H; -P.Array{region1}.H(i,:)],[P.Array{region2}.K; -P.Array{region1}.K(i)],Options);
            if R>Options.abs_tol,
                i1(i)=0;
            end
        end
    end
    % positions of non-valid constraints
    nonvalid{region1}=find(~i1);
    nNV(region1)=length(nonvalid{region1});
    
    Ae=[Ae; P.Array{region1}.H(find(i1),:)];
    Be=[Be; P.Array{region1}.K(find(i1),1)];
end

if ~isempty(Ae) & ~isempty(Be),
    % envelope may be unbounded, in that case we intersect it with the infinity box
    dimbox = size(Ae, 2);
    Ae = [Ae; eye(dimbox); -eye(dimbox)];
    Be = [Be; ones(dimbox*2, 1) * Options.infbox];
    R=polytope(Ae,Be);
else
    % envelope is R^n, but keep dimension of the input polytope.
    nx = dimension(P);
    R=polytope(eye(nx), repmat(Inf, nx, 1), 1, 1, zeros(nx, 1), Inf);    
end
return;
