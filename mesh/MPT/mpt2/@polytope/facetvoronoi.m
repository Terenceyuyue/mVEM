function [Pvor, Find] = facetvoronoi(P,Idx)
%FACETVORONOI Computes an equivalent of voronoi diagrams for facets
%
% Pvor = facetvoronoi(P)
% Pvor = facetvoronoi(P,ind)
% [Pvor, Find] = facetvoronoi(P)
% [Pvor, Find] = facetvoronoi(P,ind)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Given a convex polytope P={x | H*x<=K}, the "voronoi diagram for facets of P"
% consists of convex sets "Pvor" for which it holds that any point inside of
% Pvor(i) is closer to face "i" of polytope "P" than to any other face.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P      - polytope object
% ind    - optional: to specify a subset of facets, i.e. Pvor(i) is the
%          set of points that are closer to facet ind(i) than to the other
%          facets in ind(:).
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Pvor   - voronoi partition of facets of P
% Find   - vector of the same length as "Pvor", the element "Find(i)" denotes to
%          which facet of "P" does "Pvor(i)" correspond to.
%
% see also MPT_VORONOI
%

% Copyright is with the following author(s):
%
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

global mptOptions

error(nargchk(1,2,nargin));

if ~isstruct(mptOptions),
    mpt_error;
end

if ~isempty(P.Array),
    error('Polytope arrays not supported by this function.');
end
if ~isfulldim(P),
    error('Input polytope must be fully dimensional.');
end
if ~isnormal(P),
    P = normalize(P);
end
if ~isminrep(P),
    P = reduce(P);
end

Pvor = mptOptions.emptypoly;
Find = [];
Hn = P.H;
Kn = P.K;
[nc, nx] = size(Hn);

if nargin<2 | isempty(Idx)
    Idx = 1:nc;
end



for ii = Idx,
    ind = setdiff(Idx, ii);

    hi = Hn(ii, :);
    ki = Kn(ii);
    Hj = Hn(ind, :);
    Kj = Kn(ind);
    
    % distance from face "i" is given by -(h_i*x - k_i)
    % distance from face "j" is given by -(h_j*x - k_j)
    % note that we take the "negative" distance, because a point "x" lies inside
    % of the convex polytope H*x <= K
    %
    % then the "facet voronoi diagram" for facet "i" can be computed as:
    %
    %  x \in P
    %  -(h_i*x - k_i) <= -(h_j*x - k_j)   \forall j \notequal i

    Q = polytope([Hn; -(repmat(hi, length(Idx)-1, 1)-Hj)], [Kn; Kj-ki]);
    if isfulldim(Q),
        Find = [Find ii];
    end
    Pvor = [Pvor Q];
end
