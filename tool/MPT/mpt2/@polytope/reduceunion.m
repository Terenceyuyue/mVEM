function [Pret, keep] = reduceunion(P,Options)
% REDUCEUNION Removes redundant elements from a polytope array
%
% [Pret,kept] = reduceunion(P,Options),
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% given a (possibly non-convex) union of polytopes P, removes
% all elements of P which are completely covered by the other regions
%
% i.e. assume we have 3 polytopes:
%      p1 = polytope([0 -1; -1 0; 0 1; 1 0],[1;0;1;1])
%      p2 = polytope([0 -1; -1 0; 0 1; 1 0],[1;1;1;0])
%      p3 = polytope([0 -1; -1 0; 0 1; 1 0],[1;1;1;1])
%
% then if Pu=[p1 p2 p3], this function removes polytopes p1 and p2, since they are completely
% covered by a larger polytope (p3 in this case)
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%   P           -   polytope array
%   Options     -   will be used when calling regiondiff
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%
%  Pret    -  reduced polytope array
%  kept    -  0/1 vector of length(P) which stores a 0 at index i if polytope
%             P(i) is redundant and 1 at index i if P(i) is not redundant.
%

% (C) 2007 Michal Kvasnica, Slovak University of Technology in Bratislava
%          michal.kvasnica@stuba.sk
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

if nargin < 2
    Options = [];
end
Options = [];

lenP=length(P.Array);
if lenP == 0
    Pret = P;
    regions_kept = 1;
    keep = 1;
    return 
end

% options for fast subset check
Options = mpt_defaultOptions(Options, ...
    'reduce', 0, ...
    'constructset', 0, ...
    'reduce_output', 0 );

% sort regions by size for better performance.
% reasoning: at each step we need to check whether polytope P(i) is covered
%            by all other remaining polytopes. it is more likely that small
%            regions will be kicked out at the beginning, reducing the size
%            of the polytope array we need to check P(i) against.
[xcheb, rcheb] = chebyball(P);
[rcheb, index] = sort(rcheb); 
keep=ones(1, lenP);

index = index(:);
for selected = index',
    regions_kept = find(keep);
    to_check = setdiff(regions_kept, selected);
    selected_region = P.Array{selected};
    if isempty(to_check)
        % we kicked out all other regions, return the current one
        Pret = selected_region;
        return
    end
    Q = P; Q.Array = P.Array(to_check);
    R = regiondiff(selected_region, Q, Options);
    if ~isfulldim(R)
        % region "selected" is covered by the remaining regions
        keep(selected) = 0;
    end
end

%write ouput
regions_kept = find(keep);
if isempty(regions_kept)
    Pret=polytope;
else
    Pret = polytope(P.Array(regions_kept));
end
