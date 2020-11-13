function [adjList] = mpt_buildAdjacency(Pi,Options)
%MPT_BUILDADJACENCY Builds adjacency information
%
% adjList = mpt_buildAdjacency(Pi)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% For a collection of non-overlaping polytopes Pi generates adjacency list.
% Fast and dirty implementation.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pi                 - polytope array (must be non-overlapping!)
% Options.abs_tol    - absolute tolerance (default is 1e-9)
% Options.remove_inf - if 1, forced cleaning of adjList data from numerical 
%                      inaccuracies, one can remove columns in neighbor list 
%                      that consist only of -Inf or 0. (default 0)
%                      [this can always be performed if Pi is convex]
%
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% adjList          - list of adjacent regions, where adjList{i}(j,:) is
%                    the list of neighbor regions indices at facet j of
%                    region i. "-Inf" denotes that the facet is at the
%                    boundary of feasibility
%                   
%

% Copyright is with the following author(s):
%
%(c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%         fjc@control.ee.ethz.ch
%(c) 2005 Miroslav Baric, Automatic Control Laboratory, ETH Zurich,
%         baric@control.ee.ethz.ch

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

% Motivated by: fjc, 2005.

if nargin < 2,
    Options = [];
end
if ( ~isfield(Options,'abs_tol') ),
    Options.abs_tol = 1e-9;
end
if ( ~isfield(Options,'remove_inf') ),
    Options.remove_inf = 0;
end
absTol   = min(Options.abs_tol,1e-9);   
stepSize = Options.abs_tol * 10;

if ~isa(Pi, 'polytope'),
    error('First input must be a polytope object.');
end
nPoly    = length(Pi);
adjList  = cell(nPoly,1);
dimPoly  = dimension(Pi(1));

for polyIdx = 1:nPoly,
   [Hpoly,Kpoly] = double(Pi(polyIdx));
   [xCheby,rCheby] = chebyball(Pi(polyIdx));
   nFacets = size(Hpoly,1);
   
   for faceIdx = 1:nFacets,
       % find the facet representation
       %
       hFacet = Hpoly(faceIdx,:); % this is normalized vector
       kFacet = Kpoly(faceIdx);
       xFacet = xCheby + (kFacet - hFacet*xCheby) * hFacet';
       N0 = null(hFacet);
       H0 = Hpoly * N0;
       K0 = Kpoly - Hpoly * xFacet;
       Pfacet = polytope(H0,K0);  % this is facet representation
       
       % get the list of previously discovered neighbors and substract it
       % from the facet
       %
       neighList = adjList{polyIdx};
       faceList = [];
       if ( size(neighList,1)>=faceIdx ),
           faceList = neighList(faceIdx,:);
           faceList = faceList(find(faceList > 0)); % remove free slots
           for nzIdx = 1:length(faceList),
               % ok, this is MESSY: find the neighboring facet, project
               % it to the facet we're crossing and find the remaining
               % part
               %
               [Hneigh,Kneigh] = double(Pi(faceList(nzIdx)));
               Kproj = Kneigh - Hneigh * xFacet;
               Hproj = Hneigh * N0;
               Pproj = polytope(Hproj,Kproj);
               
               lenPfacet = length(Pfacet);
               for ii = 1:lenPfacet,
                   auxPoly = Pfacet(ii) \ Pproj;
                   Pfacet(ii)  = auxPoly(1);
                   for jj = 2:length(auxPoly),
                       Pfacet(end+1) = auxPoly(jj);
                   end
               end
           end
       end

       for extraIdx = 1:length(Pfacet),
           isBorder = false;
           while ( isfulldim(Pfacet(extraIdx)) ),
               % find the point on what remaining of the original facet and cross
               % the facet at that point
               %
               [xStart,dummyR] = chebyball(Pfacet(extraIdx));
               % perturb the initial point a little bit
               %
               randPoint = randn(dimPoly-1,1);
               pertVect = randPoint - xStart;
               pertVect = pertVect / norm(pertVect);
               xStart = xStart + pertVect * dummyR / 2;
               xStart = xFacet + N0 * xStart;   % lift xStart
               xBeyond = xStart + stepSize * hFacet';

               % find a new neighbor
               %
               [isin,inwhich,closest] = isinside(Pi,xBeyond,Options);
               if ( ~isin ),
                   faceList(end+1) = -Inf; % edge of the feasible area
                   isBorder = true;
               else
                   % we should get only one new neighbor. Since we're doing a
                   % random perturbation of the crossing point it is HIGHLY
                   % unlikely thet we hit the border of two regions on the other
                   % side. So, no additional check. Also, no check if the "new"
                   % neighbor is actually the old one. This shouldn't happen.
                   %
                   newNeighIdx = inwhich(1);
                   faceList(end+1) = newNeighIdx;

                   % update adjacency information for the new neighbor
                   %
                   otherList = adjList{newNeighIdx};
                   % find the neighboring facet
                   %
                   [Hneigh,Kneigh] = double(Pi(newNeighIdx));
                   neighSlacks = Kneigh - Hneigh * xStart;
                   [minSlack,otherFaceIdx] = min(abs(neighSlacks));
                   if ( size(otherList,1)<otherFaceIdx )
                       otherList(otherFaceIdx,1) = polyIdx;
                   else
                       freeIdx = find(otherList(otherFaceIdx,:) == 0);
                       if ( ~isempty(freeIdx) ),
                           otherList(otherFaceIdx,freeIdx(1)) = polyIdx;
                       else
                           otherList(otherFaceIdx,end+1) = polyIdx;
                       end
                   end
                   adjList{newNeighIdx} = otherList;
               end
               % and .. that's it, just store the updated adajcency info into the
               % big list and continue
               %
               adjList{polyIdx}(faceIdx,1:length(faceList)) = faceList;
               if ( isBorder ),
                   break;
               end

               % take the covered part off and proceed
               %
               [Hneigh,Kneigh] = double(Pi(newNeighIdx));
               Kproj = Kneigh - Hneigh * xFacet;
               Hproj = Hneigh * N0;
               Pproj = polytope(Hproj,Kproj);
               lenPfacet = length(Pfacet);
               auxPoly = Pfacet(extraIdx) \ Pproj;
               Pfacet(extraIdx) = auxPoly(1);
               for jj = 2:length(auxPoly),
                   Pfacet(end+1) = auxPoly(jj);
               end
           end  % end WHILE
       end % end FOR extraIdx
   end  % end FOR (facets)
end  % end FOR (polytopes)


% clean neighborlist data from numerical inaccuracies
% if Pi is convex, one can remove columns in neighbor list that consist only of
% -Inf or 0
if Options.remove_inf
    status_convex = 1;
else
    [status_convex] = isconvex(Pi);
end

if status_convex
    for kk=1:length(adjList)
        for ll=1:length(adjList{kk}(1,:))
            if all(adjList{kk}(:,ll)<=0)
                adjList{kk}(:,ll) = [];
            end
        end
    end
end
 
return
