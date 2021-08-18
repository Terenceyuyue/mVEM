function [edget,elem2edge,ep,et] = triaconn2(elem)
%TRIACONN2 edge-centred connectivity for a conforming 2-sim-
%plex triangulation.
%   [EE,TE,EP,ET] = TRIACONN2(TT) returns the edge-based ad-
%   jacency for a mesh of 2-simlexes (triangles). EE = [P1,
%   P2] is the set of unique 1-simplexes (edges) in the mesh
%   TT. TE = [E1,E2,E3] is the tria-to-edge adjacency, repr-
%   esenting the unique edges within each 2-simplex in TT.
%   [EP,ET] is the edge-to-tria adjacency, where the set of
%   2-simplexes adjacent to a given edge EI is ET(EP(EI,1):
%   EP(EI,2)). TT may be topologically non-manifold.

% 给定一条边序号 ei, ep(ei,1):ep(ei,2)

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 11/12/2014

NT = size(elem,1);
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[edget, ~, totalJ] = unique(totalEdge, 'rows');
elem2edge = reshape(totalJ,NT,3);

if (nargout <= 2), return; end

%------------------- edge-to-tria indexing: m trias per edge
NE = size(edget,1);
ep = zeros(NE,2);
%-------------------------------------- count trias per edge
ep(:,1) = +1;
ep(:,2) = accumarray(elem2edge(:),+1);

%-------------------------------------- init sparse indexing
ep(:,2) = cumsum(ep(:,2));
ep(2:NE,1) = ep(1:NE-1,2)+1 ;

%---------------------------- assemble edge-to-tria indexing
et = zeros(ep(NE,2),1); tp = ep(:,1);
for iel = 1 : NT
    %--------------------------------------------- adj. edge
    ei = elem2edge(iel,1);
    ej = elem2edge(iel,2);
    ek = elem2edge(iel,3);
    %--------------------------------------------- push tria
    et(tp(ei))=iel; tp(ei) = tp(ei)+1;
    et(tp(ej))=iel; tp(ej) = tp(ej)+1;
    et(tp(ek))=iel; tp(ek) = tp(ek)+1;    
end
end

