function aux = auxgeometry3(node3,elem3)
%auxgeometry3 gets geometry data of polyhedral mesh
%
% Copyright (C) Terence Yu.

NT = size(elem3,1);
if ~iscell(elem3) % tetrahderal mesh: transform to cell
    elemTet = cell(NT,1);
    for iel = 1:NT
        tet = elem3(iel,:);
        face = [tet(:,[2 4 3]);tet(:,[1 3 4]);tet(:,[1 4 2]);tet(:,[1 2 3])];
        elemTet{iel} = mat2cell(face, ones(4,1), 3);  
    end
    elem3 = elemTet; 
end

centroid3 = zeros(NT,3); 
diameter3 = zeros(NT,1); volume = zeros(NT,1);
for iel = 1:NT
    % faces
    elemf = elem3{iel}; 
    % triangulation of faces
    [Tri,index3,TriLocal] = faceTriangulation(elemf);
    % centroid
    centroid3(iel,:) = polycentroid3(node3,Tri);
    % diameter
    V = node3(index3,:);
    diameter3(iel) = max(pdist(V));  
    % volume: divided into tetrahedrons by connecting centroid with face
    % triangulation
    if size(Tri,1)>4
        % triangulation of the polyhedron
        nodeTet = [V; centroid3];
        elemTet = [TriLocal, (max(TriLocal(:))+1)*ones(size(TriLocal,1),1)];
    else % tetrahedron
        nodeTet = V;
        elemTet = 1:4;
    end
    volume(iel) = sum(simplexvolume(nodeTet,elemTet));
end

aux.node3 = node3; aux.elem3 = elem3;
aux.centroid3 = centroid3;
aux.volume = volume;
aux.diameter3 = diameter3;