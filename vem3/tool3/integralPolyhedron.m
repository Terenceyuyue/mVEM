function Int = integralPolyhedron(fun,n,node3,elemf)
% Approximate integrals in a polyhedral domain
% n: n-th order quadrature rule
% fun: one or more anonymous functions, e.g. 
%    fun = @(x,y,z) [f1(x,y,z), f2(x,y,z)]
%
% Copyright (C)  Terence Yu.

%% Triangulation and geo of the polyhedon
% triangulation of faces
[Tri,index3,TriLocal] = faceTriangulation(elemf);
V = node3(index3,:);
if size(Tri,1)>4
    % centroid
    centroid3 = polycentroid3(node3,Tri);
    % triangulation of the polyhedron    
    nodeTet = [V; centroid3];
    elemTet = [TriLocal, (max(TriLocal(:))+1)*ones(size(TriLocal,1),1)];
else % tetrahedron
    nodeTet = V; 
    elemTet = 1:4;
end
% volumes of all tetrahedrons
volume = simplexvolume(nodeTet,elemTet);

%% Guass-Quadrature
[lambda,weight] = quadpts3(n); 
weight = weight(:)'; % must be a row vector
nQuad = length(weight);
Int = 0;
for p = 1:nQuad
    pz = lambda(p,1)*nodeTet(elemTet(:,1),:) ...
        +lambda(p,2)*nodeTet(elemTet(:,2),:) ...
        +lambda(p,3)*nodeTet(elemTet(:,3),:) ...
        +lambda(p,4)*nodeTet(elemTet(:,4),:);
    fp = fun(pz(:,1),pz(:,2),pz(:,3));
    Int = Int + weight(p)*fp;
end
Int = sum(repmat(volume,1,size(Int,2)).*Int, 1);