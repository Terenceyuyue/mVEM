function Int = integralPolygonTri(fun,n,verts)

nv = size(verts,1);

% nodeT, elemT
if nv>3 % triangulation
    edges = [1:nv; [2:nv,1]]';
    DT = delaunayTriangulation(verts, edges);
    nodeT = DT.Points;
    elemT = DT.ConnectivityList;
else
    nodeT = verts;  elemT = 1:3;
end
% area of all triangles
area = simplexvolume(nodeT,elemT); 
% integral
[lambda,weight] = quadpts(n);
Int = 0; NT = size(elemT,1);
for p = 1:length(weight)
    pxy = lambda(p,1)*nodeT(elemT(:,1),:) ...
        + lambda(p,2)*nodeT(elemT(:,2),:) ...
        + lambda(p,3)*nodeT(elemT(:,3),:);
    for iel = 1:NT
        fp = fun(pxy(iel,1),pxy(iel,2));
        Int = Int + weight(p)*area(iel)*fp;
    end
end
Int