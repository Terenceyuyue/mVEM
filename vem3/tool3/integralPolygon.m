function Int = integralPolygon(fun,n,verts)
% Approximate integrals in a polygonal domain
% n: n-th order quadrature rule
% fun: one or more anonymous functions, e.g. 
%    fun = @(x,y) [f1(x,y), f2(x,y)]
%
% Copyright (C)  Terence Yu.

% centroid
nv = size(verts,1);
verts1 = verts([2:nv,1],:);
area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
ar = 0.5*abs(sum(area_components));
centroid = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
% nodeT, elemT
if nv>3
    nodeT = [verts; centroid];
    elemT = [(nv+1)*ones(nv,1),(1:nv)',[2:nv,1]'];
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