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
nodeT = [verts; centroid];
elemT = [(nv+1)*ones(nv,1),(1:nv)',[2:nv,1]'];
% integral
[lambda,weight] = quadpts(n);
NT = size(elemT,1); Int = 0;
for iel = 1:NT
    vT = nodeT(elemT(iel,:),:);
    area = 0.5*abs(det([[1;1;1],vT]));    
    xy = lambda*vT;
    for p = 1:size(xy,1)
        f = fun(xy(p,1),xy(p,2));
        Int = Int + area*weight(p)*f;
    end
end