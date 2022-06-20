function Int = integralPolygon3(fun,n,P)
% Approximate integrals in a polygonal domain in 3-D
% n: n-th order quadrature rule
% fun: one or more anonymous functions, e.g. 
%    fun = @(x,y,z) [f1(x,y,z), f2(x,y,z)]
%
% Copyright (C)  Terence Yu.

% triangulation
nv = size(P,1);
tri = [ones(nv-2,1), (2:nv-1)', (3:nv)'];
% area of triangles
d12 = P(tri(:,2),:)-P(tri(:,1),:);
d13 = P(tri(:,3),:)-P(tri(:,1),:);
normal = mycross(d12,d13,2);
area = sqrt(sum(normal.^2,2))/2;
% Gauss quadrature points and weights
[lambda,weight] = quadpts(n);
% integral
Int = zeros(nv-2,1);
for p = 1:length(weight)
    pxyz = lambda(p,1)*P(tri(:,1),:) ...
        + lambda(p,2)*P(tri(:,2),:) ...
        + lambda(p,3)*P(tri(:,3),:);
    f = fun(pxyz(:,1), pxyz(:,2), pxyz(:,3));
    Int = Int + weight(p)*f;
end
Int = sum(area.*Int);