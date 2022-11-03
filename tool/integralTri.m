function Int = integralTri(fun,n,nodeT,elemT)
%integralTri approximates integrals in a polygonal domain with trianguation (nodeT,elemT).
% n: n-th order quadrature rule
% fun: one or more anonymous functions, e.g. fun = @(x,y) [f1(x,y), f2(x,y)]

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