function [Pc,Area,Err] = PolyMesher_VoroCentroid(P,node,elem)
% Compute the centroid of Voronoi cell
NT = size(P,1);
Pc = zeros(NT,2); A = zeros(NT,1);
for iel = 1:NT    
    vx = node(elem{iel},1);  vy = node(elem{iel},2); nv = length(elem{iel});
    vxS = vx([2:nv,1]); vyS = vy([2:nv,1]);
    temp = vx.*vyS - vy.*vxS;
    A(iel) = 0.5*sum(temp);
    Pc(iel,:) = 1/(6*A(iel))*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
end
Area = sum(abs(A));
Err = sqrt(sum((A.^2).*sum((Pc-P).^2,2)))*NT/Area^1.5;