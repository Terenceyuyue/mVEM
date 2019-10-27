function [node,elem] = PolyMesher(Domain,NT,MaxIter)

% ------- Generate intial pointset --------
P = PolyMesher_init_Pointset(Domain,NT);

% ------- Lloyd's iteration --------
Iter = 0; Err = 1; Tol = 5e-6; 
Pc = P;
while(Iter<=MaxIter && Err>Tol)
    % Lloyd's update
    P = Pc;
    % Compute the reflection pointset
    R_P = PolyMesher_Reflect(P,Domain);
    % Compute the centroid of Voronoi cell
    [Pc,Err,node,elem] = PolyMesher_VoroCentroid(P,R_P);
    Iter = Iter+1;
end

% ------------ node, elem -------------------
[node,elem] = PolyMesher_rm_smalledge(NT,node,elem);

