function [node,elem] = PolyMesher(Domain,MaxIter,varargin)

% ----------------- Generate intial pointset ----------------
nvar = length(varargin);
switch nvar
    case 1
        NT = varargin{1};
        P = PolyMesher_init_Pointset(Domain,NT);
    case 2
        Nx = varargin{1}; Ny = varargin{2};
        P = PolyMesher_init_Pointset(Domain,Nx,Ny);
        NT = size(P,1);
end

% -------------------- Lloyd's iteration ---------------
Iter = 0; Err = 1; Tol = 5e-6;
BdBox = Domain('BdBox'); 
Area = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
while(Iter<=MaxIter && Err>Tol)
    % Compute the reflection pointset
    R_P = PolyMesher_Reflect(P,Domain,Area);
    % Construct the Voronoi diagram
    [node,elem] = voronoin([P;R_P],{'Qbb','Qz'});
    % Compute the centroid of Voronoi cell (Lloyd's update)
    [P,Area,Err] = PolyMesher_VoroCentroid(P,node,elem);
    Iter = Iter+1;
end

% ----------- Remove small edges to obtain (node,elem) -------------
[node,elem] = PolyMesher_rm_smalledge(NT,node,elem);

