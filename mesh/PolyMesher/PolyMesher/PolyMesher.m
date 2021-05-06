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
end

PFix = Domain('PFix');   P = [PFix; P];  
nFix = size(PFix,1);     NT = size(P,1);

% -------------------- Lloyd's iteration ---------------
Iter = 0; Err = 1; Tol = 1e-4;
BdBox = Domain('BdBox'); 
Area = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
while(Iter<=MaxIter && Err>Tol)
    % Compute the reflection pointset
    R_P = PolyMesher_Reflect(P,Domain,Area);
    % Construct the Voronoi diagram
    [node,elem] = voronoin([P;R_P],{'Qbb','Qz'});
    % Compute the centroid of Voronoi cell (Lloyd's update)
    [P,Area,Err] = PolyMesher_VoroCentroid(P,node,elem);
    % Include the fixed seeds
    P(1:nFix,:) = PFix;
    Iter = Iter+1;
    if mod(Iter,10)==0
        fprintf('Iter: %3d   Error: %1.3e\n',Iter,Err);
    end
end
if NT<=1000
    clf; showmesh(node,elem(1:NT)); 
end

% ----------- Remove small edges to obtain (node,elem) -------------
[node,elem] = PolyMesher_rm_smalledge(node,elem(1:NT));
