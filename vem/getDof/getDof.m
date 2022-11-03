function chie = getDof(node,elem,pde,kOrder)

if nargin==3, kOrder=1; end

uexact = pde.uexact;

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge;
% numbers
NT = size(elem,1);

%% k = 1
if kOrder==1
    chie = uexact(node);
end

%% k = 2
if kOrder==2
    % chi1: evaluation at all vertices
    chi1 = uexact(node);
    % chi2: evaluation at all mid-edge points
    z1 = node(edge(:,1),:); z2 = node(edge(:,2),:); zc = (z1+z2)/2;
    chi2 = uexact(zc);
    % chi3: moments on element K    
    chi3 = zeros(NT,1);
    fun = @(x,y) uexact([x,y]);
    for iel = 1:NT
        index = elem{iel};     Nv = length(index);
        nodeT = [node(index,:); centroid(iel,:)];
        elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
        chi3(iel) = 1/area(iel)*integralTri(fun,3,nodeT,elemT);
    end
    chie = [chi1; chi2; chi3];
end

%% k = 3
if kOrder==3
    % interior Gauss-Lobatto points
    r = [-1, -1/sqrt(5), 1/sqrt(5), 1]; % [-1,1]
    r = (r(2:3)+1)/2; % [0,1]: gives the ratios of interior nodes
    z1 = node(edge(:,1),:); z2 = node(edge(:,2),:); 
    za = z1 + r(1)*(z2-z1); 
    zb = z1 + r(2)*(z2-z1);
    
    % evaluation at all vertices and interior Gauss-Lobatto points
    chi1 = uexact(node);    
    chi2 = uexact(za);
    chi3 = uexact(zb);
    
    % moments on K
    chiRes = zeros(NT,3);
    for iel = 1:NT
        % element information
        index = elem{iel};  Nv = length(index);
        xK = centroid(iel,1); yK = centroid(iel,2);  hK = diameter(iel);
        nodeT = [node(index,:);centroid(iel,:)]; % triangulation of K
        elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
        
        % scaled monomials
        m1 = @(x,y) 1+0*x;
        m2 = @(x,y) (x-xK)/hK;
        m3 = @(x,y) (y-yK)/hK;
        m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];        
        
        % moments on K
        fun = @(x,y) repmat(uexact([x,y]),1,3).*m(x,y);
        chiRes(iel,:) = 1/area(iel)*integralTri(fun,4,nodeT,elemT);
    end
    
    % chie
    chie = [chi1; chi2; chi3; chiRes(:)];
end


