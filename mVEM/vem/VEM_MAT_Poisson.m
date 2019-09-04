function [D,Bs,G,Gs,H] = VEM_MAT_Poisson(aux)

% aux = auxgeometry(node,elem);

node = aux.node; elem = aux.elem;
elemCentroid = aux.elemCentroid;
diameter = aux.diameter;
NT = size(elem,1); 

D = cell(NT,1); 
% B = cell(NT,1);  % it is not used in the computation
Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1); H = cell(NT,1);
for iel = 1:NT
    index = elem{iel};     Nv = length(index);    
    xK = elemCentroid(iel,1); yK = elemCentroid(iel,2); hK = diameter(iel);
    
    m1 = @(x,y) 1 * ones(size(x));
    m2 = @(x,y) (x-xK)./hK;
    m3 = @(x,y) (y-yK)./hK;
    
    m = @(x,y) [m2(x,y), m3(x,y)];  % m2,m3  
    x = node(index,1); y = node(index,2);
    
    % D
    D1 = zeros(Nv,3); D1(:,1) = 1;  
    D1(:,2:3) = m(x,y);   D{iel} = D1;
    
    % B, Bs, G, Gs
    rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
    Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]];
    normVec = 0.5*[y(rotid2) - y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
    B1 = Gradm*normVec; B1s = B1; B1s(1,:) = 1/Nv;
    % B{iel} = B1;        
    Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1; 
    
    % H
    nodeT = [node(index,:);elemCentroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    m = @(x,y) [m1(x,y).*m1(x,y), m1(x,y).*m2(x,y), m1(x,y).*m3(x,y), ...
                m2(x,y).*m1(x,y), m2(x,y).*m2(x,y), m2(x,y).*m3(x,y), ...
                m3(x,y).*m1(x,y), m3(x,y).*m2(x,y), m3(x,y).*m3(x,y)];    
    H1 = zeros(3,3);
    H1(:) = integralTri(m,2,nodeT,elemT); % n = 2   
    H{iel} = H1; 
end