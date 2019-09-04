function [D,Bs,G,Gs,H0,C0] = VEM_MAT_elasticity_tensor(aux)

% aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
elemCentroid = aux.elemCentroid;
diameter = aux.diameter; area = aux.area;
NT = size(elem,1); 

D = cell(NT,1);
% B = cell(NT,1);  % not used in the computation
Bs = cell(NT,1);
G = cell(NT,1); Gs = cell(NT,1);
H0 = cell(NT,1); C0 = cell(NT,1);
for iel = 1:NT
    index = elem{iel};     Nv = length(index);
    xK = elemCentroid(iel,1); yK = elemCentroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    m2 = @(x,y) (x-xK)./hK;  m3 = @(x,y) (y-yK)./hK;
    m = @(x,y) [m2(x,y), m3(x,y)];  % m2,m3
    
    % D
    D1 = zeros(Nv,3); D1(:,1) = 1;
    D1(:,2:3) = m(x,y);   D1 = blkdiag(D1,D1);  D{iel} = D1;
    
    % H0,C0
    rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
    normVec = 0.5*[y(rotid2) - y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
    H0{iel} = area(iel); C01 = reshape(normVec',1,[]); 
    C0{iel} = C01;
    
    % B,Bs,G,Gs
    E = zeros(6,4); % E = [E11,E12,E21,E22]
    E(2,1) = 1/hK;  E([3,5],[2,3]) = 1/hK; E(6,4) = 1/hK;
    B1 = [E(:,1)*C01(1:Nv)+E(:,2)*C01(Nv+1:end), ...
        E(:,3)*C01(1:Nv)+E(:,4)*C01(Nv+1:end)];
    B0 = zeros(3,2*Nv); B0(1,1:Nv) = 1; B0(2,Nv+1:end) = 1;
    B0(3,1:Nv) = -y;    B0(3,Nv+1:end) = x;
    % B{iel} = B1;
    B1s = B1; B1s([1,3,4],:) = B0;   Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;   
end
