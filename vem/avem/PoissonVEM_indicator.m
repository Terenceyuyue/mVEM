function eta = PoissonVEM_indicator(node,elem,uh,info,pde)
% This function returns the local error indicator of Poisson equation in 2-D.
% 
% Copyright (C) Terence Yu.

%% Pis and chi
Ph = info.Ph; D = info.D;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s

%% auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
centroid = aux.centroid;  diameter = aux.diameter;
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge; 
elem2edge = auxT.elem2edge; edge2elem = auxT.edge2elem;
% number
NT = size(elem,1); NE = size(edge,1);

%% elementwise residuals
eta1 = zeros(NT,1);  eta2 = zeros(NT,1);  eta3 = zeros(NT,1);  % squared
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    % scaled monomials
    m1 = @(x,y)  1+0*x;
    m2 = @(x,y) (x-xK)./hK;
    m3 = @(x,y) (y-yK)./hK;
    mm = @(x,y) [m1(x,y).*m1(x,y), m1(x,y).*m2(x,y), m1(x,y).*m3(x,y), ...
                 m2(x,y).*m1(x,y), m2(x,y).*m2(x,y), m2(x,y).*m3(x,y), ...
                 m3(x,y).*m1(x,y), m3(x,y).*m2(x,y), m3(x,y).*m3(x,y)];
    % coefficient of elliptic projection
    a = Ph{iel}*chi{iel}; % Pis = Ph{iel}
    % eta1
    H = zeros(3,3);
    H(:) = integralTri(mm,3,nodeT,elemT); % n = 3
    fm = @(x,y) repmat(pde.f([x,y]),1,3).*[m1(x,y),m2(x,y),m3(x,y)];
    f = integralTri(fm,3,nodeT,elemT); f = f(:);
    c = H\f;
    Pif = @(x,y) c(1)*m1(x,y)+c(2)*m2(x,y)+c(3)*m3(x,y);
    eta1f = @(x,y) (pde.f([x,y])-Pif(x,y)).^2;
    eta1(iel) = hK^2*integralTri(eta1f,3,nodeT,elemT);
    % eta2
    eta2f = @(x,y) Pif(x,y).^2;  % Laplace(Pi(uh)) = 0
    eta2(iel) = hK^2*integralTri(eta2f,3,nodeT,elemT);
    % eta3
    eta3(iel) = sum((chi{iel}-D{iel}*a).^2);
end
elemRes = eta1 + eta2 + eta3;

%% elementwise edge jumps
% scaled norm vectors he*ne
e = node(edge(:,2),:)-node(edge(:,1),:);
Ne = [-e(:,2), e(:,1)];
% edgeJump
edgeJump = zeros(NE,1);
for s = 1:NE
    % element information
    k1 = edge2elem(s,1);  k2 = edge2elem(s,2); 
    hK1 = diameter(k1);   hK2 = diameter(k2);
    % coefficient of elliptic projection
    gradmk1 = [0 0; 1/hK1 0; 0 1/hK1];
    gradmk2 = [0 0; 1/hK2 0; 0 1/hK2];
    % grad of Pi(uh)
    ak1 = Ph{k1}*chi{k1};  ak2 = Ph{k2}*chi{k2};
    gradLu = ak1'*gradmk1; gradRu = ak2'*gradmk2;
    % jump of grad(pi(uh))
    Jumpu = gradLu-gradRu; ce = 0.5;
    if k1==k2, ce = 0; end    
    % edgeJump
    edgeJump(s) = ce*(dot(Jumpu,Ne(s,:)))^2;
end
% elemJump
elemJump = cellfun(@(ie) sum(edgeJump(ie)), elem2edge);

%% Local error indicator
eta = sqrt(elemRes + elemJump);