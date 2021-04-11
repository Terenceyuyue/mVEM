function eta = PoissonVEM_indicator(node,elem,uh,info,pde)
% This function returns the local error indicator of solving the Poisson equation 
% using virtual element method in V1
% 
% See also PoissonVEM.m
%
% Copyright (C) Terence Yu.

%% Pis and chi
% Pis
Ph = info.Ph; D = info.D;
% elementwise numerical d.o.f.s
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); 
% coefficient of elliptic projection: Ph{iel}*chi{iel}
a = cellfun(@mtimes, Ph, chi, 'UniformOutput', false); 

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
    eta3(iel) = sum((chi{iel}-D{iel}*a{iel}).^2);
end
elemRes = eta1 + eta2 + eta3;

%% elementwise edge jumps
% scaled norm vectors he*ne
e = node(edge(:,2),:)-node(edge(:,1),:);
Ne = [-e(:,2), e(:,1)];
% edgeJump
edgeJump = zeros(NE,1);
gradm = @(hK) [0 0; 1/hK 0; 0 1/hK];    
for s = 1:NE
    % element information
    k1 = edge2elem(s,1);  k2 = edge2elem(s,2);
    if k1==k2, continue; end  % initialized as zero
    h1 = diameter(k1);   h2 = diameter(k2);
    % grad of Pi(uh)
    gradLu = a{k1}'*gradm(h1); gradRu = a{k2}'*gradm(h2);
    % jump of grad(Pi(uh))
    Jumpu = gradLu-gradRu;  
    % edgeJump
    edgeJump(s) = 0.5*(dot(Jumpu,Ne(s,:)))^2; % ce = 0.5
end
% elemJump
elemJump = cellfun(@(ie) sum(edgeJump(ie)), elem2edge);

%% Local error indicator
eta = sqrt(elemRes + elemJump);