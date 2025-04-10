function [eta,etaStab] = IPVEM_indicator(node,elem,uh,info,pde,bdStruct)
% This function returns the local error indicator of solving the biharmonic equation 
% using interior penalty virtual element method in the lowest-order case: k = 2.
% 
% See also IPVEM.m
%
% Copyright (C) Terence Yu.

%% Pis and chi
% Pis
Ph1 = info.Ph(:,1);  % H1
Dh = info.Dh;   Hh = info.Hh;
% elementwise numerical d.o.f.s
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false);
% coefficient of elliptic projection: Ph{iel}*chi{iel}
a = cellfun(@mtimes, Ph1, chi, 'UniformOutput', false); 
% penalty parameter
cp = info.cp;

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
Nm = 6;

%% elementwise residuals
[eta4,eta5,eta6] = deal(zeros(NT,1));  % squared
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    % scaled monomials
    m1 = @(x,y) 1+0*x;                     
    m2 = @(x,y) (x-xK)/hK;                 
    m3 = @(x,y) (y-yK)/hK;                 
    m4 = @(x,y) (x-xK).^2/hK^2;            
    m5 = @(x,y) (x-xK).*(y-yK)/hK^2;      
    m6 = @(x,y) (y-yK).^2/hK^2;            
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)];
    
    % eta4
    dofK = chi{iel};   
    dofPi1 = Dh{iel}*a{iel}; % Pi_1*chi = D*Pi_{s,1}*chi = D*a1
    eta4(iel) = hK^(-2)*sum((dofK-dofPi1).^2);

    % eta5
    H = Hh{iel};
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    f = integralTri(fm,5,nodeT,elemT);
    c = H\f(:);
    Pif = @(x,y) m(x,y)*c;
    eta5f = @(x,y) (pde.f([x,y])-Pif(x,y)).^2;
    eta5(iel) = hK^4*integralTri(eta5f,5,nodeT,elemT);

    % eta6 (\Delta^2 \Pi_h uh = 0)
    eta6f = @(x,y) Pif(x,y).^2;
    eta6(iel) = hK^4*integralTri(eta6f,5,nodeT,elemT);
end
elemRes = eta4 + eta5 + eta6;

%% elementwise edge jumps
% scaled norm vectors he*ne
bdEdgeIdx = bdStruct.bdEdgeIdx;
bdEdge = bdStruct.bdEdge;
edge(bdEdgeIdx,:) = bdEdge; % adjustment on the boundary
v12 = node(edge(:,1),:)-node(edge(:,2),:);
Ne = [-v12(:,2),v12(:,1)];
he = vecnorm(v12,2,2);
ne = Ne./he;
% Gauss-Lobatto: Simpson for k = 2
ng = 3;
w = [1/6, 4/6, 1/6];

[eta1,eta2,eta3] = deal(zeros(NE,1));  % squared
for ie = 1:NE
    % --------- info of left and right elements -----------
    k1 = edge2elem(ie,1);  k2 = edge2elem(ie,2);
    xK1 = centroid(k1,1);  yK1 = centroid(k1,2); 
    xK2 = centroid(k2,1);  yK2 = centroid(k2,2); 
    hK1 = diameter(k1);    hK2 = diameter(k2);
    nx = ne(ie,1);  ny = ne(ie,2);
    z1 = node(edge(ie,1),:);  z2 = node(edge(ie,2),:);  ze = (z1+z2)/2;
    xg = [z1(1),ze(1),z2(1)]; yg = [z1(2),ze(2),z2(2)]; % Gauss-Lobatto

    % --------- scaled monomials --------
    mx = @(x,y,xK,yK,hK) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x];
    my = @(x,y,xK,yK,hK) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2]; 
    mxx = @(x,y,xK,yK,hK) [0*x, 0*x, 0*x, 2/hK^2+0*x, 0*x, 0*x];
    mxy = @(x,y,xK,yK,hK) [0*x, 0*x, 0*x, 0*x, 1/hK^2+0*x, 0*x];
    myy = @(x,y,xK,yK,hK) [0*x, 0*x, 0*x, 0*x, 0*x, 2/hK^2+0*x];

    % eta1    
    ak1 = a{k1};  ak2 = a{k2};
    for s = 1:ng
        x = xg(s); y = yg(s);
        J1 = mx(x,y,xK1,yK1,hK1)*ak1*nx + my(x,y,xK1,yK1,hK1)*ak1*ny;
        J2 = mx(x,y,xK2,yK2,hK2)*ak2*nx + my(x,y,xK2,yK2,hK2)*ak2*ny;
        Jn = J1-J2;
        eta1(ie) = eta1(ie) + cp(ie)*w(s)*Jn^2;
    end

    % eta2
    for s = 1:ng
        x = xg(s); y = yg(s);
        Jxx1 = mxx(x,y,xK1,yK1,hK1)*ak1*nx*nx; 
        Jxx2 = mxx(x,y,xK2,yK2,hK2)*ak2*nx*nx;
        Jxy1 = mxy(x,y,xK1,yK1,hK1)*ak1*nx*ny;
        Jxy2 = mxy(x,y,xK2,yK2,hK2)*ak2*nx*ny;
        Jyy1 = myy(x,y,xK1,yK1,hK1)*ak1*ny*ny;
        Jyy2 = myy(x,y,xK2,yK2,hK2)*ak2*ny*ny;
        Jnn = (Jxx1+2*Jxy1+Jyy1) - (Jxx2+2*Jxy2+Jyy2);
        eta2(ie) = eta2(ie) + he(ie)^2*w(s)*Jnn^2;
    end

    % eta3 = 0    
end
edgeJump = eta1 + eta2 + eta3;
% elemJump
elemJump = cellfun(@(ie) sum(edgeJump(ie)), elem2edge);

%% Local error indicator
eta = sqrt(elemRes + elemJump);
etaStab = sqrt(elemRes - eta4 + elemJump);