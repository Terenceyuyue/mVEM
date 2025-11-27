function eta = IPVEM_Zhao_indicator(node,elem,uh,info,pde,bdStruct)
% This function returns the local error indicator of solving the biharmonic equation 
% using interior penalty virtual element method in the lowest-order case: k = 2.
% 
% See also IPVEM_Zhao.m
%
% Copyright (C) Terence Yu.

%% Pis and chi
% Pis
Ph1 = info.Ph(:,1);  % H1
Ph2 = info.Ph(:,2);  % H2
Dh = info.Dh;   Hh = info.Hh;
% elementwise numerical d.o.f.s
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false);
% coefficient of elliptic projection: Ph{iel}*chi{iel}
a1 = cellfun(@mtimes, Ph1, chi, 'UniformOutput', false); 
a2 = cellfun(@mtimes, Ph2, chi, 'UniformOutput', false); 
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
    dofPi2 = Dh{iel}*a2{iel}; % Pi_2*chi = D*Pi_{s,2}*chi = D*a2
    eta4(iel) = hK^(-2)*sum((dofK-dofPi2).^2);

    % eta5
    H = Hh{iel};
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    f = integralTri(fm,5,nodeT,elemT);
    c = H\f(:);
    Pif = @(x,y) m(x,y)*c;
    eta5f = @(x,y) (pde.f([x,y])-Pif(x,y)).^2;
    eta5(iel) = hK^4*integralTri(eta5f,5,nodeT,elemT);

    % eta6
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
% Gauss-Lobatto for integral: not for the space
ng = 4;
r = [-1, -1/sqrt(5), 1/sqrt(5), 1]; % [-1,1]
r = (r+1)/2; % [0,1]: gives the ratios
w = [1/12, 5/12, 5/12, 1/12];
cr = -1/2+r;

[eta1,eta2,eta3] = deal(zeros(NE,1));  % squared
for ie = 1:NE
    % --------- info of left and right elements -----------
    k1 = edge2elem(ie,1);  k2 = edge2elem(ie,2);
    xK1 = centroid(k1,1);  yK1 = centroid(k1,2); 
    xK2 = centroid(k2,1);  yK2 = centroid(k2,2); 
    hK1 = diameter(k1);    hK2 = diameter(k2);
    nx = ne(ie,1);  ny = ne(ie,2);
    z1 = node(edge(ie,1),:);  z2 = node(edge(ie,2),:);  ze = (z1+z2)/2;
    % Gauss-Lobatto for space
    x = [z1(1),ze(1),z2(1)]; 
    y = [z1(2),ze(2),z2(2)]; 
    % Gauss-Lobatto for integral
    za = z1+r(2)*(z2-z1); 
    zb = z1+r(3)*(z2-z1);
    xg = [z1(1),za(1),zb(1),z2(1)]; 
    yg = [z1(2),za(2),zb(2),z2(2)]; 

    % --------- scaled monomials --------
    mx = @(x,y,xK,yK,hK) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x];
    my = @(x,y,xK,yK,hK) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2]; 
    mxx = @(x,y,xK,yK,hK) [0*x, 0*x, 0*x, 2/hK^2+0*x, 0*x, 0*x];
    mxy = @(x,y,xK,yK,hK) [0*x, 0*x, 0*x, 0*x, 1/hK^2+0*x, 0*x];
    myy = @(x,y,xK,yK,hK) [0*x, 0*x, 0*x, 0*x, 0*x, 2/hK^2+0*x];

    % eta0: neglected when consider the H2-semi-norm    

    % eta1  
    % Pi_{1,s}*chi
    a1k1 = a1{k1};  a1k2 = a1{k2};
    % \partial_{n_e} (Pi_1 uh) = (\partial_{n_e} m')*Pi_{1,s}*chi
    gn1 = @(x,y,xK,yK,hK) mx(x,y,xK,yK,hK)*a1k1*nx + my(x,y,xK,yK,hK)*a1k1*ny; 
    gn2 = @(x,y,xK,yK,hK) mx(x,y,xK,yK,hK)*a1k2*nx + my(x,y,xK,yK,hK)*a1k2*ny;
    % expansion coefficients for space: k=2
    mm = zeros(4,4); % (k+1)-th order polynomials: a0,...,a_{k+1}
    for i = 0:3
        for j = i:3
            cij = i+j+1;
            mm(i+1,j+1) = he(ie)/cij*(1/2^cij - 1/(-2)^cij);
            mm(j+1,i+1) = mm(i+1,j+1);
        end
    end
    Aa = zeros(4,4);    
    Aa(1,:) = cr(1).^(0:3);  % zi
    Aa(2,:) = cr(end).^(0:3); % z_{i+1}
    Aa(3:4,:) = mm(1:2,:);
    [b1,b2] = deal(zeros(4,1));     
    b1(1) = gn1(x(1),y(1),xK1,yK1,hK1);  % zi
    b1(2) = gn1(x(3),y(3),xK1,yK1,hK1); % z_{i+1}
    b2(1) = gn2(x(1),y(1),xK2,yK2,hK2);  % zi
    b2(2) = gn2(x(3),y(3),xK2,yK2,hK2); % z_{i+1}
    for i = 1:3
        gni1 = gn1(x(i),y(i),xK1,yK1,hK1);
        gni2 = gn2(x(i),y(i),xK2,yK2,hK2);
        b1(3) = b1(3) + he(ie)*w(i)*gni1;
        b1(4) = b1(4) + he(ie)*w(i)*cr(i)*gni1;
        b2(3) = b2(3) + he(ie)*w(i)*gni2;
        b2(4) = b2(4) + he(ie)*w(i)*cr(i)*gni2;
    end
    c1 = Aa\b1;  c2 = Aa\b2;
    % compute eta1    
    for s = 1:ng
        Jn = cr(s).^(0:3)*(c1-c2);
        eta1(ie) = eta1(ie) + cp(ie)*w(s)*Jn^2;
    end

    % eta2
    a2k1 = a2{k1};  a2k2 = a2{k2};
    for s = 1:ng
        x = xg(s); y = yg(s);
        Jxx1 = mxx(x,y,xK1,yK1,hK1)*a2k1*nx*nx; 
        Jxx2 = mxx(x,y,xK2,yK2,hK2)*a2k2*nx*nx;
        Jxy1 = mxy(x,y,xK1,yK1,hK1)*a2k1*nx*ny;
        Jxy2 = mxy(x,y,xK2,yK2,hK2)*a2k2*nx*ny;
        Jyy1 = myy(x,y,xK1,yK1,hK1)*a2k1*ny*ny;
        Jyy2 = myy(x,y,xK2,yK2,hK2)*a2k2*ny*ny;
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