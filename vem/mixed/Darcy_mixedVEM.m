function [uh,ph,info] = Darcy_mixedVEM(node,elem,pde,bdStruct)
%Darcy_mixedVEM solves Darcy problem using the mixed virtual element method
% in the lowest order
%
%     -div(K*grad(p)) = f   in Omega
%     Neumann boundary condition   K*grad(u)*n=g on \Gamma
%
%  The mixed formulation is: 
%      u = K*grad(p), div u = -f  in \Omega
%      u*n = g                    on \Gamma
%  The velocity u = K*grad(p) is approximated by the lowest order
%  H(div)-conforming virtual element (k = 1)
%      V_k(K) = { v\in H(div;K) \cap H(rot;K): 
%                   v\dot n|_e\in P_k(e), e \subsect \partial K,
%                   div(v)|_K\in P_{k-1}(K), rot(v)|_K\in P_{k-1}(K) }
%  and p by piecewise constant element (k = 1)
%      Q_k(K) = P_{k-1}(K).
%
%   References
%   F. Brezzi, R.S. Falk and L.D. Marini, "Basic principles of mixed
%   virtual element methods", ESAIM Math. Model. Numer. Anal.,
%   Vol 48. No 4., pp. 1227--1240, 2014.
%
% Copyright (C)  Terence Yu. 

K = pde.K;  % coefficient matrix

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
NT = size(elem,1); NE = size(edge,1);
NNdofA = 2*NE+NT;  NNdofB = NT;  NNdof = NNdofA + NNdofB;
Nm = 6;   Nmh = Nm-1;

%% Get elementwise signs of basis functions
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
sgnBase = cell(NT,1);
for iel = 1:NT
    index = elem{iel};   Nv = length(index);   NdofA = 2*Nv+1;
    sgnedge = sign(diff(index([1:Nv,1])));
    id = elem2edge{iel}; sgnbd = E(id); sgnedge(sgnbd) = 1;
    sgnelem = ones(NdofA,1); sgnelem(1:Nv) = sgnedge; 
    sgnBase{iel} = sgnelem;
end

%% Compute and assemble the linear system
elemLen = cellfun('length',elem); 
nnzA = sum((2*elemLen+1).^2);  nnzB = sum((2*elemLen+1)*1); 
iiA = zeros(nnzA,1); jjA = zeros(nnzA,1); ssA = zeros(nnzA,1);
iiB = zeros(nnzB,1); jjB = zeros(nnzB,1); ssB = zeros(nnzB,1);
elemb = zeros(NNdofB,1);  Fb = zeros(NNdofB,1); 
idA = 0; idB = 0;  ib = 0;
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
for iel = 1:NT
    % ------- element information --------
    index = elem{iel};     Nv = length(index);   NdofA = 2*Nv+1;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); 
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    Te = [-Ne(:,2), Ne(:,1)]; % he*te
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % ------- scaled monomials --------     
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;
    m4 = @(x,y) (x-xK).^2./hK^2; 
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2; 
    m6 = @(x,y) (y-yK).^2./hK^2;
    % \hat{m}_a = K*grad(hK*m_{a+1})
    mh1 = @(x,y) [K(1,1)+0*x; K(2,1)+0*x];
    mh2 = @(x,y) [K(1,2)+0*x; K(2,2)+0*x];
    mh3 = @(x,y) [2*K(1,1)*m2(x,y); 2*K(2,1)*m2(x,y)];
    mh4 = @(x,y) [K(1,1)*m3(x,y)+K(1,2)*m2(x,y); K(2,1)*m3(x,y)+K(2,2)*m2(x,y)];
    mh5 = @(x,y) [2*K(1,2)*m3(x,y); 2*K(2,2)*m3(x,y)];
    mh = @(x,y) [mh1(x,y), mh2(x,y), mh3(x,y), mh4(x,y), mh5(x,y)];  
    
    % -------- transition matrix ----------
    D = zeros(NdofA,Nmh); 
    for i = 1:Nv % loop of edges
        % v at z_i, z_{i+1} for v = [mK1,...,mK5]
        va = mh(x(v1(i)),y(v1(i))); vb = mh(x(v2(i)),y(v2(i)));
        % chi_i, i = 1,...,Nv
        D(i,:) =  1/2*Ne(i,:)*(va+vb);
        % chi_{Nv+i}, i = 1,...,Nv
        D(Nv+i,:) =  1/12*Ne(i,:)*(vb-va);
        % chi_{2Nv+1}
        D(end,:) = D(end,:) + 1/2*Te(i,:)*(va+vb);
    end
    
    % --------- elliptic projection -----------
    % first term
    m = @(x,y) [m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)]; % m_{a+1},...
    ci = zeros(1,NdofA); ci(1:Nv) = 1/area(iel);
    Intm = integralTri(m,3,nodeT,elemT);
    I1 = Intm'*ci;
    % second term
    rij = zeros(NdofA,Nv); rij(1:Nv,:) = eye(Nv);
    sij = zeros(NdofA,Nv); sij(Nv+1:2*Nv,:) = eye(Nv);
    rija = rij - 6*sij;
    rije = rij;
    rijb = rij + 6*sij;
    I2 = zeros(Nmh,NdofA);
    for j = 1:Nv  
        ma = m(x(v1(j)),y(v1(j)));
        me = m(xe(j),ye(j));
        mb = m(x(v2(j)),y(v2(j)));
        rja = rija(:,j)';  rje = rije(:,j)';  rjb = rijb(:,j)';
        I2 = I2 + 1/6*(ma'*rja + 4*me'*rje + mb'*rjb);
    end
    % Bs=B, Gs=G: no constraint
    Bs = hK*(-I1+I2); Gs = Bs*D;    
    
    % -------- sign matrix and sign vector -------
    sgnelem = sgnBase{iel}; 
    sgnK = sgnelem*sgnelem'; 
    
    % ------------- stiffness matrix -------------
    % projection
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));    
    % stiffness matrix A
    AK  = Pis'*Gs*Pis + norm(inv(K),'fro')*(I-Pi)'*(I-Pi);  % G = Gs
    AK = AK.*sgnK;
    AK = reshape(AK',1,[]); % straighten as row vector for easy assembly 
    % stiffness matrix B
    BK = zeros(NdofA,1); BK(1:Nv) = 1;
    BK = BK.*sgnelem;
    BK = reshape(BK',1,[]); % straighten as row vector for easy assembly 
    
    % -------- load vector f ----------
    fxy = @(x,y) pde.f([x,y]); % f(p) = f([x,y])
    rhs = integralTri(fxy,3,nodeT,elemT); rhs = rhs';
    fK = -rhs;
    
    % ------ assembly index for bilinear forms --------
    NdofA = 2*Nv+1; NdofB = 1;
    indexDofA = [elem2edge{iel},elem2edge{iel}+NE,iel+2*NE];
    indexDofB = iel;
    iiA(idA+1:idA+NdofA^2) = reshape(repmat(indexDofA, NdofA,1), [], 1);
    jjA(idA+1:idA+NdofA^2) = repmat(indexDofA(:), NdofA, 1);
    ssA(idA+1:idA+NdofA^2) = AK(:);
    idA = idA + NdofA^2;
    iiB(idB+1:idB+NdofA*NdofB) = reshape(repmat(indexDofA, NdofB,1), [], 1);
    jjB(idB+1:idB+NdofA*NdofB) = repmat(indexDofB(:), NdofA, 1);
    ssB(idB+1:idB+NdofA*NdofB) = BK(:);
    idB = idB + NdofA*NdofB;
    
    % ------- assembly index for rhs -------
    elemb(ib+1:ib+NdofB) = indexDofB(:);
    Fb(ib+1:ib+NdofB) = fK(:);
    ib = ib + NdofB;
    
    % ------- matrix for error evaluation -------
    sgnPis = repmat(sgnelem',size(Pis,1),1);
    Ph{iel} = sgnPis.*Pis; 
    elem2dof{iel} = indexDofA;
end
A = sparse(iiA,jjA,ssA,NNdofA,NNdofA);
B = sparse(iiB,jjB,ssB,NNdofA,NNdofB);
FB = accumarray(elemb,Fb,[NNdofB 1]);
d = area;  % for Lagrange multiplier 

%% Get block linear system
kk = sparse(NNdof+1,NNdof+1);  ff = zeros(NNdof+1,1);
kk(1:NNdofA,1:NNdofA) = A;
kk(1:NNdofA, (1:NNdofB)+NNdofA) = B;
kk((1:NNdofB)+NNdofA, 1:NNdofA) = B';
kk((1:NNdofB)+NNdofA, end) = d;
kk(end, (1:NNdofB)+NNdofA) = d';
ff((1:NNdofB)+NNdofA) = FB;

%% Apply Dirichlet boundary conditions 
% bdDof, freeDof
bdEdge = bdStruct.bdEdge;  bdEdgeIdx = bdStruct.bdEdgeIdx;
id = [bdEdgeIdx; bdEdgeIdx+NE];
isBdDof = false(NNdof+1,1); isBdDof(id) = true;
bdDof = (isBdDof); freeDof = (~isBdDof);
% bdval
u = pde.uexact;
z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:); ze = (z1+z2)./2;
e = z1-z2;  % e = z2-z1
Ne = [-e(:,2),e(:,1)];  
chi1 = 1/6*sum((u(z1)+4*u(ze)+u(z2)).*Ne,2); % u*n = g
chi2 = 1/12*sum((u(z2)-u(z1)).*Ne,2);
bdval = [chi1;chi2]; 
% sol
sol = zeros(NNdof+1,1); 
sol(bdDof) = bdval;
ff = ff - kk*sol;

%% Set solver
sol(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
uh = sol(1:NNdofA); % u = [u1,u2]
ph = sol(NNdofA+1:end-1);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; %info.freeDof = freeDof;