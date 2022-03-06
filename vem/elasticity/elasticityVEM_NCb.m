function [u,info] = elasticityVEM_NCb(node,elem,pde,bdStruct)
%elasticityVEM_NC solves linear elasticity equation of tensor form using
% virtual element method in the local nonconforming space of form (k = 1):
%
%    V_k(K) = { v \in H^1(K):  \Delta v = P_{k-2}(K),
%               grad(v)*n|_e \in P_{k-1}(e) for interior edges,
%               v|_e \in P_k(e) for boundary edges  }
% Note that the boundary edges are treated continuously.
%
% The problem is
%
%       u = [u1, u2]
%       -div sigma(u) = f in \Omega,
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D,
%       Neumann boundary condition  sigma(u)n = [g1_N, g2_N] on \Gamma_N
%
% Copyright (C)  Terence Yu, Li Jia.

%% PDE data
f1xy = @(x,y) pde.f([x,y])*[1;0];
f2xy = @(x,y) pde.f([x,y])*[0;1];
para = pde.para;

%% Get auxiliary data
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge;  edge = auxT.edge;
% auxgeometry
aux = auxgeometry(node,elem);
centroid = aux.centroid;  diameter = aux.diameter;  area = aux.area;
% numbers
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Get the connection number of d.o.f.s
% Boundary condition
bdEdgeIdx = bdStruct.bdEdgeIdx;
bdEdge = bdStruct.bdEdge;
% Global logical vector of edges
isDofEdge = true(1,NE); % true for moment of edges
isDofEdge(bdEdgeIdx) = false;
% Global logical vector of vertices
isDofVertice = false(1,N);
isDofVertice(bdEdge(:)) = true; % true for all boundary vertices
% Global logical vector of edges and vertices
isDof = [isDofEdge, isDofVertice];
NNdof = sum(isDof);  % number of d.o.f.s of scalar case
NNdof2 = 2*NNdof; % number of d.o.f.s of vectorial case
% Connection number of d.o.f.s in form of edges and vertices
idDof = zeros(1,NE+N); % set the redundant as zero
idDof(isDof) = 1:NNdof;

%% Compute projection matrices
%% Also derive load vector and local-global index
D = cell(NT,1);  Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1);
H0 = cell(NT,1); C0 = cell(NT,1);
elem2dof = cell(NT,1); % local-global index
belem = cell(NT,1);  % load vector
NdofElem = zeros(NT,1); % numbers of local d.o.f.s on each cell

for iel = 1:NT
    
    % ------------- element information -------------
    % commonly used information
    index = elem{iel};  indexEdge = elem2edge{iel};
    Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1];  % starting and ending
    Ne = [y(v2)-y(v1), x(v1)-x(v2)];
    he = sqrt(Ne(:,1).^2 + Ne(:,2).^2);
    Te = [-Ne(:,2), Ne(:,1)];
    % for integration over auxiliary triangles
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    % scaled monomials
    m1 = @(x,y) 1 + 0*x; m2 = @(x,y) (x-xK)./hK; m3 = @(x,y) (y-yK)./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    
    % ------------ local logical vectors of d.o.f.s ---------------
    % isEdge, isVertice, isLocalDof
    isEdge = isDofEdge(indexEdge);
    isVertice = false(1,Nv);
    isVertice(~isEdge) = true; % left vertices
    isVertice(~isEdge([Nv,1:Nv-1])) = true; % shift to right
    isLocalDof = [isEdge, isVertice];
    NdofElem(iel) = sum(isLocalDof);
    
    % ---------------- transition matrix D ------------------
    D1 = zeros(2*Nv, 3);   % edges, vertices
    v0e = v1(isEdge); v1e = v2(isEdge);
    D1(isLocalDof,:) = [ 0.5*( m(x(v0e),y(v0e)) + m(x(v1e),y(v1e)) ); ...
        m(x(isVertice),y(isVertice)); ];
    D1 = D1(isLocalDof,:);   D1 = blkdiag(D1,D1);
    D{iel} = D1;
    
    % -------------- integration over edges ------------------
    % C0
    C01x = zeros(1,2*Nv);  C01y = zeros(1,2*Nv);
    % B0 for constraints
    B0x = zeros(3,2*Nv);  B0y = zeros(3,2*Nv);
    % fK
    fK = zeros(2*Nv,1);
    f1int = integralTri(f1xy,3,nodeT,elemT);
    f2int = integralTri(f2xy,3,nodeT,elemT);
    for i = 1:Nv   % integrating basis functions over edges
        % basis for edges
        phie = zeros(1,Nv); phie(i) = 1*isEdge(i); % moment values of all basis functions on ei
        phinxe = Ne(i,1)*phie;  phinye = Ne(i,2)*phie;
        phitxe = Te(i,1)*phie;  phitye = Te(i,2)*phie;
        % basis for vertices
        phivm = zeros(1,Nv); % average is initialized as zero
        if ~isEdge(i)  % no moment on ei
            phivm(v1(i)) = 0.5*(1 + 0); 
            phivm(v2(i)) = 0.5*(0 + 1);
        end
        phinxvm = Ne(i,1)*phivm;  phinyvm = Ne(i,2)*phivm;
        phitxvm = Te(i,1)*phivm;  phityvm = Te(i,2)*phivm;
        % C0
        C01x = C01x + [phinxe, phinxvm];
        C01y = C01y + [phinye, phinyvm];
        % B0
        B0x(1,:) = B0x(1,:) + [phitxe, phitxvm];
        B0y(1,:) = B0y(1,:) + [phitye, phityvm];
        B0x(2,:) = B0x(2,:) + he(i)*[phie, phivm];
        % right-hand side
        fK = fK + [phie, phivm]';
    end
    
    % ------------------ H0,C0 ---------------------
    H0{iel} = area(iel);
    C01x = C01x(isLocalDof); C01y = C01y(isLocalDof);
    C01 = [C01x, C01y];
    C0{iel} = C01;
    
    % ---------------- B,Bs,G,Gs -------------------
    % B0
    B0y(3,:) = B0x(2,:); % second constraint - row 3
    B0 = [B0x(:,isLocalDof), B0y(:,isLocalDof)];
    % B
    E = zeros(6,4);
    E(2,1) = 1/hK;  E([3,5],[2,3]) = 1/(2*hK); E(6,4) = 1/hK;
    B1 = [E(:,1)*C01x+E(:,2)*C01y, E(:,3)*C01x+E(:,4)*C01y];
    % Bs
    B1s = B1;
    B1s([1,3,4],:) = B0;   Bs{iel} = B1s;
    % G,Gs
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;
    
    % ------------- load vector ------------------
    f1K = f1int/Nv*fK(isLocalDof);
    f2K = f2int/Nv*fK(isLocalDof);
    fK = [f1K; f2K];
    belem{iel} = fK'; % straighten as row vector for easy assembly
    
    % ------------- local-global index -------------
    indexDof = [indexEdge, index+NE];
    elem2 = idDof(indexDof(isLocalDof)); % connection number
    elem2dof{iel} = [elem2, elem2+NNdof];
    
end

%% Get elementwise stiffness matrix
ABelem = cell(NT,1);
Ph = cell(NT,1); % matrix for error evaluation
for iel = 1:NT
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));
    Pi0s = H0{iel}\C0{iel};
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi); AK = 2*para.mu*AK;
    BK = para.lambda*Pi0s'*H0{iel}*Pi0s;
    ABelem{iel} = reshape(AK'+BK',1,[]); % straighten
    % matrix for error evaluation
    Ph{iel} = Pis;
end

%% Assemble stiffness matrix and load vector
nnz = sum((2*NdofElem).^2);
dofNum = unique(NdofElem);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
id = 0; ff = zeros(NNdof2,1);
for Ndof = dofNum(:)' % only valid for row vector
    
    % assemble the matrix
    idv = find(NdofElem == Ndof); % find polygons with Nv vertices
    NTv = length(idv); % number of elements with Nv vertices
    
    elem2dofv = cell2mat(elem2dof(idv)); % elem2dof
    K = cell2mat(ABelem(idv)); F = cell2mat(belem(idv));
    Ndof2 = 2*Ndof;
    ii(id+1:id+NTv*Ndof2^2) = reshape(repmat(elem2dofv, Ndof2,1), [], 1);
    jj(id+1:id+NTv*Ndof2^2) = repmat(elem2dofv(:), Ndof2, 1);
    ss(id+1:id+NTv*Ndof2^2) = K(:);
    id = id + NTv*Ndof2^2;
    
    % assemble the vector
    ff = ff +  accumarray(elem2dofv(:),F(:),[NNdof2 1]);
end
kk = sparse(ii,jj,ss,NNdof2,NNdof2);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN;
if ~isempty(bdEdgeN)
    g_N = pde.g_N;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    Sig1 = g_N(z1); Sig2 = g_N(z2);
    F11 = sum(Ne.*Sig1(:,[1,3]),2)./2; F12 = sum(Ne.*Sig2(:,[1,3]),2)./2; % g1
    F21 = sum(Ne.*Sig1(:,[3,2]),2)./2; F22 = sum(Ne.*Sig2(:,[3,2]),2)./2;
    FN = [F11,F12,F21,F22];
    IdxN = idDof(bdEdgeN+NE); % connection number
    ff = ff + accumarray([IdxN(:); IdxN(:)+NNdof], FN(:),[NNdof2 1]);
end

%% Lagrange multiplier for pure traction problem
bdEdgeD = bdStruct.bdEdgeD; 
nd = 0; % initialized as zero
if isempty(bdEdgeD)
    % information of boundary
    nd = 3;
    bdEdge = bdStruct.bdEdge;
    z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:);
    Te = z2-z1;
    % d1,d2
    he = sqrt(sum((z2-z1).^2,2));
    F = 0.5*[he,he];
    id = idDof(bdEdge(:)+NE);  % connection number
    d1 = accumarray(id(:), F(:), [NNdof2 1]); 
    d2 = accumarray(id(:)+NNdof, F(:), [NNdof2 1]);
    % d3
    F1 = 0.5*[Te(:,1), Te(:,1)];
    F2 = 0.5*[Te(:,2), Te(:,2)];
    d3 = accumarray([id(:); id(:)+NNdof], [F1(:);F2(:)], [NNdof2 1]);
    % kkd, ffd
    kkd = sparse(NNdof2+nd, NNdof2+nd);
    kkd(1:NNdof2,1:NNdof2) = kk;
    kkd(1:NNdof2, NNdof2+(1:nd)) = [d1,d2,d3];
    kkd(NNdof2+(1:nd), 1:NNdof2) = [d1';d2';d3'];
    ffd = zeros(NNdof2+nd,1);
    ffd(1:NNdof2) = ff;
    % kk, ff
    kk = kkd;  ff = ffd;
end

%% Apply Dirichlet boundary conditions
bdEdgeD = bdStruct.bdEdgeD;
g_D = pde.g_D;
bdNodeIdxD = unique(bdEdgeD);
bdNodeIdxCont = idDof(bdNodeIdxD+NE); % connection number
isBdNode = false(NNdof2+nd,1);
isBdNode([bdNodeIdxCont, bdNodeIdxCont+NNdof]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdNodeIdxD,:);
u = zeros(NNdof2+nd,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
u = u(1:NNdof2);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk(1:NNdof2,1:NNdof2); %info.DofI = freeDof;
