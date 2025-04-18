function [u,info] = elasticityVEM_reducedIntegration_VI(node,elem,pde,bdString,constraintType,refineType)
%elasticityVEM_reducedIntegration_VI solves linear elasticity equation of tensor form using
% conforming virtual element method based on the reduced integration technique in the lowest 
% order case.
%
% The problem is
%
%       u = [u1, u2]
%       -div sigma(u) = f in \Omega,
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D,
%       Neumann boundary condition  sigma(u)n = [g1_N, g2_N] on \Gamma_N
%
%  Frictional boundary condition: Neumann is replaced by this condition
%
% Copyright (C)  Terence Yu.

para = pde.para;
%% Information of coarse mesh
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge;  
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
% numbers
N = size(node,1);  NT = size(elem,1);  NE = size(edge,1); 
switch refineType
    case 1   % center --- midpoints
        NNdof = N+NE+NT;   
    case 2   % connecting midpoints
        NNdof = N+NE;
    case 3   % two edges for hanging node 
        NNdof = N+NE;
end

%% Assemble L2 div projection on coarse mesh
elemLen = cellfun('length',elem); 
nnz = sum((2*(2*elemLen)).^2);
ii0 = zeros(nnz,1); jj0 = zeros(nnz,1); ss0 = zeros(nnz,1);
ia = 0;
for iel = 1:NT
    % element information
    index = elem{iel};  indexEdge = auxT.elem2edge{iel}; 
    Nv = length(index);
    x = node(index,1); y = node(index,2); %size:Nv*1
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    
    % L2 div projection
    H0K = aux.area(iel);
    C0x = zeros(1,2*Nv);  C0y = zeros(1,2*Nv);
    p1 = 1:Nv; p2 = [Nv,1:Nv-1];
    C0x(p1) = 1/4*(Ne(p2,1) + Ne(p1,1)); % vertices
    C0x(p1+Nv) = 1/2*Ne(:,1);            % mid-points
    C0y(p1) = 1/4*(Ne(p2,2) + Ne(p1,2));
    C0y(p1+Nv) = 1/2*Ne(:,2);
    C0K = [C0x, C0y];     
    Pi0sK = H0K\C0K; 
    BK = Pi0sK'*H0K*Pi0sK;
    
    % local to global index
    Ndof = 2*Nv; Ndof2 = 2*Ndof;
    indexDof = [index, indexEdge+N, index+NNdof, indexEdge+N+NNdof];

    % assembly index
    ii0(ia+1:ia+Ndof2^2) = reshape(repmat(indexDof, Ndof2, 1), [], 1);
    jj0(ia+1:ia+Ndof2^2) = repmat(indexDof(:), Ndof2, 1);
    ss0(ia+1:ia+Ndof2^2) = reshape(BK',[],1);
    ia = ia + Ndof2^2;
end

%% Information of the refined mesh
% refine the mesh
[node,elem] = PolyMeshUniformRefine(node,elem,refineType);
% boundary information
bdStruct = setboundaryPiecewise(node,elem,bdString);
bdEdgeType = bdStruct.bdEdgeType;
bdEdgeC = bdEdgeType{1}; 
bdEdgeD = bdEdgeType{2}; 
bdEdgeN = [];
if length(bdEdgeType)==3
    bdEdgeN = bdEdgeType{3};
end
% auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
N = size(node,1);  NT = size(elem,1);

%% Assemble elliptic projection on the refined mesh
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum((2*elemLen));
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % ------------- element information -------------
    index = elem{iel}; Nv = length(index);
    xT = centroid(iel,1); yT = centroid(iel,2); hT = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    Te = [-Ne(:,2), Ne(:,1)];
    elem1 = [v1(:), v2(:)]; 
    
    % --------------- scaled monomials --------------
    m1 = @(x,y) 1 + 0*x;
    m2 = @(x,y) (x-xT)./hT;
    m3 = @(x,y) (y-yT)./hT;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    
    % ---------------- transition matrix ------------------
    D1 = m(x,y);
    D = blkdiag(D1,D1);
    
    % --------- L2 projection -----------    
    p2 = [Nv,1:Nv-1];
    phinx = 0.5*(Ne(p2,1) + Ne(:,1))';
    phiny = 0.5*(Ne(p2,2) + Ne(:,2))';
    C01 = [phinx, phiny];
    
    % --------- elliptic projection -----------
    % E = [E11,E12,E21,E22]
    E = zeros(6,4);
    E(2,1) = 1/hT;  E([3,5],[2,3]) = 1/(2*hT); E(6,4) = 1/hT;
    % B
    B = [E(:,1)*C01(1:Nv)+E(:,2)*C01(Nv+1:end), ...
        E(:,3)*C01(1:Nv)+E(:,4)*C01(Nv+1:end)];
    
    % B0 for constraint
    B0 = zeros(3,2*Nv);
    switch constraintType
        case 1  % \sum_{i=1}^{Nv} (\boldsymbol{\phi}(zi),\boldsymbol{p}(zi))
            B0(1,1:Nv) = 1;     B0(2,Nv+1:end) = 1;
            B0(3,1:Nv) = -y;    B0(3,Nv+1:end) = x;
        case 2  % \int_{\partial K} \boldsymbol{\phi} \cdot \boldsymbol{p} ds
            F = 1/2*ones(2*Nv,1); % v1,v2
            phi = accumarray(elem1(:), F(:), [Nv 1]);
            F = 1/2*[x(v1); x(v2)];
            phix = accumarray(elem1(:), F(:), [Nv 1]);
            F = 1/2*[y(v1); y(v2)];
            phiy = accumarray(elem1(:), F(:), [Nv 1]);
            B0(1,1:Nv) = phi; B0(2,Nv+1:end) = phi;
            B0(3,:) = [-phiy', phix'];
        case 3 % int_K \nabla \times v dx,  int_{\partial K} v ds
            % first constraint
            F = 1/2*[(1*Te); (1*Te)]; % [t1, t2]
            B0(1,:) = accumarray([elem1(:);elem1(:)+Nv], F(:), [2*Nv 1]);
            % second constraint
            F = 1/2*ones(2*Nv,1); % v1,v2
            phi = accumarray(elem1(:), F(:), [Nv 1]);
            B0(2,1:Nv) = phi;  B0(3,Nv+1:end) = phi;
    end
    
    % Bs
    Bs = B;   Bs([1,3,4],:) = B0;
    % G,Gs
    G = B*D;     Gs = Bs*D;    
    % local stiffness matrix
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));
    AT  = Pis'*G*Pis + (I-Pi)'*(I-Pi);   
    
    % --------- load vector -----------
    fT = pde.f1(centroid(iel,:))*area(iel)/Nv;
    fT = repmat(fT,Nv,1);
    
    % --------- assembly index for ellptic projection -----------
    indexDof = [index, index+N];  % local to global index
    Ndof = Nv; Ndof2 = 2*Ndof;
    ii(ia+1:ia+Ndof2^2) = reshape(repmat(indexDof, Ndof2, 1), [], 1);
    jj(ia+1:ia+Ndof2^2) = repmat(indexDof(:), Ndof2, 1);
    ss(ia+1:ia+Ndof2^2) = reshape(AT',1,[]);
    ia = ia + Ndof2^2;
    
    % --------- assembly index for right hand side -----------
    elemb(ib+1:ib+Ndof2) = indexDof(:);
    Fb(ib+1:ib+Ndof2) = fT(:);
    ib = ib + Ndof2;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kk = sparse([ii;ii0], [jj;jj0], [2*para.mu*ss; para.lambda*ss0], 2*N,2*N);
ff = accumarray(elemb,Fb,[2*N 1]);

%% Assemble Neumann boundary conditions
if ~isempty(bdEdgeN)
    g_N = pde.g_N;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    ne = [-e(:,2),e(:,1)]; % scaled ne
    Sig1 = g_N(z1); Sig2 = g_N(z2);    % Sig = [sig11,sig22,sig12]
    F11 = sum(ne.*Sig1(:,[1,3]),2)./2; % g1
    F12 = sum(ne.*Sig2(:,[1,3]),2)./2; % g1
    F21 = sum(ne.*Sig1(:,[3,2]),2)./2;
    F22 = sum(ne.*Sig2(:,[3,2]),2)./2;
    FN = [F11,F12,F21,F22];
    ff = ff + accumarray([bdEdgeN(:); bdEdgeN(:)+N], FN(:),[2*N 1]);
end

%% Assemble frictional boundary conditions
g = pde.g_C;
z1 = node(bdEdgeC(:,1),:); z2 = node(bdEdgeC(:,2),:);
gC = [g(z1);g(z2)];  gC = max(gC);
e = z2-z1;
he = sqrt(e(:,1).^2+e(:,2).^2);
Fri = 0.5*gC*[he,he];
w1 = accumarray(bdEdgeC(:), Fri(:),[2*N 1]);
w2 = zeros(size(w1));
w = [w1; w2];

%% Apply Dirichlet boundary conditions
% u = gD on \Gamma_D
% u_\tau = 0 on \Gamma_C
g_D = pde.g_D;
bdNodeIdxD = unique(bdEdgeD);
bdNodeIdxFri = unique(bdEdgeC);
bdNodeIdx = unique([bdNodeIdxD; bdNodeIdxFri]);
isBdNode = false(2*N,1);
isBdNode([bdNodeIdx, bdNodeIdx+N]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeIdx = node(bdNodeIdx,:);
u = zeros(2*N,1); uIdx = g_D(nodeIdx); u(bdDof) = uIdx(:);

%% Get minimization problem
%fun = @(v)  0.5*v'*A*v - b'*v + w'*abs(v);
A = speye(2*N,2*N);
A(freeDof,freeDof) = kk(freeDof,freeDof);
b = u;
b(freeDof) = ff(freeDof) - kk(freeDof,freeDof)*u(freeDof);

%% Set solver
% schur 
idFri = unique(bdEdgeC(:));
idFri = [idFri; idFri+N];
idFriRes = setdiff(1:2*N,idFri); 
idFriRes = idFriRes(:);
A11 = A(idFri,idFri);  
A12 = A(idFri,idFriRes);
A22 = A(idFriRes,idFriRes);
b1 = b(idFri);  b2 = b(idFriRes);
AA = A11-A12*(A22\A12');
bb = b1-A12*(A22\b2);
ww = w(idFri);
% Fixed point algorithm based on proximity operator
v = cvx2cvx(AA,bb,ww);
vRes = A22\(b2-A12'*v);
u = zeros(2*N,1);
u(idFri) = v;
u(idFriRes) = vRes;

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; %info.DofI = freeDof;
info.node = node; info.elem = elem;
