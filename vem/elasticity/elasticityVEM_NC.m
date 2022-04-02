function [u,info] = elasticityVEM_NC(node,elem,pde,bdStruct)
%elasticityVEM_NC solves linear elasticity equation of tensor form using
% the standard nonconforming virtual element method in the local space of form (k = 1):
%
%    V_k(K) = { v \in H^1(K):  \Delta v = P_{k-2}(K),
%               grad(v)*n|_e \in P_{k-1}(e) for edges }
%
% The problem is
%
%       u = [u1, u2]
%       -div sigma(u) = f in \Omega,
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D,
%       Neumann boundary condition  sigma(u)n = [g1_N, g2_N] on \Gamma_N
%
% Copyright (C)  Terence Yu.

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
NT = size(elem,1); NE = size(edge,1);
NNdof = NE; NNdof2 = 2*NNdof;

%% Compute projection matrices
%% Also derive load vector and local-global index
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(2*elemLen);
elemf = zeros(nnz,1); Ff = zeros(nnz,1);
ia = 0; ib = 0;

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
    
    % ---------------- transition matrix D ------------------
    D = 0.5*( m(x(v1),y(v1)) + m(x(v2),y(v2)) );
    D = blkdiag(D,D);
    
    % ------------------ H0,C0 ---------------------
    H0 = area(iel);
    C01x = Ne(:,1)';  C01y = Ne(:,2)';
    C0 = [C01x, C01y];
    
    % ---------------- B,Bs,G,Gs -------------------
    % B0 for constraints
    B0x = zeros(3,Nv);    B0y = zeros(3,Nv);
    B0x(1,:) = Te(:,1)';  B0x(2,:) = he';
    B0y(2,:) = Te(:,2)';      
    B0y(3,:) = B0x(2,:); % second constraint - row 3
    B0 = [B0x, B0y];
    % B
    E = zeros(6,4);
    E(2,1) = 1/hK;  E([3,5],[2,3]) = 1/(2*hK); E(6,4) = 1/hK;
    B = [E(:,1)*C01x+E(:,2)*C01y, E(:,3)*C01x+E(:,4)*C01y];
    % Bs
    Bs = B;
    Bs([1,3,4],:) = B0;   
    % G,Gs
    G = B*D;     Gs = Bs*D;
    
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));
    Pi0s = H0\C0;
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi); 
    AK = 2*para.mu*AK;
    BK = para.lambda*Pi0s'*H0*Pi0s;
    AB = reshape(AK'+BK',1,[]);
    
    % ------------- local load vector ------------------
    fK = ones(Nv,1);
    f1int = integralTri(f1xy,3,nodeT,elemT);
    f2int = integralTri(f2xy,3,nodeT,elemT);
    f1K = f1int/Nv*fK;
    f2K = f2int/Nv*fK;
    fK = [f1K; f2K];
    
    % --------- assembly index for ellptic projection -----------
    indexDof = [indexEdge, indexEdge+NE];  Ndof = length(indexDof);
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemf(ib+1:ib+Ndof) = indexDof(:);
    Ff(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,NNdof2,NNdof2);
ff = accumarray(elemf,Ff,[NNdof2 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN; bdEdgeIdxN = bdStruct.bdEdgeIdxN;
if ~isempty(bdEdgeN)
    g_N = pde.g_N;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    Sig1 = g_N(z1); Sig2 = g_N(z2);
    F11 = sum(Ne.*Sig1(:,[1,3]),2);    F12 = sum(Ne.*Sig2(:,[1,3]),2); % g1 at endpoints
    F21 = sum(Ne.*Sig1(:,[3,2]),2);    F22 = sum(Ne.*Sig2(:,[3,2]),2);
    F1 = (F11+F12)/2;  F2 = (F21+F22)/2;
    FN = [F1, F2];
    ff = ff + accumarray([bdEdgeIdxN(:); bdEdgeIdxN(:)+NNdof], FN(:),[NNdof2 1]);
end

%% Lagrange multiplier for pure traction problem
bdEdgeD = bdStruct.bdEdgeD; 
nd = 0; % initialized as zero
if isempty(bdEdgeD)
    % information of boundary
    nd = 3;
    bdEdge = bdStruct.bdEdge;
    bdEdgeIdx = bdStruct.bdEdgeIdx;
    z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:);
    Te = z2-z1;
    % d1,d2
    d1 = zeros(NNdof2,1); d1(bdEdgeIdx) = 1;
    d2 = zeros(NNdof2,1); d2(bdEdgeIdx+NNdof) = 1;
    % d3
    d3 = zeros(NNdof2,1); 
    d3([bdEdgeIdx;bdEdgeIdx+NNdof]) = Te(:);
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
bdEdgeD = bdStruct.bdEdgeD;  bdEdgeIdxD = bdStruct.bdEdgeIdxD;
g_D = pde.g_D;
isBdNode = false(NNdof2+nd,1);
isBdNode([bdEdgeIdxD, bdEdgeIdxD+NNdof]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
z1 = node(bdEdgeD(:,1),:);  z2 = node(bdEdgeD(:,2),:);
uD = 0.5*(g_D(z1) + g_D(z2));
u = zeros(NNdof2+nd,1); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
u = u(1:NNdof2);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk(1:NNdof2,1:NNdof2); %info.DofI = freeDof;