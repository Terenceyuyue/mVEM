function [u,info] = elasticityVEM_Navier(node,elem,pde,bdStruct)
%elasticityVEM_Navier solves linear elasticity equation of Navier form 
% using the conforming virtual element method in lowest order case
%
%       u = [u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%
% Copyright (C)  Terence Yu.

%% Get auxiliary data
para = pde.para;
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
N = size(node,1);  NT = size(elem,1);
Nm = 3;

%% Compute and assemble linear system 
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(2*elemLen);
elemf = zeros(nnz,1); Ff = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % -------- element information ----------
    index = elem{iel};     Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1];
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    
    % -------- scaled monomials ----------
    m1 = @(x,y) 1+0*x;                gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;           gradm2 = @(x,y) [1+0*x, 0+0*x]./hK;
    m3 = @(x,y) (y-yK)./hK;           gradm3 = @(x,y) [0+0*x, 1+0*x]./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    Gradmc = {gradm1, gradm2, gradm3};
    
    % ---------- transition matrix ------------
    D = m(x,y);
    
    % --------- elliptic projection -----------
    % first term  = 0
    B = zeros(Nm, Nv);
    % second term
    elem1 = [v1(:), v2(:)];
    for im = 1:Nm
        gradmc = Gradmc{im};
        F1 = 0.5*sum(gradmc(x(v1), y(v1)).*Ne, 2);
        F2 = 0.5*sum(gradmc(x(v2), y(v2)).*Ne, 2);
        F = [F1, F2];
        B(im, :) = accumarray(elem1(:), F(:), [Nv 1]);
    end
    % constraint
    Bs = B; Bs(1,:) = 1/Nv;
    % consistency relation
    G = B*D;     Gs = Bs*D;
    
    % --------- L2 projection -----------
    H0 = area(iel);
    C0 = zeros(1,2*Nv);
    F = 1/2*[(1*Ne); (1*Ne)]; % [n1, n2]
    C0(:) = accumarray([elem1(:);elem1(:)+Nv], F(:), [2*Nv 1]);
    
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));
    Pi0s = H0\C0;
    % Stiffness matrix
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi); AK = para.mu*blkdiag(AK,AK);
    BK = (para.mu+para.lambda)*Pi0s'*H0*Pi0s;
    AB = reshape(AK'+BK',1,[]); % straighten
    
     % --------- load vector -----------
     fK = pde.f(centroid(iel,:))*area(iel)/Nv; fK = repmat(fK,Nv,1);

     % --------- assembly index for ellptic projection -----------
    indexDof = [index, index+N];  Ndof = length(indexDof);
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemf(ib+1:ib+Ndof) = indexDof(:);
    Ff(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = blkdiag(Pis,Pis);
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,2*N,2*N);
ff = accumarray(elemf,Ff,[2*N 1]);

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;
bdNodeIdxD = bdStruct.bdNodeIdxD;
isBdNode = false(2*N,1); isBdNode([bdNodeIdxD;bdNodeIdxD+N]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdNodeIdxD,:);
u = zeros(2*N,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; %info.DofI = freeDof;