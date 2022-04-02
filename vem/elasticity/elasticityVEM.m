function [u,info] = elasticityVEM(node,elem,pde,bdStruct,constrainttype)
%elasticityVEM solves linear elasticity equation of tensor form using 
% virtual element method in V1
%
%       u = [u1, u2]
%       -div sigma(u) = f in \Omega,
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D,
%       Neumann boundary condition  sigma(u)n = [g1_N, g2_N] on \Gamma_N
%
% Copyright (C)  Terence Yu.

%% Input check
if ~exist('constrainttype','var') || isempty(constrainttype), constrainttype = 1; end

%% Get auxiliary data
para = pde.para;
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
N = size(node,1);  NT = size(elem,1);

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
    % ------------- element information -------------
    index = elem{iel};     Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    Te = [-Ne(:,2), Ne(:,1)];
    
    % --------------- scaled monomials --------------
    m1 = @(x,y) 1 + 0*x; 
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    
    % ---------------- transition matrix ------------------    
    D = m(x,y);   
    D = blkdiag(D,D);  
    
    % --------- L2 projection -----------
    H0 = area(iel);
    elem1 = [v1(:), v2(:)];
    C0 = zeros(1,2*Nv);
    F = 1/2*[(1*Ne); (1*Ne)]; % [he*n1, he*n2]
    C0(:) = accumarray([elem1(:);elem1(:)+Nv], F(:), [2*Nv 1]);
    
    % --------- elliptic projection -----------
    % E = [E11,E12,E21,E22]
    E = zeros(6,4);
    E(2,1) = 1/hK;  E([3,5],[2,3]) = 1/(2*hK); E(6,4) = 1/hK;
    % B
    B = [E(:,1)*C0(1:Nv)+E(:,2)*C0(Nv+1:end), ...
          E(:,3)*C0(1:Nv)+E(:,4)*C0(Nv+1:end)];
    
    % B0 for constraint
    B0 = zeros(3,2*Nv);
    switch constrainttype
        case 1  % sum_i (v(zi),p(zi))
            B0(1,1:Nv) = 1;     B0(2,Nv+1:end) = 1;
            B0(3,1:Nv) = -y;    B0(3,Nv+1:end) = x;
        case 2  % int_{\partial K} v.p ds
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
    Bs = B; Bs([1,3,4],:) = B0;   
    % consistency relation
    G = B*D;     Gs = Bs*D;
    
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));
    Pi0s = H0\C0;
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi); 
    AK = 2*para.mu*AK;
    BK = para.lambda*Pi0s'*H0*Pi0s;
    AB = reshape(AK'+BK',1,[]);
    
    % --------- local load vector -----------
    fK = pde.f(centroid(iel,:))*area(iel)/Nv; 
    fK = repmat(fK,Nv,1);
    
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
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,2*N,2*N);
ff = accumarray(elemf,Ff,[2*N 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN;
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
    d1 = accumarray(bdEdge(:), F(:), [2*N 1]);
    d2 = accumarray(bdEdge(:)+N, F(:), [2*N 1]);
    % d3
    F1 = 0.5*[Te(:,1), Te(:,1)];
    F2 = 0.5*[Te(:,2), Te(:,2)];
    d3 = accumarray([bdEdge(:); bdEdge(:)+N], [F1(:);F2(:)], [2*N 1]);
    % kkd, ffd
    kkd = sparse(2*N+3, 2*N+3);
    kkd(1:2*N,1:2*N) = kk;
    kkd(1:2*N, 2*N+(1:nd)) = [d1,d2,d3];
    kkd(2*N+(1:nd), 1:2*N) = [d1';d2';d3'];
    ffd = zeros(2*N+nd,1);
    ffd(1:2*N) = ff;
    % kk, ff
    kk = kkd;  ff = ffd;
end

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;
bdNodeIdxD = bdStruct.bdNodeIdxD;
isBdNode = false(2*N+nd,1); isBdNode([bdNodeIdxD;bdNodeIdxD+N]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdNodeIdxD,:);
u = zeros(2*N+nd,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
u = u(1:2*N);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk(1:2*N,1:2*N);  %info.DofI = freeDof;