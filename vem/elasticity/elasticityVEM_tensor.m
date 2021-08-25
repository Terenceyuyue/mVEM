function [u,info] = elasticityVEM_tensor(node,elem,pde,bdStruct,constrainttype)
%elasticityVEM_tensor solves linear elasticity equation of tensor
%form using virtual element method in V1
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

%% Compute projection matrices
D = cell(NT,1);
% B = cell(NT,1);  % not used in the computation
Bs = cell(NT,1);
G = cell(NT,1); Gs = cell(NT,1);
H0 = cell(NT,1); C0 = cell(NT,1);
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
    D1 = m(x,y);   
    D1 = blkdiag(D1,D1);  
    D{iel} = D1;
    
    % --------- L2 projection -----------
    H0{iel} = area(iel);
    elem1 = [v1(:), v2(:)];
    C01 = zeros(1,2*Nv);
    F = 1/2*[(1*Ne); (1*Ne)]; % [n1, n2]
    C01(:) = accumarray([elem1(:);elem1(:)+Nv], F(:), [2*Nv 1]);
    C0{iel} = C01;
    
    % --------- elliptic projection -----------
    % E = [E11,E12,E21,E22]
    E = zeros(6,4);
    E(2,1) = 1/hK;  E([3,5],[2,3]) = 1/(2*hK); E(6,4) = 1/hK;
    % B1
    B1 = [E(:,1)*C01(1:Nv)+E(:,2)*C01(Nv+1:end), ...
          E(:,3)*C01(1:Nv)+E(:,4)*C01(Nv+1:end)];
    
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
    B1s = B1; B1s([1,3,4],:) = B0;   
    Bs{iel} = B1s;
    % G,Gs
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;
 end

%% Get elementwise stiffness matrix and load vector
ABelem = cell(NT,1); belem = cell(NT,1);
Ph = cell(NT,1); % matrix for error evaluation
for iel = 1:NT
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));
    Pi0s = H0{iel}\C0{iel};
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi); 
    AK = 2*para.mu*AK;
    BK = para.lambda*Pi0s'*H0{iel}*Pi0s;
    ABelem{iel} = reshape(AK'+BK',1,[]); % straighten
    % Load vector
    Nv = length(elem{iel});
    fK = pde.f(centroid(iel,:))*area(iel)/Nv; 
    fK = repmat(fK,Nv,1);
    belem{iel} = fK(:)'; % straighten
    % matrix for L2 and H1 error evaluation
    Ph{iel} = Pis;
end

%% Assemble stiffness matrix and load vector
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnz = sum((2*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
id = 0; ff = zeros(2*N,1);
elem2dof = cell(NT,1);
for Nv = vertNum(:)' % only valid for row vector
    
    % assemble the matrix
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices
    
    elemNv = cell2mat(elem(idNv)); % elem
    K = cell2mat(ABelem(idNv)); F = cell2mat(belem(idNv));
    elem2 = [elemNv, elemNv+N];
    Ndof = Nv; Ndof2 = 2*Ndof;
    ii(id+1:id+NTv*Ndof2^2) = reshape(repmat(elem2, Ndof2, 1), [], 1);
    jj(id+1:id+NTv*Ndof2^2) = repmat(elem2(:), Ndof2, 1);
    ss(id+1:id+NTv*Ndof2^2) = K(:);
    id = id + NTv*Ndof2^2;
    
    % assemble the vector
    ff = ff +  accumarray(elem2(:),F(:),[2*N 1]);
    
    % elementwise global index
    elem2dof(idNv) = mat2cell(elem2, ones(NTv,1), Ndof2);
end
kk = sparse(ii,jj,ss,2*N,2*N);

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

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;
bdNodeIdx = bdStruct.bdNodeIdx;
isBdNode = false(2*N,1); isBdNode([bdNodeIdx;bdNodeIdx+N]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdNodeIdx,:);
u = zeros(2*N,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; %info.DofI = freeDof;