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
    % element information
    index = elem{iel};     Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    
    % D
    m1 = @(x,y) 1 + 0*x; m2 = @(x,y) (x-xK)./hK; m3 = @(x,y) (y-yK)./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    D1 = m(x,y);   D1 = blkdiag(D1,D1);  D{iel} = D1;
    
    % H0,C0
    rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
    normVec = 0.5*[y(rotid2) - y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
    H0{iel} = area(iel); C01 = reshape(normVec',1,[]);
    C0{iel} = C01;
    
    % ------- B,Bs,G,Gs ----------
    % E = [E11,E12,E21,E22]
    E = zeros(6,4);
    E(2,1) = 1/hK;  E([3,5],[2,3]) = 1/hK; E(6,4) = 1/hK;
    % phinx, phiny
    nx = ne(:,1); ny = ne(:,2);
    id = [Nv,1:Nv-1];
    phinx = 0.5*(nx(id).*he(id) + nx.*he)';
    phiny = 0.5*(ny(id).*he(id) + ny.*he)';
    % B1
    B1 = [E(:,1)*phinx+E(:,2)*phiny, E(:,3)*phinx+E(:,4)*phiny];
    % constraint
    B0 = zeros(3,2*Nv);
    switch constrainttype
        case 1  % sum_i (v(zi),p(zi))
            B0(1,1:Nv) = 1; B0(2,Nv+1:end) = 1;
            B0(3,1:Nv) = -y;    B0(3,Nv+1:end) = x;
        case 2  % int_{\partial K} v.p ds
            base = eye(Nv);
            phi = zeros(1,Nv); phix = zeros(1,Nv); phiy = zeros(1,Nv);
            for i = 1:Nv % loop of edges
                base1 = base(i,v1); base2 = base(i,v2); basec = (base1+base2)./2;
                x1 = x(v1(i)); x2 = x(v2(i)); xc = (x1+x2)/2;
                y1 = y(v1(i)); y2 = y(v2(i)); yc = (y1+y2)/2;
                phi = phi + he(i)*basec;
                phix = phix + he(i)/6*(base1*x1 + 4*basec*xc + base2*x2);
                phiy = phiy + he(i)/6*(base1*y1 + 4*basec*yc + base2*y2);
            end
            B0(1,1:Nv) = phi; B0(2,Nv+1:end) = phi;
            B0(3,:) = [-phiy, phix];
        case 3 % int_K \nabla \times v dx,  int_{\partial K} v ds
            % first constraint
            te = [-ne(:,2), ne(:,1)];
            tx = te(:,1);  ty = te(:,2);
            phitx = 0.5*(tx(id).*he(id) + tx.*he)';
            phity = 0.5*(ty(id).*he(id) + ty.*he)';
            B0(1,:) = [phitx, phity];
            % second constraint
            phi = 0.5*(he(id) + he)';
            B0(2,1:Nv) = phi;  B0(3,Nv+1:end) = phi;
    end    
    % B{iel} = B1;
    B1s = B1; B1s([1,3,4],:) = B0;   Bs{iel} = B1s;
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
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi); AK = 2*para.mu*AK;
    BK = para.lambda*Pi0s'*H0{iel}*Pi0s;
    ABelem{iel} = reshape(AK'+BK',1,[]); % straighten
    % Load vector
    Nv = length(elem{iel});
    fK = pde.f(centroid(iel,:))*area(iel)/Nv; fK = repmat(fK,Nv,1);
    belem{iel} = fK(:)'; % straighten
    % matrix for L2 and H1 error evaluation
    Ph{iel} = Pis;
end

%% Assemble stiffness matrix and load vector
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnz = sum(4*elemLen.^2);
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
    ii(id+1:id+NTv*Ndof2^2) = reshape(repmat(elem2, Ndof2,1), [], 1);
    jj(id+1:id+NTv*Ndof2^2) = repmat(elem2(:), Ndof2, 1);
    ss(id+1:id+NTv*Ndof2^2) = K(:);
    id = id + NTv*Ndof2^2;
    
    % assemble the vector
    ff = ff +  accumarray(elem2(:),F(:),[2*N 1]);
    
    % elementwise global indices
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
    Sig1 = g_N(z1); Sig2 = g_N(z2);
    F11 = sum(ne.*Sig1(:,[1,3]),2)./2; F12 = sum(ne.*Sig2(:,[1,3]),2)./2; % g1
    F21 = sum(ne.*Sig1(:,[3,2]),2)./2; F22 = sum(ne.*Sig2(:,[3,2]),2)./2;
    FN = [F11,F12,F21,F22];
    ff = ff + accumarray([bdEdgeN(:); bdEdgeN(:)+N], FN(:),[2*N 1]);
end

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;
bdNodeIdx = bdStruct.bdNodeIdx;
isBdNode = false(N,1); isBdNode(bdNodeIdx) = true;
bdNode = find(isBdNode); freeNode = find(~isBdNode);
nodeD = node(bdNode,:);
bdDof = [bdNode; bdNode+N]; freeDof = [freeNode;freeNode+N];
u = zeros(2*N,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
solver = 'amg';
if 2*N < 2e3, solver = 'direct'; end
% solve
switch solver
    case 'direct'
        u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
    case 'amg'
        option.solver = 'CG';
        u(freeDof) = amg(kk(freeDof,freeDof),ff(freeDof),option);
end

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; info.freeDof = freeDof;