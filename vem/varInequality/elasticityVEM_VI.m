function [u,info] = elasticityVEM_VI(node,elem,pde,bdStruct)
%elasticityVEM_VI solves linear elasticity equation of tensor form using
% conforming virtual element method in the lowest order case.
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

%% PDE data and boundary informaiton
para = pde.para;
bdEdgeType = bdStruct.bdEdgeType;
bdEdgeC = bdEdgeType{1}; 
bdEdgeD = bdEdgeType{2}; 
bdEdgeN = [];
if length(bdEdgeType)==3
    bdEdgeN = bdEdgeType{3};
end

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
% number
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
    
    % ---------------- D ------------------
    m1 = @(x,y) 1 + 0*x; m2 = @(x,y) (x-xK)./hK; m3 = @(x,y) (y-yK)./hK;
    D1 = [m1(x,y), m2(x,y), m3(x,y)];   
    D1 = blkdiag(D1,D1);  
    D{iel} = D1;
    
    % ---------------- H0,C0 ---------------- 
    id = [Nv,1:Nv-1];
    phinx = 0.5*(Ne(id,1) + Ne(:,1))';
    phiny = 0.5*(Ne(id,2) + Ne(:,2))';
    H0{iel} = area(iel); 
    C01 = [phinx, phiny];
    C0{iel} = C01;
    
    % ---------------- B,Bs,G,Gs -------------------
    % E = [E11,E12,E21,E22]
    E = zeros(6,4);
    E(2,1) = 1/hK;  E([3,5],[2,3]) = 1/(2*hK); E(6,4) = 1/hK;
    % B1
    B1 = [E(:,1)*phinx+E(:,2)*phiny, E(:,3)*phinx+E(:,4)*phiny];
    
    % B0 for constraint
    B0 = zeros(3,2*Nv);    
    B0(1,1:Nv) = 1;     B0(2,Nv+1:end) = 1;
    B0(3,1:Nv) = -y;    B0(3,Nv+1:end) = x;
   
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
    fK = pde.f1(centroid(iel,:))*area(iel)/Nv; 
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
if ~isempty(bdEdgeN)
    g_N = pde.g_N;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    Sig1 = g_N(z1); Sig2 = g_N(z2);
    F11 = sum(Ne.*Sig1(:,[1,3]),2)./2; F12 = sum(Ne.*Sig2(:,[1,3]),2)./2; % g1
    F21 = sum(Ne.*Sig1(:,[3,2]),2)./2; F22 = sum(Ne.*Sig2(:,[3,2]),2)./2;
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
info.kk = kk; info.freeDof = freeDof;
