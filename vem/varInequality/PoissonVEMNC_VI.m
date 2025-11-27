function [u,info] = PoissonVEMNC_VI(node,elem,pde,bdStruct)
%PoissonVEMNC_VI solves the simplified friction problem using nonconforming
% virtual element method in the local space of form (k = 1):
%
%    V_k(K) = { v \in H^1(K):  \Delta v = P_{k-2}(K),
%               grad(v)*n|_e \in P_{k-1}(e) for interior edges,
%               v|_e \in P_k(e) for boundary edges  }
%
% The problem is
%
%  D.E.
%     -\Delta u + cu = f,  in Omega,
%  Frictional boundary conditions
%     |grad(u)*n| <= g,  u*grad(u)*n + g*|u| = 0 on \Gamma_C
%  Dirichlet boundary conditions
%     u = g_D    on \Gamma_D
%  \partial(\Omega) = \Gamma_C + \Gamma_D
%
%   References
%   F. Feng, W. M. Han and J. G. Huang, "Virtual element methods for
%   elliptic variational inequalities of the second kind", J. Sci. Comput., 
%   Vol 80. No 1., pp. 60¨C80, 2019.  (See Eq. (4.6) there)
%
% Copyright (C)  Terence Yu, Li Jia and Fang Feng.

fxy = @(x,y) pde.f([x,y]);

%% Get auxiliary data
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge;  edge = auxT.edge;
% auxgeometry
aux = auxgeometry(node,elem);
centroid = aux.centroid;  diameter = aux.diameter;
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
NNdof = sum(isDof);  % number of d.o.f.s
% Connection number of d.o.f.s in form of edges and vertices
idDof = zeros(1,NE+N); % set the redundant as zero
idDof(isDof) = 1:NNdof;

%% Compute projection matrices
%% Also derive load vector and local-global index
D = cell(NT,1);  Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1);
H = cell(NT,1);
elem2dof = cell(NT,1); % local-global index
belem = cell(NT,1);  % load vector
NdofElem = zeros(NT,1); % numbers of local d.o.f.s on each cell
for iel = 1:NT
    
    % element information
    index = elem{iel};  indexEdge = elem2edge{iel};
    Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2);
    hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    % scaled monomials
    m1 = @(x,y) 1 + 0*x; m2 = @(x,y) (x-xK)./hK; m3 = @(x,y) (y-yK)./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]];
    % Ne, he
    v0 = 1:Nv;  v1 = [2:Nv,1];  % starting and ending
    Ne = [y(v1)-y(v0), x(v0)-x(v1)];
    he = sqrt(Ne(:,1).^2 + Ne(:,2).^2);
    
    % local logical vectors of d.o.f.s
    % isEdge, isVertice, isLocalDof
    isEdge = isDofEdge(indexEdge);
    isVertice = false(1,Nv);
    isVertice(~isEdge) = true; % left vertices
    isVertice(~isEdge([Nv,1:Nv-1])) = true; % shift to right
    isLocalDof = [isEdge, isVertice];
    NdofElem(iel) = sum(isLocalDof);
    
    % D
    D1 = zeros(2*Nv, 3);   % edges, vertices
    v0e = v0(isEdge); v1e = v1(isEdge);
    D1(isLocalDof,:) = [ 0.5*( m(x(v0e),y(v0e)) + m(x(v1e),y(v1e)) ); ...
        m(x(isVertice),y(isVertice)); ];
    D1 = D1(isLocalDof,:);
    D{iel} = D1;
    
    % B, Bs, fK
    B1 = zeros(3,2*Nv);  % initialized as zero
    B1s = zeros(3,2*Nv);
    fK = zeros(2*Nv,1);
    fint = integralTri(fxy,3,nodeT,elemT);
    for i = 1:Nv
        % gradm*Ne
        gi = sum(Gradm.*repmat(Ne(i,:),3,1), 2);
        % for edges
        phie = zeros(1,Nv); phie(i) = 1*isEdge(i);
        Gradmphie =  gi*phie;
        % for vertices
        phivm = zeros(1,Nv); % average is initialized as zero
        if ~isEdge(i)  % no moment on ei
            phivm(v0(i)) = 0.5*(1 + 0);  % i-th base
            phivm(v1(i)) = 0.5*(0 + 1); % (i+1)-th base
        end
        Gradmphiv = gi*phivm;
        B1 = B1 + [Gradmphie, Gradmphiv];
        % add constraint
        B1s(1,:) = B1s(1,:) + he(i)*[phie, phivm];
        % right-hand side
        fK = fK + [phie, phivm]';
    end
    B1s(2:end,:) = B1(2:end,:);
    B1 = B1(:,isLocalDof);  B1s = B1s(:,isLocalDof);
    Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;
    fK = fint/Nv*fK(isLocalDof);
    belem{iel} = fK'; % straighten as row vector for easy assembly
    
    % H
    mm = @(x,y) [m1(x,y).*m1(x,y), m1(x,y).*m2(x,y), m1(x,y).*m3(x,y), ...
        m2(x,y).*m1(x,y), m2(x,y).*m2(x,y), m2(x,y).*m3(x,y), ...
        m3(x,y).*m1(x,y), m3(x,y).*m2(x,y), m3(x,y).*m3(x,y)];
    H1 = zeros(3,3);
    H1(:) = integralTri(mm,3,nodeT,elemT); % n = 3
    H{iel} = H1;
    
    % elem2dof
    indexDof = [indexEdge, index+NE];
    elem2dof{iel} = idDof(indexDof(isLocalDof)); % connection number
end

%% Get elementwise stiffness matrix
ABelem = cell(NT,1);
Ph = cell(NT,1); % matrix for error evaluation
for iel = 1:NT
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pis'*H{iel}*Pis + pde.c*hK^2*(I-Pi)'*(I-Pi);  % Pis = Pi0s;
    ABelem{iel} = reshape(AK'+BK',1,[]); % straighten as row vector for easy assembly    
    % matrix for error evaluation
    Ph{iel} = Pis;
end

%% Assemble stiffness matrix and load vector
nnz = sum(NdofElem.^2); dofNum = unique(NdofElem);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
id = 0; ff = zeros(NNdof,1);
for Ndof = dofNum(:)' % only valid for row vector
    
    % assemble the matrix
    idv = find(NdofElem == Ndof); % find polygons with Nv vertices
    NTv = length(idv); % number of elements with Nv vertices
    
    elem2dofv = cell2mat(elem2dof(idv)); % elem2dof
    K = cell2mat(ABelem(idv)); F = cell2mat(belem(idv));
    ii(id+1:id+NTv*Ndof^2) = reshape(repmat(elem2dofv, Ndof,1), [], 1);
    jj(id+1:id+NTv*Ndof^2) = repmat(elem2dofv(:), Ndof, 1);
    ss(id+1:id+NTv*Ndof^2) = K(:);
    id = id + NTv*Ndof^2;
    
    % assemble the vector
    ff = ff +  accumarray(elem2dofv(:),F(:),[NNdof 1]);
end
kk = sparse(ii,jj,ss,NNdof,NNdof);

%% Assemble frictional boundary conditions
bdEdgeFri = bdStruct.bdEdgeN;
if ~isempty(bdEdgeFri)
    g = pde.g_C;
    z1 = node(bdEdgeFri(:,1),:); z2 = node(bdEdgeFri(:,2),:);
    gC = max([g(z1); g(z2)]);
    e = z2-z1;
    he = sqrt(e(:,1).^2+e(:,2).^2);
    Fri = 0.5*gC*[he,he];
    IdxFri = idDof(bdEdgeFri+NE); % connection number
    w = accumarray(IdxFri(:), Fri(:),[NNdof 1]);
end

%% Apply Dirichlet boundary conditions
bdEdgeD = bdStruct.bdEdgeD;
g_D = pde.g_D;
bdNodeIdxD = unique(bdEdgeD);
bdNodeIdxCont = idDof(bdNodeIdxD+NE); % connection number
isBdNode = false(NNdof,1);
isBdNode(bdNodeIdxCont) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdNodeIdxD,:);
u = zeros(NNdof,1); uD = g_D(nodeD); u(bdDof) = uD(:);

%% Get minimization problem
%fun = @(v)  0.5*v'*A*v - b'*v + w'*abs(v);
A = speye(NNdof,NNdof);
A(freeDof,freeDof) = kk(freeDof,freeDof);
b = u;
b(freeDof) = ff(freeDof) - kk(freeDof,freeDof)*u(freeDof);

%% Set solver
% Note: It can be directly solved by: u = cvx2cvx(A,b,w)
% schur 
idFri = unique(bdEdgeFri(:)) + NE;
idFri = idDof(idFri);  idFri = idFri(:);  % connection number of frictional nodes
idFriRes = setdiff(1:NNdof,idFri); % connection number of other nodes
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
u = zeros(NNdof,1);
u(idFri) = v;
u(idFriRes) = vRes;
u(bdDof) = uD(:);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; info.freeDof = freeDof;
info.isDof = isDof;