function [u,info] = PoissonVEM_vec(node,elem,pde,bdStruct)
%PoissonVEM_vec solves Poisson equation using virtual element method in V1
%with a vectorized implementation.
%
%     -\Delta u + cu = f,  in Omega
%     Dirichlet boundary condition u=g_D on \Gamma_D,
%     Neumann boundary condition   grad(u)*n=g_N on \Gamma_N
%
% Copyright (C)  Terence Yu.

%% Get auxiliary data
% auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
% number
N = size(node,1);  NT = size(elem,1);
elemLen = cellfun('length',elem);
vertNum = unique(elemLen);

%% Triangulation
nodeTri = [node; aux.centroid];
% elemTri: [a_i, a_{i+1}, centroid]

%% Quadrature points
[lambda,weight] = quadpts(2);

%% Compute projection matrices
Dc = cell(NT,1); 
Bc = cell(NT,1);  Bsc = cell(NT,1);
Hc = cell(NT,1);
for Nv = vertNum(:)'
    % information of polygons with Nv vertices
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv);  % number of elements with Nv vertices
    elemNv = cell2mat(elem(idNv)); % elem
    x = reshape(node(elemNv,1),NTv,Nv);
    y = reshape(node(elemNv,2),NTv,Nv);
    xK = aux.centroid(idNv,1); 
    yK = aux.centroid(idNv,2);
    hK = aux.diameter(idNv);    
    v1 = 1:Nv;  v2 = [2:Nv,1];
    rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
    Nrotx = y(:,rotid2) - y(:,rotid1);
    Nroty = x(:,rotid1)-x(:,rotid2); 
    
    % DoF values of scaled monomials
    xKv = repmat(xK,1,Nv);  yKv = repmat(yK,1,Nv);  
    hKv = repmat(hK,1,Nv);
    m = zeros(3*NTv, Nv);
    m(1:3:end,:) = 1 + 0*x;  
    m(2:3:end,:) = (x-xKv)./hKv; 
    m(3:3:end,:) = (y-yKv)./hKv;     
    
    % transition matrix
    Dc(idNv) = mat2cell(m', Nv, 3*ones(1,NTv));  
    
    % elliptic projection
    B = zeros(3*NTv, Nv);
    mx = [0*x, 1./hKv, 0*x];  
    my = [0*x, 0*x, 1./hKv];
    for i = 1:3
        col = (i-1)*Nv+1:i*Nv;
        B(i:3:end,:) = 0.5*(mx(:,col).*Nrotx + my(:,col).*Nroty);
    end
    Bc(idNv) = mat2cell(B, 3*ones(NTv,1), Nv);
    B(1:3:end,:) = 1/Nv;  % constraint
    Bsc(idNv) = mat2cell(B, 3*ones(NTv,1), Nv); 
    
    % L2 projection
    H = zeros(3*NTv, 3);
    for s = 1:Nv  % loop over triangles of each element
        elemTri = [elemNv(:,v1(s)), elemNv(:,v2(s)), N+idNv];
        areaTri = simplexvolume(nodeTri,elemTri);
        for p = 1:length(weight)
            pxy = lambda(p,1)*nodeTri(elemTri(:,1),:) ...
                + lambda(p,2)*nodeTri(elemTri(:,2),:) ...
                + lambda(p,3)*nodeTri(elemTri(:,3),:);
            x = pxy(:,1);  y = pxy(:,2);
            mp = [ones(NTv,1), (x-xK)./hK, (y-yK)./hK];
            for i = 1:3   
                H(i:3:end,:) = H(i:3:end,:) ...
                    + weight(p)*repmat(areaTri.*mp(:,i),1,3).*mp;  
            end
        end
    end
    Hc(idNv) = mat2cell(H, 3*ones(NTv,1), 3); 
end

%% Assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  
nnz = sum(elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
nnz = sum(elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT    
    % matrices
    index = elem{iel};  Nv = length(index);
    hK = aux.diameter(iel);
    D = Dc{iel};  
    B = Bc{iel};  Bs = Bsc{iel};
    G = B*D;      Gs = Bs*D;
    H = Hc{iel};

    % local stiffness matrix
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pis'*H*Pis + pde.c*hK^2*(I-Pi)'*(I-Pi);  % Pis = Pi0s;
    AB = reshape(AK'+BK',1,[]);

    % load vector
    fK = Pis'*[pde.f(aux.centroid(iel,:))*aux.area(iel);0;0];

    % assembly index for ellptic projection
    indexDof = index;  Ndof = Nv;  % local to global index
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;

    % assembly index for right hand side
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;

    % matrix for L2 and H1 error evaluation
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,N,N);
ff = accumarray(elemb,Fb,[N 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN;
if ~isempty(bdEdgeN)
    g_N = pde.Du;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    F1 = 0.5*sum(Ne.*g_N(z1),2);
    F2 = 0.5*sum(Ne.*g_N(z2),2);
    FN = [F1,F2];
    ff = ff + accumarray(bdEdgeN(:), FN(:),[N 1]);
end

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;  bdNodeIdxD = bdStruct.bdNodeIdxD;
isBdNode = false(N,1); isBdNode(bdNodeIdxD) = true;
bdDof = find(isBdNode); freeDof = find(~isBdNode);
nodeD = node(bdDof,:);
u = zeros(N,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
solver = 'amg';
if N < 2e3, solver = 'direct'; end
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
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm
info.D = Dc;