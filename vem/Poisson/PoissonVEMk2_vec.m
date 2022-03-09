function [u,info] = PoissonVEMk2_vec(node,elem,pde,bdStruct)
%PoissonVEMk2_vec solves Poisson equation using virtual element method in V2
%with a vectorized implementation.
%
%     -\Delta u + cu = f,  in Omega
%     Dirichlet boundary condition u=g_D on \Gamma_D, 
%     Neumann boundary condition   grad(u)*n=g_N on \Gamma_N
%
% Copyright (C)  Terence Yu. 

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+NE+NT;  % total dof number
Nm = 6;
elemLen = cellfun('length',elem);
vertNum = unique(elemLen);

%% Triangulation
nodeTri = [node; aux.centroid];
% elemTri: [a_i, a_{i+1}, centroid]

%% Quadrature points
[lambda,weight] = quadpts(3);

%% Compute projection matrices
Dc = cell(NT,1); 
Bc = cell(NT,1);  Bsc = cell(NT,1);
Hc = cell(NT,1);  fm = zeros(NT,Nm);
for Nv = vertNum(:)'
    % ---- information of polygons with Nv vertices ----
    idNv = find(elemLen == Nv);  NTv = length(idNv);  
    elemNv = cell2mat(elem(idNv)); % elem
    xK = aux.centroid(idNv,1); yK = aux.centroid(idNv,2);
    hK = aux.diameter(idNv);
    x = reshape(node(elemNv,1),NTv,Nv);
    y = reshape(node(elemNv,2),NTv,Nv);
    v1 = 1:Nv;  v2 = [2:Nv,1];
    xe = (x(:,v1)+x(:,v2))/2;  % mid-edge points
    ye = (y(:,v1)+y(:,v2))/2; 
    Nex = y(:,v2)-y(:,v1);  Ney = x(:,v1)-x(:,v2);
    
    % ---- DoF values of scaled monomials ----
    % the moments are given by L2 projection
    Ndof = 2*Nv+1;
    xKv = repmat(xK,1,Nv);  yKv = repmat(yK,1,Nv);  
    hKv = repmat(hK,1,Nv);
    m = zeros(Nm*NTv, Ndof);
    m(1:Nm:end,1:2*Nv) = [1+0*x, 1+0*xe];  
    m(2:Nm:end,1:2*Nv) = [(x-xKv)./hKv, (xe-xKv)./hKv]; 
    m(3:Nm:end,1:2*Nv) = [(y-yKv)./hKv, (ye-yKv)./hKv]; 
    m(4:Nm:end,1:2*Nv) = [(x-xKv).^2./hKv.^2, (xe-xKv).^2./hKv.^2];  
    m(5:Nm:end,1:2*Nv) = [(x-xKv).*(y-yKv)./hKv.^2, (xe-xK).*(ye-yKv)./hKv.^2]; 
    m(Nm:Nm:end,1:2*Nv) = [(y-yKv).^2./hKv.^2, (ye-yKv).^2./hKv.^2]; 
    
    % ---- L2 projection ----
    H = zeros(Nm*NTv, Nm);  rhs = zeros(NTv,Nm);
    for s = 1:Nv  % loop over triangles of each element
        elemTri = [elemNv(:,v1(s)), elemNv(:,v2(s)), N+idNv];
        areaTri = simplexvolume(nodeTri,elemTri);
        for p = 1:length(weight)
            pxy = lambda(p,1)*nodeTri(elemTri(:,1),:) ...
                + lambda(p,2)*nodeTri(elemTri(:,2),:) ...
                + lambda(p,3)*nodeTri(elemTri(:,3),:);
            xp = pxy(:,1);  yp = pxy(:,2);
            mp = [1+0*xp, (xp-xK)./hK, (yp-yK)./hK, ...
                (xp-xK).^2./hK.^2, (xp-xK).*(yp-yK)./hK.^2, (yp-yK).^2./hK.^2];
            for i = 1:Nm                   
                H(i:Nm:end,:) = H(i:Nm:end,:) ...
                    + weight(p)*repmat(areaTri.*mp(:,i),1,Nm).*mp;
            end
            rhs = rhs + ...
                weight(p)*repmat(areaTri.*pde.f([xp,yp]),1,Nm).*mp;
        end
    end
    Hc(idNv) = mat2cell(H, Nm*ones(NTv,1), Nm); 
    m(:,end) = H(:,1)./kron(aux.area(idNv),ones(Nm,1)); % moments of m
    fm(idNv,:) = rhs;
    
    % ---- transition matrix ----  
    Dc(idNv) = mat2cell(m', Ndof, Nm*ones(1,NTv));  
    
    % ---- elliptic projection ----
    B = zeros(Nm*NTv, Ndof);
    % first term
    B(4:Nm:end,end) = -2./hK.^2.*aux.area(idNv); 
    B(Nm:Nm:end,end) = B(4:Nm:end,end);
    % second term
    mx = [0*x, 1./hKv+0*x, 0*x, 2*(x-xKv)./hKv.^2, (y-yKv)./hKv.^2, 0*x];  
    my = [0*x, 0*x, 1./hKv+0*x, 0*x, (x-xKv)./hKv.^2, 2*(y-yKv)./hKv.^2];  
    mxe = [0*x, 1./hKv+0*x, 0*x, 2*(xe-xKv)./hKv.^2, (ye-yKv)./hKv.^2, 0*x];  
    mye = [0*x, 0*x, 1./hKv+0*x, 0*x, (xe-xKv)./hKv.^2, 2*(ye-yKv)./hKv.^2]; 
    for i = 1:Nm
        col = (i-1)*Nv+1:i*Nv;
        B(i:Nm:end,v1) = B(i:Nm:end,v1) ...
            + 1/Nm*( mx(:,col).*Nex + my(:,col).*Ney ); % zi
        B(i:Nm:end,Nv+1:2*Nv) = B(i:Nm:end,Nv+1:2*Nv) ...
            + 4/Nm*( mxe(:,col).*Nex + mye(:,col).*Ney ); % z_{i+1/2}
        B(i:Nm:end,v2) = B(i:Nm:end,v2) ...
            + 1/Nm*( mx(:,col(v2)).*Nex + my(:,col(v2)).*Ney ); % z_{i+1}
    end
    Bc(idNv) = mat2cell(B, Nm*ones(NTv,1), Ndof);
    % constraint
    B(1:Nm:end,end) = aux.area(idNv);
    Bsc(idNv) = mat2cell(B, Nm*ones(NTv,1), Ndof); 
end

%% Assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen+1).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(2*elemLen+1);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % matrices
    index = elem{iel};  indexEdge = elem2edge{iel};   
    Nv = length(index); Ndof = 2*Nv+1;
    hK = aux.diameter(iel);
    D = Dc{iel};
    B = Bc{iel}; Bs = Bsc{iel};
    G = B*D;     Gs = Bs*D;
    H = Hc{iel};        
    
    % local stiffness matrix
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi)); 
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pis'*H*Pis + pde.c*hK^2*(I-Pi)'*(I-Pi);  % Pis = Pi0s;
    AB = reshape(AK'+BK',1,[]); 
    
    % load vector
    rhs = fm(iel,:);
    fK = Pis'*rhs(:);     

    % assembly index for ellptic projection
    indexDof = [index,indexEdge+N,iel+N+NE];  % local to global index
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
kk = sparse(ii,jj,ss,NNdof,NNdof);
ff = accumarray(elemb,Fb,[NNdof 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN; bdEdgeIdxN = bdStruct.bdEdgeIdxN;
if ~isempty(bdEdgeN)
    g_N = pde.Du;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:); zc = (z1+z2)/2;
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    F1 = sum(Ne.*g_N(z1),2)/Nm;
    F2 = sum(Ne.*g_N(z2),2)/Nm;
    F3 = sum(Ne.*g_N(zc),2)*4/Nm;
    FN = [F1,F2,F3];
    bdEdgeN2 = [bdEdgeN, bdEdgeIdxN+N];
    ff = ff + accumarray(bdEdgeN2(:), FN(:),[NNdof 1]);
end

%% Apply Dirichlet boundary conditions
% boundary information
g_D = pde.g_D;  
bdNodeIdxD = bdStruct.bdNodeIdxD; 
bdEdgeD = bdStruct.bdEdgeD;
bdEdgeIdxD = bdStruct.bdEdgeIdxD; 
id = [bdNodeIdxD; bdEdgeIdxD+N];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
pD = node(bdNodeIdxD,:); uD = g_D(pD);
% mid-edge on the boundary
z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:); zc = (z1+z2)./2;
uDc = g_D(zc);
% rhs
u = zeros(NNdof,1); u(bdDof) = [uD; uDc];
ff = ff - kk*u;

%% Set solver
solver = 'amg';
if NNdof < 2e3, solver = 'direct'; end
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