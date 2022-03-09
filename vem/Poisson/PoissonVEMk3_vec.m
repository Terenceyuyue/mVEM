function [u,info] = PoissonVEMk3_vec(node,elem,pde,bdStruct)
%PoissonVEMk3_vec solves Poisson equation using virtual element method in W3
%with a vectorized implementation.
%
%     -\Delta u + cu = f,  in Omega
%     Dirichlet boundary condition u=g_D on \Gamma_D,
%     Neumann boundary condition   grad(u)*n=g_N on \Gamma_N
%
% Copyright (C)  Terence Yu.

%% Gauss-Lobatto weights and points
r = [-1, -1/sqrt(5), 1/sqrt(5), 1]; % [-1,1]
r = (r(2:3)+1)/2; % [0,1]: gives the ratios of interior nodes
w = [1/6, 5/6, 5/6, 1/6];

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+2*NE+3*NT;  % total dof number
Nm = 10; % number of scaled monomials
elemLen = cellfun('length',elem);
vertNum = unique(elemLen);

%% Triangulation
nodeTri = [node; aux.centroid];
% elemTri: [a_i, a_{i+1}, centroid]

%% Quadrature points
[lambda,weight] = quadpts(4);

%% Compute projection matrices
Dc = cell(NT,1); 
Bc = cell(NT,1);  Bsc = cell(NT,1);
Hc = cell(NT,1);  fm = zeros(NT,3);
for Nv = vertNum(:)'
    % ---- information of polygons with Nv vertices ----
    idNv = find(elemLen == Nv);  NTv = length(idNv);  
    elemNv = cell2mat(elem(idNv)); % elem
    xK = aux.centroid(idNv,1); yK = aux.centroid(idNv,2);
    hK = aux.diameter(idNv);
    x = reshape(node(elemNv,1),NTv,Nv);
    y = reshape(node(elemNv,2),NTv,Nv);
    v1 = 1:Nv;  v2 = [2:Nv,1];
    xa = x(:,v1)+r(1)*(x(:,v2)-x(:,v1)); % Gauss-Lobatto
    ya = y(:,v1)+r(1)*(y(:,v2)-y(:,v1));
    xb = x(:,v1)+r(2)*(x(:,v2)-x(:,v1));
    yb = y(:,v1)+r(2)*(y(:,v2)-y(:,v1));
    Nex = y(:,v2)-y(:,v1);  Ney = x(:,v1)-x(:,v2);
    
    % ---- DoF values of scaled monomials ----
    % the moments are given by L2 projection
    Ndof = 3*Nv+3;
    xKv = repmat(xK,1,Nv);  yKv = repmat(yK,1,Nv);  
    hKv = repmat(hK,1,Nv);
    m = zeros(Nm*NTv, Ndof);
    m(1:Nm:end,1:3*Nv) = [1+0*x, 1+0*xa, 1+0*xb];  
    m(2:Nm:end,1:3*Nv) = [x-xKv, xa-xKv, xb-xKv]./repmat(hKv,1,3); 
    m(3:Nm:end,1:3*Nv) = [y-yKv, ya-yKv, yb-yKv]./repmat(hKv,1,3); 
    m(4:Nm:end,1:3*Nv) = [(x-xKv).^2, (xa-xKv).^2, (xb-xKv).^2]./repmat(hKv.^2,1,3);  
    m(5:Nm:end,1:3*Nv) = [(x-xKv).*(y-yKv), (xa-xK).*(ya-yKv), (xb-xK).*(yb-yKv)]./repmat(hKv.^2,1,3); 
    m(6:Nm:end,1:3*Nv) = [(y-yKv).^2, (ya-yKv).^2, (yb-yKv).^2]./repmat(hKv.^2,1,3); 
    m(7:Nm:end,1:3*Nv) = [(x-xKv).^3, (xa-xKv).^3, (xb-xKv).^3]./repmat(hKv.^3,1,3); 
    m(8:Nm:end,1:3*Nv) = [(x-xKv).^2.*(y-yKv), (xa-xKv).^2.*(ya-yKv), (xb-xKv).^2.*(yb-yKv)]./repmat(hKv.^3,1,3); 
    m(9:Nm:end,1:3*Nv) = [(x-xKv).*(y-yKv).^2, (xa-xKv).*(ya-yKv).^2, (xb-xKv).*(yb-yKv).^2]./repmat(hKv.^3,1,3); 
    m(Nm:Nm:end,1:3*Nv) = [(y-yKv).^3, (ya-yKv).^3, (yb-yKv).^3]./repmat(hKv.^3,1,3);
    
    % ---- L2 projection ----
    H = zeros(Nm*NTv, Nm);  rhs = zeros(NTv,3);
    for s = 1:Nv  % loop over triangles of each element
        elemTri = [elemNv(:,v1(s)), elemNv(:,v2(s)), N+idNv];
        areaTri = simplexvolume(nodeTri,elemTri);
        for p = 1:length(weight)
            pxy = lambda(p,1)*nodeTri(elemTri(:,1),:) ...
                + lambda(p,2)*nodeTri(elemTri(:,2),:) ...
                + lambda(p,3)*nodeTri(elemTri(:,3),:);
            xp = pxy(:,1);  yp = pxy(:,2);
            mp = [1+0*xp, (xp-xK)./hK, (yp-yK)./hK, ...
                [(xp-xK).^2, (xp-xK).*(yp-yK), (yp-yK).^2]./repmat(hK.^2,1,3), ...
                [(xp-xK).^3, (xp-xK).^2.*(yp-yK), (xp-xK).*(yp-yK).^2, (yp-yK).^3]./repmat(hK.^3,1,4)
                ];
            for i = 1:Nm                   
                H(i:Nm:end,:) = H(i:Nm:end,:) ...
                    + weight(p)*repmat(areaTri.*mp(:,i),1,Nm).*mp;
            end
            rhs = rhs + ...
                weight(p)*repmat(areaTri.*pde.f([xp,yp]),1,3).*mp(:,1:3);
        end
    end
    Hc(idNv) = mat2cell(H, Nm*ones(NTv,1), Nm); 
    m(:,end-2:end) = H(:,1:3)./repmat(kron(aux.area(idNv),ones(Nm,1)),1,3); % moments of m
    fm(idNv,:) = rhs;
    
    % ---- transition matrix ----  
    Dc(idNv) = mat2cell(m', Ndof, Nm*ones(1,NTv));
    
    % ---- elliptic projection ----
    B = zeros(Nm*NTv, Ndof);
    % first term
    temp = 1./hK.^2.*aux.area(idNv);
    B(4:Nm:end,3*Nv+1) = -2*temp; 
    B(6:Nm:end,3*Nv+1) = -2*temp;
    B(7:Nm:end,3*Nv+2) = -6*temp;
    B(8:Nm:end,3*Nv+3) = -2*temp;
    B(9:Nm:end,3*Nv+2) = -2*temp;
    B(Nm:Nm:end,3*Nv+3) = -6*temp;
    % second term
    mx = [0*x, 1./hKv+0*x, 0*x, 2*(x-xKv)./hKv.^2, (y-yKv)./hKv.^2, 0*x, ...
          3*(x-xKv).^2./hKv.^3, 2*(x-xKv).*(y-yKv)./hKv.^3, (y-yKv).^2./hK.^3, 0*x];  
    my = [0*x, 0*x, 1./hKv+0*x, 0*x, (x-xKv)./hKv.^2, 2*(y-yKv)./hKv.^2, ...
        0*x, (x-xKv).^2./hKv.^3, 2*(x-xKv).*(y-yKv)./hK.^3, 3*(y-yKv).^2./hKv.^3]; 
    mxa = [0*x, 1./hKv+0*x, 0*x, 2*(xa-xKv)./hKv.^2, (ya-yKv)./hKv.^2, 0*x, ...
          3*(xa-xKv).^2./hKv.^3, 2*(xa-xKv).*(ya-yKv)./hKv.^3, (ya-yKv).^2./hK.^3, 0*x];  
    mya = [0*x, 0*x, 1./hKv+0*x, 0*x, (xa-xKv)./hKv.^2, 2*(ya-yKv)./hKv.^2, ...
        0*x, (xa-xKv).^2./hKv.^3, 2*(xa-xKv).*(ya-yKv)./hK.^3, 3*(ya-yKv).^2./hKv.^3];
    mxb = [0*x, 1./hKv+0*x, 0*x, 2*(xb-xKv)./hKv.^2, (yb-yKv)./hKv.^2, 0*x, ...
          3*(xb-xKv).^2./hKv.^3, 2*(xb-xKv).*(yb-yKv)./hKv.^3, (yb-yKv).^2./hK.^3, 0*x];  
    myb = [0*x, 0*x, 1./hKv+0*x, 0*x, (xb-xKv)./hKv.^2, 2*(yb-yKv)./hKv.^2, ...
        0*x, (xb-xKv).^2./hKv.^3, 2*(xb-xKv).*(yb-yKv)./hK.^3, 3*(yb-yKv).^2./hKv.^3];
    for i = 1:Nm
        col = (i-1)*Nv+1:i*Nv;
        B(i:Nm:end,1:Nv) = B(i:Nm:end,1:Nv) ...
            + w(1)/2*( mx(:,col).*Nex + my(:,col).*Ney );  
        B(i:Nm:end,Nv+1:2*Nv) = B(i:Nm:end,Nv+1:2*Nv) ...
            + w(2)/2*( mxa(:,col).*Nex + mya(:,col).*Ney ); 
        B(i:Nm:end,2*Nv+1:3*Nv) = B(i:Nm:end,2*Nv+1:3*Nv) ...
            + w(3)/2*( mxb(:,col).*Nex + myb(:,col).*Ney ); 
        B(i:Nm:end,v2) = B(i:Nm:end,v2) ...
            + w(4)/2*( mx(:,col(v2)).*Nex + my(:,col(v2)).*Ney ); 
    end
    Bc(idNv) = mat2cell(B, Nm*ones(NTv,1), Ndof);
    % constraint
    B(1:Nm:end,3*Nv+1) = aux.area(idNv);
    Bsc(idNv) = mat2cell(B, Nm*ones(NTv,1), Ndof); 
end

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum((3*elemLen+3).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(3*elemLen+3);
elemf = zeros(nnz,1); Ff = zeros(nnz,1);
ia = 0; ib = 0;
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
for iel = 1:NT
    % information
    index = elem{iel};   indexEdge = elem2edge{iel};
    Nv = length(index);  Ndof = 3*Nv+3;
    hK = diameter(iel);
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    
    % elliptic projection
    D = Dc{iel};
    B = Bc{iel}; Bs = Bsc{iel};
    G = B*D;     Gs = Bs*D;
    
    % L2 projection
    H = Hc{iel};        
    DD = zeros(3*Nv+Nm, Nm); % elliptic projection of lifting problem
    DD(1:Ndof,:) = D;
    DD(3*Nv+4:end,:) = 1/area(iel)*H(4:end,:);
    BBs = zeros(Nm,3*Nv+Nm);
    BBs(:,1:Ndof) = Bs;
    GGs = BBs*DD;
    PPis = GGs\BBs;   PPi = DD*PPis;
    Dof = PPi(Ndof+1:3*Nv+Nm, 1:Ndof);    
    C = zeros(Nm,Ndof); %
    C(1,3*Nv+1) = area(iel);
    C(2,3*Nv+2) = area(iel);
    C(3,3*Nv+3) = area(iel);
    C(4:end,:) = area(iel)*Dof;    
    
    % L2 projection for rhs
    Hf = H(1:3,1:3);
    Cf = C(1:3,:);
    
    % local stiffness matrix 
    Pis = Gs\Bs; Pi  = D*Pis; I = eye(size(Pi)); % elliptic
    Pi0s = H\C;  Pi0 = D*Pi0s; % L2
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pi0s'*H*Pi0s + pde.c*hK^2*(I-Pi0)'*(I-Pi0);
    AB = reshape(AK'+BK',1,[]); 
    
    % load vector
    rhs = fm(iel,:);
    Pifs = Hf\Cf;
    fK = Pifs'*rhs(:);
    
    % local to global index
    % sgnBase
    sgnBase = sign(index(v2)-index(v1));
    id = elem2edge{iel}; sgnbd = E(id); sgnBase(sgnbd) = 1;
    sgnBase(sgnBase==-1) = 0;
    elema = indexEdge+N*sgnBase+(N+NE)*(~sgnBase);
    elemb = indexEdge+(N+NE)*sgnBase +N*(~sgnBase);
    indexDof = [index, elema, elemb,...
        iel+N+2*NE, iel+N+2*NE+NT, iel+N+2*NE+2*NT];
    
    % assembly index for ellptic projection
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
    
    % assembly index for right hand side 
    elemf(ib+1:ib+Ndof) = indexDof(:);
    Ff(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % matrix for L2 and H1 error evaluation
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,NNdof,NNdof);
ff = accumarray(elemf,Ff,[NNdof 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN; bdEdgeIdxN = bdStruct.bdEdgeIdxN;
if ~isempty(bdEdgeN)
    g_N = pde.Du;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    za = z1+r(1)*(z2-z1); zb = z1+r(2)*(z2-z1); % Gauss-Lobatto
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    F1 = w(1)*sum(Ne.*g_N(z1),2)/2;
    F2 = w(4)*sum(Ne.*g_N(z2),2)/2;
    Fa = w(2)*sum(Ne.*g_N(za),2)/2;
    Ff = w(3)*sum(Ne.*g_N(zb),2)/2;
    FN = [F1,F2,Fa,Ff];
    bdEdgeN2 = [bdEdgeN, bdEdgeIdxN+N,  bdEdgeIdxN+N+NE];
    ff = ff + accumarray(bdEdgeN2(:), FN(:),[NNdof 1]);
end

%% Apply Dirichlet boundary conditions
% boundary information
g_D = pde.g_D;
bdNodeIdx = bdStruct.bdNodeIdx; bdEdgeD = bdStruct.bdEdgeD;
bdEdgeIdxD = bdStruct.bdEdgeIdxD;
bdEdgeIdxa = bdEdgeIdxD + N;
bdEdgeIdxb = bdEdgeIdxD + N+NE;
id = [bdNodeIdx; bdEdgeIdxa; bdEdgeIdxb];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
pD = node(bdNodeIdx,:); uD = g_D(pD);
% ai, bi on the boundary
z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:);
za = z1+r(1)*(z2-z1); zb = z1+r(2)*(z2-z1); % Gauss-Lobatto
uDa = g_D(za); uDb = g_D(zb);
% rhs
u = zeros(NNdof,1); u(bdDof) = [uD; uDa; uDb];
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
% For k = 3 we need to remove all boundary d.o.f.s
bdEdgeIdx = bdStruct.bdEdgeIdx; 
bdEdgeIdxa = bdEdgeIdx + N;
bdEdgeIdxb = bdEdgeIdx + N+NE;
id = [bdNodeIdx; bdEdgeIdxa; bdEdgeIdxb];
isBdNode = false(NNdof,1); isBdNode(id) = true;
DofI = (~isBdNode);  % interior d.o.f.s
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; info.DofI = DofI;  % d.o.f.s for computing errors in the energy norm