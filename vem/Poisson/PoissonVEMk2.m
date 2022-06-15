function [u,info] = PoissonVEMk2(node,elem,pde,bdStruct)
%PoissonVEMk2 solves Poisson equation using virtual element method in V2
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
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
Nm = 6;
NNdof = N+NE+NT;  % total dof number

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen+1).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(2*elemLen+1);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % ------- element information ----------
    index = elem{iel};  indexEdge = elem2edge{iel};   
    Nv = length(index); Ndof = 2*Nv+1;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); 
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    
    % ------------ scaled monomials --------
    m1 = @(x,y) 1+0*x;                    gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;               gradm2 = @(x,y) [1+0*x, 0+0*x]./hK;
    m3 = @(x,y) (y-yK)./hK;               gradm3 = @(x,y) [0+0*x, 1+0*x]./hK;
    m4 = @(x,y) (x-xK).^2/hK^2;           gradm4 = @(x,y) [2*(x-xK), 0+0*x]./hK^2;
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2;     gradm5 = @(x,y) [(y-yK), (x-xK)]./hK^2;
    m6 = @(x,y) (y-yK).^2./hK^2;          gradm6 = @(x,y) [0+0*x, 2*(y-yK)]./hK^2;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)]; 
    mc = {m1,m2,m3,m4,m5,m6};
    Gradmc = {gradm1, gradm2, gradm3, gradm4, gradm5, gradm6};
    
    % -------- transition matrix ----------
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    D = zeros(Ndof,Nm); 
    D(1:2*Nv,:) = [m(x,y); m(xe,ye)];  % vertices and mid-edge points
    D(end,:) = 1/area(iel)*integralTri(m,3,nodeT,elemT);
    
    % --------- elliptic projection -----------
    % first term
    Lapm = zeros(6,1); Lapm([4,6]) = 2/hK^2;
    Dof = [zeros(1,2*Nv), area(iel)];
    I1 = Lapm*Dof;
    % second term
    I2 = zeros(Nm, Ndof);
    elem1 = [v1(:), v2(:), v1(:)+Nv];
    for im = 1:Nm
        gradmc = Gradmc{im};
        F1 = 1/6*sum(gradmc(x(v1), y(v1)).*Ne, 2);
        F2 = 1/6*sum(gradmc(x(v2), y(v2)).*Ne, 2);
        F3 = 4/6*sum(gradmc(xe, ye).*Ne, 2);
        F = [F1, F2, F3];
        I2(im, :) = accumarray(elem1(:), F(:), [Ndof 1]);
    end
    B = -I1 + I2;
    % constraint
    Bs = B;  Bs(1,:) = Dof;  
    % consistency relation
    G = B*D;     Gs = Bs*D;
    
    % --------- L2 projection -----------
    H = zeros(Nm,Nm);
    for i = 1:Nm
        fun = @(x,y) repmat(mc{i}(x,y),1,Nm).*m(x,y);
        H(i,:) = integralTri(fun,3,nodeT,elemT);
    end
    
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi)); 
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pis'*H*Pis + pde.c*hK^2*(I-Pi)'*(I-Pi);  % Pis = Pi0s;
    AB = reshape(AK'+BK',1,[]); 
    
    % --------- load vector -----------
    fm = @(x,y) repmat(pde.f([x,y]),1,6).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,3,nodeT,elemT); rhs = rhs';
    fK = Pis'*rhs;     

    % --------- assembly index for ellptic projection -----------
    indexDof = [index,indexEdge+N,iel+N+NE];  % local to global index
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
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
    F1 = sum(Ne.*g_N(z1),2)/6;
    F2 = sum(Ne.*g_N(z2),2)/6;
    F3 = sum(Ne.*g_N(zc),2)*4/6;
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