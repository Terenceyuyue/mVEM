function [u,info] = PoissonVEMk3(node,elem,pde,bdStruct)
%PoissonVEMk3 solves Poisson equation using virtual element method in V3/W3
%
%     -\Delta u = f,  in Omega
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
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+2*NE+3*NT;  % total dof number
Nm = 10; % number of scaled monomials

%% Compute projection matrices
D = cell(NT,1);  Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1);
H = cell(NT,1);  C = cell(NT,1);
Hf = cell(NT,1); Cf = cell(NT,1);

sgnelem = cell(NT,1); % elementwise signs of edges
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
for iel = 1:NT
    % ----------- element information ----------
    index = elem{iel};     Nv = length(index);     Ndof = 3*Nv+3;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); z = [x,y];
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    nodeT = [node(index,:);centroid(iel,:)]; % triangulation of K
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    za = z(v1,:)+r(1)*(z(v2,:)-z(v1,:)); % Gauss-Lobatto
    zb = z(v1,:)+r(2)*(z(v2,:)-z(v1,:));
    % ------------- scaled monomials -------------
    m1 = @(x,y) 1+0*x;
    m2 = @(x,y) (x-xK)/hK;
    m3 = @(x,y) (y-yK)/hK;
    m4 = @(x,y) (x-xK).^2/hK^2;
    m5 = @(x,y) (x-xK).*(y-yK)/hK^2;
    m6 = @(x,y) (y-yK).^2/hK^2;
    m7 = @(x,y) (x-xK).^3/hK^3;
    m8 = @(x,y) (x-xK).^2.*(y-yK)/hK^3;
    m9 = @(x,y) (x-xK).*(y-yK).^2/hK^3;
    m10 = @(x,y) (y-yK).^3/hK^3;
    gradm1 = @(x,y) [0+0*x, 0+0*x];
    gradm2 = @(x,y) [1/hK+0*x, 0+0*x];
    gradm3 = @(x,y) [0+0*x, 1/hK+0*x];
    gradm4 = @(x,y) [2*(x-xK)/hK^2, 0+0*x];
    gradm5 = @(x,y) [(y-yK)/hK^2,  (x-xK)/hK^2];
    gradm6 = @(x,y) [0+0*x, 2*(y-yK)/hK^2];
    gradm7 = @(x,y) [3*(x-xK).^2/hK^3, 0+0*x];
    gradm8 = @(x,y) [2*(x-xK).*(y-yK)/hK^3, (x-xK).^2/hK^3];
    gradm9 = @(x,y) [(y-yK).^2/hK^3, 2*(x-xK).*(y-yK)/hK^3];
    gradm10 = @(x,y) [0+0*x, 3*(y-yK).^2/hK^3];
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), ...
        m6(x,y), m7(x,y), m8(x,y), m9(x,y), m10(x,y)];
    Gradm = @(x,y) [gradm1(x,y); gradm2(x,y); gradm3(x,y); gradm4(x,y); gradm5(x,y); ...
        gradm6(x,y); gradm7(x,y); gradm8(x,y); gradm9(x,y); gradm10(x,y)];
    mm = {m1,m2,m3,m4,m5,m6,m7,m8,m9,m10};
    % ---------- transition matrix ------------
    D1 = zeros(Ndof,Nm);
    D1(1:Nv,:) = m(x,y); % at zi
    D1(Nv+1:2*Nv,:) = m(za(:,1), za(:,2)); % at ai
    D1(2*Nv+1:3*Nv,:) = m(zb(:,1), zb(:,2)); % at bi
    for j = 1:3
        mjv = @(x,y) repmat(mm{j}(x,y),1,10).*m(x,y); % moment j on K
        D1(3*Nv+j,:) = 1/area(iel)*integralTri(mjv,4,nodeT,elemT);
    end
    D{iel} = D1;
    % ----------------- elliptic projection ----------------
    % first term
    I1 = zeros(Nm,Ndof);
    I1([4,6],3*Nv+1) = 2*area(iel)/hK^2;
    I1(7,3*Nv+2) = 6*area(iel)/hK^2;
    I1(8,3*Nv+3) = 2*area(iel)/hK^2;
    I1(9,3*Nv+2) = 2*area(iel)/hK^2;
    I1(10,3*Nv+3) = 6*area(iel)/hK^2;
    % second term
    I2 = zeros(Nm,Ndof);
    for j = 1:Nv % loop of edges
        % he*\partial_n(m) at zj, aj, bj, zj+1
        Dn1 = Gradm(x(j),y(j))*Ne(j,:)';
        Dn2 = Gradm(za(j,1),za(j,2))*Ne(j,:)';
        Dn3 = Gradm(zb(j,1),zb(j,2))*Ne(j,:)';
        Dn4 = Gradm(x(v2(j)),y(v2(j)))*Ne(j,:)';
        % [ei,0,0]
        e1 = zeros(1,Ndof); e1(j) = 1;
        e2 = zeros(1,Ndof); e2(Nv+j) = 1;
        e3 = zeros(1,Ndof); e3(2*Nv+j) = 1;
        e4 = zeros(1,Ndof); e4(v2(j)) = 1;
        % f
        f1 = Dn1*e1; f2 = Dn2*e2; f3 = Dn3*e3; f4 = Dn4*e4;
        I2 = I2 + (w(1)*f1+w(2)*f2+w(3)*f3+w(4)*f4);
    end
    B1 = -I1 + 1/2*I2;
    % constraint
    B1s = B1;  B1s(1,3*Nv+1) = 1;  Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;
    % ---------- elliptic projection of lifting problem ----------
    DD = zeros(3*Nv+10, 10);
    DD(1:Ndof,:) = D1;
    for j = 4:10
        mjv = @(x,y) repmat(mm{j}(x,y),1,10).*m(x,y); % moment j on K
        DD(3*Nv+j,:) = 1/area(iel)*integralTri(mjv,6,nodeT,elemT);
    end
    BBs = zeros(10,3*Nv+10);
    BBs(:,1:Ndof) = B1s;
    GGs = BBs*DD;
    PPis = GGs\BBs;   PPi = DD*PPis;
    Dof = PPi(Ndof+1:3*Nv+10, 1:Ndof);
    % ------------- L2 projection ------------
    H1 = zeros(Nm,Nm);
    for i = 1:Nm
        for j = 1:Nm
            fun = @(x,y) mm{i}(x,y).*mm{j}(x,y);
            H1(i,j) = integralTri(fun,6,nodeT,elemT);
        end
    end
    H{iel} = H1;
    C1 = zeros(Nm,Ndof);
    C1(1,3*Nv+1) = area(iel);
    C1(2,3*Nv+2) = area(iel);
    C1(3,3*Nv+3) = area(iel);
    C1(4:end,:) = area(iel)*Dof;
    C{iel} = C1;
    % ---------- L2 projection for rhs ----------
    Hf{iel} = H1(1:3,1:3);
    Cf{iel} = C1(1:3,:);
    % --------- sgnelem -----------
    sgnL = sign(index(v2)-index(v1));
    id = elem2edge{iel}; sgnbd = E(id); sgnL(sgnbd) = 1;
    sgnL(sgnL==-1) = 0;
    sgnelem{iel} = sgnL;
end

%% Get elementwise stiffness matrix and load vector
ABelem = cell(NT,1); belem = cell(NT,1);
Ph = cell(NT,1); % matrix for error evaluation
for iel = 1:NT
    % element information
    index = elem{iel};   Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    % Projection
    Pis = Gs{iel}\Bs{iel}; Pi  = D{iel}*Pis; I = eye(size(Pi)); % elliptic
    Pi0s = H{iel}\C{iel};  Pi0 = D{iel}*Pi0s; % L2
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pi0s'*H{iel}*Pi0s + pde.c*hK^2*(I-Pi0)'*(I-Pi0);
    ABelem{iel} = reshape(AK'+BK',1,[]); % straighten as row vector for easy assembly
    % Load vector (L2 projection Pi_{k-2})
    m1 = @(x,y) 1+0*x;
    m2 = @(x,y) (x-xK)/hK;
    m3 = @(x,y) (y-yK)/hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    fm = @(x,y) repmat(pde.f([x,y]),1,3).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,4,nodeT,elemT); rhs = rhs';
    Pifs = Hf{iel}\Cf{iel};
    fK = Pifs'*rhs;    
    belem{iel} = fK'; % straighten as row vector for easy assembly
    % matrix for error evaluation
    Ph{iel} = Pis;
end
clear AK BK fK rhs Pis Pi I Pi0s Pi0 D Bs G Gs H C Hf Cf;

%% Assemble stiffness matrix and load vector
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnz = sum((3*elemLen+3).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
id = 0; ff = zeros(NNdof,1);
elem2dof = cell(NT,1);
for Nv = vertNum(:)' % only valid for row vector
    
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices
    
    elemNv = cell2mat(elem(idNv)); % elem
    elem2edgeNv = cell2mat(elem2edge(idNv));
    sgnelemNv = cell2mat(sgnelem(idNv));
    elemNva = elem2edgeNv+N*sgnelemNv+(N+NE)*(~sgnelemNv);
    elemNvb = elem2edgeNv+(N+NE)*sgnelemNv +N*(~sgnelemNv);
    elem2 = [elemNv, elemNva, elemNvb,...
             idNv+N+2*NE, idNv+N+2*NE+NT, idNv+N+2*NE+2*NT];
    K = cell2mat(ABelem(idNv)); F = cell2mat(belem(idNv));
    Ndof = 3*Nv+3;
    ii(id+1:id+NTv*Ndof^2) = reshape(repmat(elem2, Ndof,1), [], 1);
    jj(id+1:id+NTv*Ndof^2) = repmat(elem2(:), Ndof, 1);
    ss(id+1:id+NTv*Ndof^2) = K(:);
    id = id + NTv*Ndof^2;
    
    % assemble the vector
    ff = ff +  accumarray(elem2(:),F(:),[NNdof 1]);
    
    % elementwise global indices
    elem2dof(idNv) = mat2cell(elem2, ones(NTv,1), Ndof);
end
kk = sparse(ii,jj,ss,NNdof,NNdof);

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
    Fb = w(3)*sum(Ne.*g_N(zb),2)/2;
    FN = [F1,F2,Fa,Fb];
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
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; info.freeDof = freeDof;
