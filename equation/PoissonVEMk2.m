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
NNdof = N+NE+NT;  % total dof number

%% Compute projection matrices
D = cell(NT,1);  Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1); H = cell(NT,1);
for iel = 1:NT
    % ------- element information ----------
    index = elem{iel};     Nv = length(index);     Ndof = 2*Nv+1;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); 
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    % scaled monomials
    m1 = @(x,y) 1+0*x; 
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;
    m4 = @(x,y) (x-xK).^2/hK^2; 
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2; 
    m6 = @(x,y) (y-yK).^2./hK^2;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)]; 
    Gradm = @(x,y) [[0,0]; [1, 0]/hK; [0, 1]/hK; [2*(x-xK), 0]/hK^2;
                    [(y-yK), (x-xK)]/hK^2; [0, 2*(y-yK)]/hK^2]; % grad m
    % -------- transition matrix ----------
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    D1 = zeros(Ndof,6); 
    D1(1:2*Nv,:) = [m(x,y); m(xe,ye)];  % vertices and mid-edge points
    D1(end,:) = 1/area(iel)*integralTri(m,3,nodeT,elemT);
    D{iel} = D1;
    % --------- elliptic projection -----------
    % first term
    Lapm = zeros(6,1); Lapm([4,6]) = 2/hK^2;
    Dof = [zeros(1,2*Nv), area(iel)];
    I1 = Lapm*Dof;
    % second term
    I2 = 0;
    for j = 1:Nv % loop of edges
        % he*\partial_n(m) at z_j, ze, z_j+1
        DnL = Gradm(x(j),y(j))*Ne(j,:)';
        Dne = Gradm(xe(j),ye(j))*Ne(j,:)';
        DnR = Gradm(x(v2(j)),y(v2(j)))*Ne(j,:)';                
        % [ei,0,0]
        e1 = zeros(1,Ndof); e1(j) = 1;
        e2 = zeros(1,Ndof); e2(Nv+j) = 1;
        e3 = zeros(1,Ndof); e3(v2(j)) = 1;
        % f
        f1 = DnL*e1; f2 = Dne*e2; f3 = DnR*e3;
        I2 = I2 + (f1+4*f2+f3);
    end
    B1 = -I1 + 1/6*I2;
    % constraint
    B1s = B1;  B1s(1,:) = Dof;  Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;
    % --------- L2 projection -----------
    m = {m1,m2,m3,m4,m5,m6};
    H1 = zeros(6,6);
    for i = 1:6
        for j = 1:6
            fun = @(x,y) m{i}(x,y).*m{j}(x,y);
            H1(i,j) = integralTri(fun,3,nodeT,elemT);
        end
    end
    H{iel} = H1;
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
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));    
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pis'*H{iel}*Pis + pde.c*hK^2*(I-Pi)'*(I-Pi);  % Pis = Pi0s;
    ABelem{iel} = reshape(AK'+BK',1,[]); % straighten as row vector for easy assembly    
    % Load vector   
    m1 = @(x,y) 1+0*x; 
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;
    m4 = @(x,y) (x-xK).^2/hK^2; 
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2; 
    m6 = @(x,y) (y-yK).^2./hK^2;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)];
    fm = @(x,y) repmat(pde.f([x,y]),1,6).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,3,nodeT,elemT); rhs = rhs';
    fK = Pis'*rhs;
    belem{iel} = fK'; % straighten as row vector for easy assembly
    % matrix for error evaluation
    Ph{iel} = Pis; 
end
clear AK BK fK Pis Pi I D Bs G Gs H;

%% Assemble stiffness matrix and load vector
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnz = sum((2*elemLen+1).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
id = 0; ff = zeros(NNdof,1); 
elem2dof = cell(NT,1);
for Nv = vertNum(:)' % only valid for row vector
    
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices
    
    elemNv = cell2mat(elem(idNv)); % elem     
    elem2edgeNv = cell2mat(elem2edge(idNv));    
    elem2 = [elemNv,elem2edgeNv+N,idNv+N+NE];
    K = cell2mat(ABelem(idNv)); F = cell2mat(belem(idNv));
    Ndof = 2*Nv+1;
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
bdNodeIdx = bdStruct.bdNodeIdx; 
bdEdgeD = bdStruct.bdEdgeD;
bdEdgeIdxD = bdStruct.bdEdgeIdxD; 
id = [bdNodeIdx; bdEdgeIdxD+N];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
pD = node(bdNodeIdx,:); uD = g_D(pD);
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
info.kk = kk; info.freeDof = freeDof;
