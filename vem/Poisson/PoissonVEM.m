function [u,info] = PoissonVEM(node,elem,pde,bdStruct)
%PoissonVEM solves Poisson equation using virtual element method in V1
%
%     -\Delta u + cu = f,  in Omega
%     Dirichlet boundary condition u=g_D on \Gamma_D, 
%     Neumann boundary condition   grad(u)*n=g_N on \Gamma_N
%
% Copyright (C)  Terence Yu. 

%% Get auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid;  diameter = aux.diameter;  area = aux.area;
N = size(node,1);  NT = size(elem,1);
Nm = 3;

%% Compute projection matrices
D = cell(NT,1);  Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1); H = cell(NT,1);
for iel = 1:NT
    % ------- element information --------
    index = elem{iel};  Nv = length(index);    
    xK = centroid(iel,1); yK = centroid(iel,2); 
    hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1]; 
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % ------- scaled monomials ----------
    m1 = @(x,y) 1+0*x;                gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;           gradm2 = @(x,y) [1+0*x, 0+0*x]./hK;
    m3 = @(x,y) (y-yK)./hK;           gradm3 = @(x,y) [0+0*x, 1+0*x]./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];   
    mc = {m1,m2,m3};  % cell 
    Gradmc = {gradm1, gradm2, gradm3};    
    
    % -------- transition matrix ----------
    D1 = m(x,y);   D{iel} = D1;   
    
    % --------- elliptic projection -----------
    % first term  = 0
    B1 = zeros(Nm, Nv);
    % second term   
    elem1 = [v1(:), v2(:)];
    for im = 1:Nm
        gradmc = Gradmc{im};
        F1 = 0.5*sum(gradmc(x(v1), y(v1)).*Ne, 2);
        F2 = 0.5*sum(gradmc(x(v2), y(v2)).*Ne, 2);
        F = [F1, F2];
        B1(im, :) = accumarray(elem1(:), F(:), [Nv 1]);
    end  
    % constraint
    B1s = B1;  B1s(1,:) = 1/Nv;  Bs{iel} = B1s;
    % consistency relation
    G{iel} = B1*D1;  Gs{iel} = B1s*D1;  
    
    % --------- L2 projection -----------   
    H1 = zeros(Nm,Nm);
    for i = 1:Nm
        for j = 1:Nm
            fun = @(x,y) mc{i}(x,y).*mc{j}(x,y);
            H1(i,j) = integralTri(fun,3,nodeT,elemT);
        end
    end
    H{iel} = H1;
end

%% Get elementwise stiffness matrix and load vector
ABelem = cell(NT,1); belem = cell(NT,1);
Ph = cell(NT,1); % matrix for error evaluation
for iel = 1:NT
    % elementwise information
    hK = diameter(iel);
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));    
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pis'*H{iel}*Pis + pde.c*hK^2*(I-Pi)'*(I-Pi);  % Pis = Pi0s;
    ABelem{iel} = reshape(AK'+BK',1,[]); % straighten as row vector for easy assembly    
    % Load vector   
    fK = Pis'*[pde.f(centroid(iel,:))*area(iel);0;0];    
    belem{iel} = fK'; % straighten as row vector for easy assembly
    % matrix for error evaluation
    Ph{iel} = Pis; 
end
clear AK BK fK Pis Pi I Bs G Gs H;

%% Assemble stiffness matrix and load vector
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnz = sum(elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
id = 0; ff = zeros(N,1); 
elem2dof = cell(NT,1);
for Nv = vertNum(:)' % only valid for row vector
    
    % assemble the matrix
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices
    
    elemNv = cell2mat(elem(idNv)); % elem     
    K = cell2mat(ABelem(idNv)); F = cell2mat(belem(idNv));
    Ndof = Nv;
    ii(id+1:id+NTv*Ndof^2) = reshape(repmat(elemNv, Ndof,1), [], 1);
    jj(id+1:id+NTv*Ndof^2) = repmat(elemNv(:), Ndof, 1);
    ss(id+1:id+NTv*Ndof^2) = K(:);
    id = id + NTv*Ndof^2;
    
    % assemble the vector
    ff = ff +  accumarray(elemNv(:),F(:),[N 1]);
    
    % elementwise global indices
    elem2dof(idNv) = mat2cell(elemNv, ones(NTv,1), Ndof);
end
kk = sparse(ii,jj,ss,N,N);

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
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
% solver = 'amg';
% if N < 2e3, solver = 'direct'; end
% % solve
% switch solver
%     case 'direct'
%         u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
%     case 'amg'
%         option.solver = 'CG';
%         u(freeDof) = amg(kk(freeDof,freeDof),ff(freeDof),option);                 
% end

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm
info.D = D;