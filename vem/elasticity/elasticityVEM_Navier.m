function [u,info] = elasticityVEM_Navier(node,elem,pde,bdStruct)
%elasticityVEM_displacement solves linear elasticity equation of Navier
%form using virtual element method in V1
%
%       u = [u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%
% Copyright (C)  Terence Yu. 

%% Get auxiliary data
para = pde.para;
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
N = size(node,1);  NT = size(elem,1);

%% Compute projection matrices
D = cell(NT,1);
% B = cell(NT,1);  % it is not used in the computation
Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1);
H0 = cell(NT,1); C0 = cell(NT,1);
for iel = 1:NT
    % element information
    index = elem{iel};     Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    % scaled monomials
    m1 = @(x,y) 1 + 0*x; m2 = @(x,y) (x-xK)./hK; m3 = @(x,y) (y-yK)./hK;
        
    % D
    D1 = [m1(x,y), m2(x,y), m3(x,y)];   
    D{iel} = D1;
    
    % B, Bs, G, Gs
    rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
    Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]];
    normVec = 0.5*[y(rotid2) - y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
    B1 = Gradm*normVec; B1s = B1; B1s(1,:) = 1/Nv;
    % B{iel} = B1;
    Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;
    
    % H0,C0
    H0{iel} = area(iel);
    C0{iel} = reshape(normVec',1,[]);
end

%% Get elementwise stiffness matrix and load vector
ABelem = cell(NT,1); belem = cell(NT,1);
Ph = cell(NT,1); % matrix for error evaluation
for iel = 1:NT
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));
    Pi0s = H0{iel}\C0{iel};
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi); AK = para.mu*blkdiag(AK,AK);
    BK = (para.mu+para.lambda)*Pi0s'*H0{iel}*Pi0s;
    ABelem{iel} = reshape(AK'+BK',1,[]); % straighten
    % Load vector
    Nv = length(elem{iel});
    fK = pde.f(centroid(iel,:))*area(iel)/Nv; fK = repmat(fK,Nv,1);
    belem{iel} = fK(:)'; % straighten
    % matrix for L2 and H1 error evaluation
    Ph{iel} = blkdiag(Pis,Pis); 
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

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;
bdNodeIdx = bdStruct.bdNodeIdx;
isBdNode = false(2*N,1); isBdNode([bdNodeIdx;bdNodeIdx+N]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdNodeIdx,:);
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