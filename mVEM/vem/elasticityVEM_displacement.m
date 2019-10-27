function u = elasticityVEM_displacement(node,elem,pde,bdStruct)

aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
elemCentroid = aux.elemCentroid; area = aux.area;

N = size(node,1);  NT = size(elem,1);
% ----------------- D,Bs,G,Gs,H0,C0 ---------------------------
[D,Bs,G,Gs,H0,C0] = VEM_MAT_elasticity_displacement(aux);

% -------- Elementwise stiffness matrix and load vector -----------
ABelem = cell(NT,1); belem = cell(NT,1);
for iel = 1:NT
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));
    Pi0s = H0{iel}\C0{iel};
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi); AK = pde.mu*blkdiag(AK,AK);
    BK = (pde.mu+pde.lambda)*Pi0s'*H0{iel}*Pi0s;
    ABelem{iel} = reshape(AK+BK,1,[]); % straighten
    % Load vector
    Nv = length(elem{iel});
    fK = pde.f(elemCentroid(iel,:))*area(iel)/Nv; fK = repmat(fK,Nv,1);
    belem{iel} = fK(:)'; % straighten
end

% --------- Assemble the stiffness matrix and load vector ---------
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnz = sum(4*elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);

id = 0; ff = zeros(2*N,1);
for Nv = vertNum(:)' % only valid for row vector
    
    % assemble the matrix
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices
    
    elemNv = cell2mat(elem(idNv)); % elem
    K = cell2mat(ABelem(idNv)); F = cell2mat(belem(idNv));
    s = 1; Ndof = Nv;
    for i = 1:2*Ndof
        for j = 1:2*Ndof
            iN = i>Ndof; jN = j>Ndof;
            i1 = i-iN*Ndof; j1 = j-jN*Ndof;
            ii(id+1:id+NTv) = elemNv(:,i1) + iN*N; % zi
            jj(id+1:id+NTv) = elemNv(:,j1) + jN*N; % zj
            ss(id+1:id+NTv) = K(:,s);
            id = id + NTv; s = s+1;
        end
    end
    
    % assemble the vector
    ff = ff +  accumarray([elemNv(:);elemNv(:)+N],F(:),[2*N 1]);
end
kk = sparse(ii,jj,ss,2*N,2*N);

% ------------ Dirichlet boundary condition ----------------
g_D = pde.g_D;  eD = bdStruct.eD;
isBdNode = false(N,1); isBdNode(eD) = true;
bdNode = find(isBdNode); freeNode = find(~isBdNode);
pD = node(bdNode,:);
bdDof = [bdNode; bdNode+N]; freeDof = [freeNode;freeNode+N];
u = zeros(2*N,1); uD = g_D(pD); u(bdDof) = uD(:);
ff = ff - kk*u;

% ------------------ Solver -------------------
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);