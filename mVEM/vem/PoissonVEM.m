function u = PoissonVEM(node,elem,pde,bdStruct)

aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
elemCentroid = aux.elemCentroid; area = aux.area;

N = size(node,1);  NT = size(elem,1);
% ----------------- D,Bs,G,Gs,H ---------------------------
[D,Bs,G,Gs,H] = VEM_MAT_Poisson(aux);

% -------- Elementwise stiffness matrix and load vector -----------
ABelem = cell(NT,1); belem = cell(NT,1);
for iel = 1:NT
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));    
    % Stiffness matrix
    AK  = Pis'*G{iel}*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pis'*H{iel}*Pis;  % Pis = Pi0s;
    ABelem{iel} = reshape(AK+BK,1,[]); % straighten as row vector for easy assembly    
    % Load vector   
    fK = Pis'*[pde.f(elemCentroid(iel,:))*area(iel);0;0];    
    belem{iel} = fK'; % straighten as row vector for easy assembly
end

% --------- Assemble the stiff matrix and load vector ---------
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnz = sum(elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);

id = 0; ff = zeros(N,1); 
for Nv = vertNum(:)' % only valid for row vector
    
    % assemble the matrix
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices
    
    elemNv = cell2mat(elem(idNv)); % elem     
    K = cell2mat(ABelem(idNv)); F = cell2mat(belem(idNv));
    s = 1; Ndof = Nv;
    for i = 1:Ndof
        for j = 1:Ndof
            ii(id+1:id+NTv) = elemNv(:,i); % zi
            jj(id+1:id+NTv) = elemNv(:,j); % zj
            ss(id+1:id+NTv) = K(:,s);
            id = id + NTv; s = s+1;
        end
    end
    
    % assemble the vector
    ff = ff +  accumarray(elemNv(:),F(:),[N 1]);
end
kk = sparse(ii,jj,ss,N,N);

% ------------ Neumann boundary condition ----------------
elemN = bdStruct.elemN;
if ~isempty(elemN)
    g_N = pde.g_N;
    z1 = node(elemN(:,1),:); z2 = node(elemN(:,2),:);
    e = z1-z2;  % e = z2-z1
    ne = [-e(:,2),e(:,1)]; % scaled ne
    F1 = sum(ne.*g_N(z1),2)./2;
    F2 = sum(ne.*g_N(z2),2)./2;
    FN = [F1,F2];
    ff = ff + accumarray(elemN(:), FN(:),[N 1]);
end

% ------------ Dirichlet boundary condition ----------------
g_D = pde.g_D;  eD = bdStruct.eD;
isBdNode = false(N,1); isBdNode(eD) = true;
bdNode = find(isBdNode); freeNode = find(~isBdNode);
pD = node(bdNode,:);
u = zeros(N,1); uD = g_D(pD); u(bdNode) = uD(:);
ff = ff - kk*u;

% ------------------ Solver -------------------
u(freeNode) = kk(freeNode,freeNode)\ff(freeNode);
