function [u,info] = PoissonVEM_VI(node,elem,pde,bdStruct)
%PoissonVEM_VI solves the simplified friction problem using conforming
% virtual element method in the lowest order case.
%
% The problem is
%
%  D.E.
%     -\Delta u + cu = f,  in Omega,
%  Frictional boundary conditions
%     |grad(u)*n| <= g,  u*grad(u)*n + g*|u| = 0 on \Gamma_C
%  Dirichlet boundary conditions
%     u = g_D    on \Gamma_D
%  \partial(\Omega) = \Gamma_C + \Gamma_D
%
%   References
%   F. Feng, W. M. Han and J. G. Huang, "Virtual element methods for
%   elliptic variational inequalities of the second kind", J. Sci. Comput., 
%   Vol 80. No 1., pp. 60¨C80, 2019.  (See Eq. (4.6) there)
%
% Copyright (C)  Terence Yu.

%% Get auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid;  diameter = aux.diameter;  area = aux.area;
N = size(node,1);  NT = size(elem,1);

%% Compute projection matrices
D = cell(NT,1);  Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1); 
H = cell(NT,1);
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);    
    xK = centroid(iel,1); yK = centroid(iel,2); 
    hK = diameter(iel);
    x = node(index,1); y = node(index,2);    
    % scaled monomials
    m1 = @(x,y) 1 + 0*x; m2 = @(x,y) (x-xK)./hK; m3 = @(x,y) (y-yK)./hK;    
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];   
    Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]];
    % D
    D1 = m(x,y);   D{iel} = D1;
    % B, Bs, G, Gs
    rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices ( z_{i-1}z_{i+1} )
    normVec = 0.5*[y(rotid2)-y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
    B1 = Gradm*normVec; % B
    B1s = B1; B1s(1,:) = 1/Nv;    Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;  
    % H
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];   
    % mm 
    mm = @(x,y) [m1(x,y).*m1(x,y), m1(x,y).*m2(x,y), m1(x,y).*m3(x,y), ...
                 m2(x,y).*m1(x,y), m2(x,y).*m2(x,y), m2(x,y).*m3(x,y), ...
                 m3(x,y).*m1(x,y), m3(x,y).*m2(x,y), m3(x,y).*m3(x,y)];
    H1 = zeros(3,3);
    H1(:) = integralTri(mm,3,nodeT,elemT); % n = 3   
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

%% Assemble frictional boundary conditions
bdEdgeFri = bdStruct.bdEdgeN;
if ~isempty(bdEdgeFri)
    g = pde.g_C;
    z1 = node(bdEdgeFri(:,1),:); z2 = node(bdEdgeFri(:,2),:);
    gC = max([g(z1); g(z2)]);
    e = z2-z1;
    he = sqrt(e(:,1).^2+e(:,2).^2);
    Fri = 0.5*gC*[he,he];
    w = accumarray(bdEdgeFri(:), Fri(:),[N 1]);
end

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;  bdNodeIdx = bdStruct.bdNodeIdx;
isBdNode = false(N,1); isBdNode(bdNodeIdx) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node(bdNodeIdx,:);
u = zeros(N,1); uD = g_D(nodeD); u(bdDof) = uD(:);

%% Get minimization problem
%fun = @(v)  0.5*v'*A*v - b'*v + w'*abs(v);
A = speye(N,N);
A(freeDof,freeDof) = kk(freeDof,freeDof);
b = u;
b(freeDof) = ff(freeDof) - kk(freeDof,freeDof)*u(freeDof);

%% Set solver
% Note: It can be directly solved by: u = cvx2cvx(A,b,w)
% schur 
idFri = unique(bdEdgeFri(:));
idFriRes = setdiff(1:N,idFri); 
A11 = A(idFri,idFri);  
A12 = A(idFri,idFriRes);
A22 = A(idFriRes,idFriRes);
b1 = b(idFri);  b2 = b(idFriRes);
AA = A11-A12*(A22\A12');
bb = b1-A12*(A22\b2);
ww = w(idFri);
% Fixed point algorithm based on proximity operator
v = cvx2cvx(AA,bb,ww);
vRes = A22\(b2-A12'*v);
u = zeros(N,1);
u(idFri) = v;
u(idFriRes) = vRes;

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; info.freeDof = freeDof;