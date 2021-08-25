function [w,info] = PlateBending_C1VEM(node,elem,pde,bdStruct,option)
%PlateBending_C0VEM solves plate bending problem using C1-VEM
%
%       -D_{ij} M_{ij}(w) = f in \Omega,
%       Dirichlet boundary condition:
%               w = g1, grad(w)n = g2    on \Gamma.
%   References
%   F. Brezzi and L.D. Marini, "Virtual Element Methods for plate being
%   problems", Comput. Methods Appl. Mech. Engrg., Vol 253., pp. 455¨C462, 2013.
%
%  Note: The problem under consideration is in the general form.
%
% Copyright (C)  Terence Yu. 

%% Get auxiliary data
para = pde.para;
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; 
% auxstructure
auxT = auxstructure(node,elem);
node2elem = auxT.node2elem;
% number
N = size(node,1); NT = size(elem,1); NNdof = 3*N;   Nm = 6;
% characteristic lengths of derivatives at nodes
hxi = cellfun(@(id) mean(diameter(id)), node2elem);

%% Compute projection matrices
D = cell(NT,1);  Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1);
for iel = 1:NT
    % element information
    index = elem{iel};     Nv = length(index);     Ndof = 3*Nv;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    hxiK = hxi(index); % characteristic length 
    % scaled monomials
    % m'
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
                (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    % Dx(m'), Dy(m')
    Dxm = @(x,y) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x];
    Dym = @(x,y) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2];
    % D 
    D1 = zeros(Ndof,Nm);
    D1(1:Nv,:) = m(x,y);
    D1(Nv+1:2*Nv,:) = repmat(hxiK,1,Nm).*Dxm(x,y); 
    D1(2*Nv+1:end,:) = repmat(hxiK,1,Nm).*Dym(x,y);
    D{iel} = D1;
    % \partial_ij (m)
    D11 = zeros(Nm,1); D11(4) = 2/hK^2;
    D12 = zeros(Nm,1); D12(5) = 1/hK^2;
    D22 = zeros(Nm,1); D22(6) = 2/hK^2;
    % Mij(m)
    M11 = -para.D*((1-para.nu)*D11 + para.nu*(D11+D22));
    M12 = -para.D*(1-para.nu)*D12;
    M22 = -para.D*((1-para.nu)*D22 + para.nu*(D11+D22));
    % Mnn(m) on e1,...,eNv
    n1 = ne(:,1);  n2 = ne(:,2);
    Mnn = M11*(n1.*n1)' + M12*(n1.*n2+n2.*n1)' + M22*(n2.*n2)';
    % Mtn(m) on e1,...,eNv
    t1 = te(:,1); t2 = te(:,2);
    Mtn = M11*(t1.*n1)' + M12*(t1.*n2+t2.*n1)' + M22*(t2.*n2)';
    % Dx(phi'), Dy(phi') at z1,...,zNv (each row)
    Dxphi = zeros(Nv,Ndof); Dyphi = zeros(Nv,Ndof);
    Dxphi(:,Nv+1:2*Nv) = eye(Nv)./repmat(hxiK,1,Nv);  
    Dyphi(:,2*Nv+1:end) = eye(Nv)./repmat(hxiK,1,Nv);
    % B, Bs, G, Gs
    B1 = zeros(Nm,Ndof);
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    for j = 1:Nv % loop of edges
        % int[\partial_n (phi')] on ej
        Dnphi1 = Dxphi(v1(j),:)*n1(j) + Dyphi(v1(j),:)*n2(j); % zj
        Dnphi2 = Dxphi(v2(j),:)*n1(j) + Dyphi(v2(j),:)*n2(j); % z_{j+1}
        nphi = 0.5*he(j)*(Dnphi1+Dnphi2);
        % Jump(m) at zj
        Jump = Mtn(:,p2(j))-Mtn(:,p1(j));
        % phi' at zj
        phi = zeros(1,Ndof);  phi(j) = 1;
        % B1 on e and at zj
        B1 = B1 - Mnn(:,j)*nphi + Jump*phi;
    end
    B1s = B1;
    % first constraint
    B1s(1,1:Nv) = 1/Nv;
    % second constraint
    for j = 1:Nv % loop of edges
        B1s(2:3,1:Nv) = te([Nv,1:Nv-1],:)' - te';
        Dnphi1 = Dxphi(v1(j),Nv+1:end)*n1(j) + Dyphi(v1(j),Nv+1:end)*n2(j); % zj
        Dnphi2 = Dxphi(v2(j),Nv+1:end)*n1(j) + Dyphi(v2(j),Nv+1:end)*n2(j); % z_{j+1}
        Nphi = 0.5*Ne(j,:)'*(Dnphi1+Dnphi2); % scaled        
        B1s(2:3,Nv+1:end) = B1s(2:3,Nv+1:end) + Nphi;
    end
    Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;
end

%% Get elementwise stiffness matrix and load vector
Aelem = cell(NT,1); belem = cell(NT,1);
Ph = cell(NT,1); % matrix for H1 and L2 error
for iel = 1:NT
    % element information
    index = elem{iel};     Nv = length(index);  Ndof = 3*Nv;
    hxiK = hxi(index); % characteristic length 
    hh = diag(repmat(hxiK.^(-2), 3, 1)); 
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi = D{iel}*Pis;   I = eye(Ndof);
    % Stiffness matrix
    %hK = diameter(iel);
    AK  = Pis'*G{iel}*Pis + (hh*(I-Pi))'*(I-Pi);
    Aelem{iel} = reshape(AK',1,[]); % straighten
    % Load vector
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
        (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,3,nodeT,elemT); rhs = rhs';
    fK = Pis'*rhs;
    belem{iel} = fK'; % straighten
    % matrix for H1 an L2 error
    Ph{iel} = Pis; 
end

%% Assemble stiffness matrix and load vector
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnz = sum((3*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
id = 0;   ff = zeros(NNdof,1);
elem2dof = cell(NT,1);
for Nv = vertNum(:)' % only valid for row vector
    
    % assemble the matrix
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices
    
    elemNv = cell2mat(elem(idNv)); % elem
    elem2 = [elemNv, elemNv+N, elemNv+2*N];
    K = cell2mat(Aelem(idNv)); F = cell2mat(belem(idNv));
    
    Ndof = 3*Nv;
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

%% Apply Dirichlet boundary conditions
% boundary information 
bdNodeIdx = bdStruct.bdNodeIdx; 
hxib = hxi(bdNodeIdx);% characteristic length on the boundary
% bdDof, freeDof
id = [bdNodeIdx; bdNodeIdx+N; bdNodeIdx+2*N];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% values on the boundary
nodeD = node(bdNodeIdx,:); g_D = pde.g_D;  wD = g_D(nodeD);  
% derivatives on the boundary with characteristic length
Dw = pde.Du;  
hDw = Dw(nodeD).*repmat(hxib,1,2);   % wxD = Dw(:,1); wyD = Dw(:,2);
w = zeros(NNdof,1); w(bdDof) = [wD; hDw(:)];
ff = ff - kk*w;

%% Set solver
% solver = option.solver;
% if NNdof < 2e3, solver = 'direct'; end
% solve
solver = 'direct';
switch solver
    case 'direct'
        w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
    case 'amg'
        option.solver = 'CG';
        w(freeDof) = amg(kk(freeDof,freeDof),ff(freeDof),option);                 
end

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; %info.freeDof = freeDof;
