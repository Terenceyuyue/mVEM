function [w,info] = PlateBending_C0VEM(node,elem,pde,bdStruct)
%PlateBending_C0VEM solves plate bending problem using C0-VEM
%
%       -D_{ij} M_{ij}(w) = f in \Omega,
%       Dirichlet boundary condition:
%               w = g1, grad(w)n = g2    on \Gamma.
%
%   References
%   J. Zhao, S. Chen and B. Zhang, "The nonconforming virtual element
%   method for plate bending problems", Math. Models Meth. Appl. Sci.,
%   Vol 26. No 9., pp. 1671¨C1687, 2016.
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
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+2*NE;   Nm = 6;

%% Compute projection matrices
D = cell(NT,1);  Bs = cell(NT,1);
G = cell(NT,1);  Gs = cell(NT,1);
for iel = 1:NT
    % element information
    index = elem{iel};     Nv = length(index);     Ndof = 3*Nv;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    % scaled monomials
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
        (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    Gradm = @(x,y) [[0+0*x,0+0*x]; [1+0*x, 0+0*x]/hK; [0+0*x, 1+0*x]/hK; ...
        [2*(x-xK), 0+0*x]/hK^2;[(y-yK), (x-xK)]/hK^2; [0+0*x, 2*(y-yK)]/hK^2]; % grad m
    % D
    D1 = zeros(Ndof,Nm);
    D1(1:2*Nv,:) = [m(x,y); m(xe,ye)];  % vertices and mid-edge points
    for i = 1:Nv  % normal moments on edges
        gradi = 0.5*(Gradm(x(v1(i)),y(v1(i)))+Gradm(x(v2(i)),y(v2(i))));
        D1(2*Nv+i,:) = Ne(i,:)*gradi';
    end
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
    % B, Bs, G, Gs
    B1 = zeros(Nm,Ndof);
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    for j = 1:Nv % loop of edges
        % nphi on ej
        nphi = zeros(1,Ndof); nphi(2*Nv+j) = 1;
        % Jump at zj
        Jump = Mtn(:,p2(j))-Mtn(:,p1(j));
        % phi at zj
        phi = zeros(1,Ndof); phi(j) = 1;
        % B1
        B1 = B1 - Mnn(:,j)*nphi + Jump*phi;
    end
    B1s = B1;
    % first constraint
    B1s(1,1:Nv) = 1/Nv; 
    % second constraint
    B1s(2:3,1:Nv) = te([Nv,1:Nv-1],:)' - te';
    B1s(2:3,2*Nv+1:end) = ne';
    Bs{iel} = B1s;
    G{iel} = B1*D1;     Gs{iel} = B1s*D1;
end

%% Get elementwise signs of basis functions
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
sgnBase = cell(NT,1);
for iel = 1:NT
    index = elem{iel};     Nv = length(index);  Ndof = 3*Nv;
    sgnedge = sign(diff(index([1:Nv,1])));
    id = elem2edge{iel}; sgnbd = E(id); sgnedge(sgnbd) = 1;
    sgnelem = ones(Ndof,1); sgnelem(2*Nv+1:end) = sgnedge;
    sgnBase{iel} = sgnelem;
end

%% Get elementwise stiffness matrix and load vector
Aelem = cell(NT,1); belem = cell(NT,1);
Ph = cell(NT,1); % matrix for H1 and L2 error
for iel = 1:NT
    % sign matrix and sign vector
    index = elem{iel};     Nv = length(index);  Ndof = 3*Nv;
    sgnelem = sgnBase{iel};
    sgnK = sgnelem*sgnelem';  sgnF = sgnelem;
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi = D{iel}*Pis;   I = eye(Ndof);
    % Stiffness matrix
    hK = diameter(iel);
    AK  = Pis'*G{iel}*Pis + hK^(-2)*(I-Pi)'*(I-Pi);
    AK = AK.*sgnK;
    Aelem{iel} = reshape(AK',1,[]); % straighten
    % Load vector
%     cK = centroid(iel,:); areaK = aux.area(iel);
%     P0 = zeros(Ndof,1); P0(1:Nv) = 1/Nv;
%     fK = pde.f(cK)*areaK*P0;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
        (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,2,nodeT,elemT); rhs = rhs';
    fK = Pis'*rhs;
    fK = fK.*sgnF;
    belem{iel} = fK'; % straighten
    % matrix for H1 an L2 error
    sgnPis = repmat(sgnelem',size(Pis,1),1);
    Ph{iel} = sgnPis.*Pis; 
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
    elem2edgeNv = cell2mat(elem2edge(idNv));
    elem2 = [elemNv, elem2edgeNv+N, elem2edgeNv+N+NE];
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
bdEdgeD = bdStruct.bdEdgeD;
g_D = pde.g_D;  Dw = pde.Du;
id = [bdNodeIdx; bdEdgeIdx+N; bdEdgeIdx+N+NE];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
nodeD = node(bdNodeIdx,:); wD = g_D(nodeD);
% mid-edge on the boundary
z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:); zc = (z1+z2)./2;
wDc = g_D(zc);
% moments on the boundary
eb = z1-z2;  % e = z2-z1
neb = [-eb(:,2),eb(:,1)];  % scaled ne
wnD = sum(1/6*(Dw(z1)+4*Dw(zc)+Dw(z2)).*neb,2);
% rhs
w = zeros(NNdof,1); w(bdDof) = [wD; wDc; wnD];
ff = ff - kk*w;

%% Set solver
% solver = 'amg';
% if NNdof < 2e3, solver = 'direct'; end
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
info.kk = kk; info.freeDof = freeDof;
