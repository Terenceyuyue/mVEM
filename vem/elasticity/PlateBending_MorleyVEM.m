function [w,info] = PlateBending_MorleyVEM(node,elem,pde,bdStruct)
%PlateBending_MorleyVEM solves plate bending problem using Morley-VEM
%
%       -D_{ij} M_{ij}(w) = f in \Omega,
%       Dirichlet boundary condition:
%               w = g1, grad(w)n = g2    on \Gamma.
%   References
%   J. Zhao, B. Zhang and S. Chen, "The Morley-Type Virtual Element
%   for Plate Bending Problems", J. Sci. Comput.,
%   Vol 76., pp. 610¨C629, 2018.
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
NNdof = N+NE;   Nm = 6;

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(2*elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
for iel = 1:NT
    % --------- element information ---------
    index = elem{iel};  indexEdge = elem2edge{iel};
    Nv = length(index);     Ndof = 2*Nv;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    
    % ------- scaled monomials ---------
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
        (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    Gradm = @(x,y) [[0,0]; [1, 0]/hK; [0, 1]/hK; [2*(x-xK), 0]/hK^2;
        [(y-yK), (x-xK)]/hK^2; [0, 2*(y-yK)]/hK^2]; % grad m
    
    % ------ transition matrix ---------
    D = zeros(Ndof,Nm);
    D(1:Nv,:) = m(x,y);  % vertices
    for i = 1:Nv
        gradi = 0.5*(Gradm(x(v1(i)),y(v1(i)))+Gradm(x(v2(i)),y(v2(i))));
        D(Nv+i,:) = Ne(i,:)*gradi';
    end
    
    % --------- elliptic projection -----------
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
    B = zeros(Nm,Ndof);
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    for j = 1:Nv % loop of edges
        % nphi on ej
        nphi = zeros(1,Ndof); nphi(Nv+j) = 1;  
        % Jump at zj
        Jump = Mtn(:,p2(j))-Mtn(:,p1(j));
        % phi at zj
        phi = zeros(1,Ndof); phi(j) = 1;    
        % B1
        B = B - Mnn(:,j)*nphi + Jump*phi;
    end
    Bs = B;
    % first constraint
    Bs(1,1:Nv) = 1/Nv;
    % second constraint
    Bs(2:3,1:Nv) = te([Nv,1:Nv-1],:)' - te';
    Bs(2:3,Nv+1:end) = ne';
    % consistency relation
    G = B*D;  Gs = Bs*D;
    
    % -------- sgnBase ----------------
    sgnedge = sign(diff(index([1:Nv,1])));
    id = indexEdge; sgnbd = E(id); sgnedge(sgnbd) = 1;
    sgnBase = ones(Ndof,1); sgnBase(Nv+1:end) = sgnedge;
    sgnK = sgnBase*sgnBase';  sgnF = sgnBase;
    
    % --------- local stiffness matrix ---------     
    Pis = Gs\Bs;   Pi = D*Pis;   I = eye(Ndof);
    AK  = Pis'*G*Pis + hK^(-2)*(I-Pi)'*(I-Pi);
    AK = AK.*sgnK;  A = reshape(AK',1,[]); 
    
    % --------- load vector -----------
%   cK = centroid(iel,:); areaK = aux.area(iel);
%   P0 = zeros(Ndof,1); P0(1:Nv) = 1/Nv;
%   fK = pde.f(cK)*areaK*P0;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
        (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,2,nodeT,elemT); rhs = rhs';
    fK = Pis'*rhs;  fK = fK.*sgnF;
    
    % --------- assembly index for ellptic projection -----------
    indexDof = [index, indexEdge+N];  
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = A(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    sgnPis = repmat(sgnBase',size(Pis,1),1);
    Ph{iel} = sgnPis.*Pis; 
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,NNdof,NNdof);
ff = accumarray(elemb,Fb,[NNdof 1]);

%% Apply Dirichlet boundary conditions
bdNodeIdx = bdStruct.bdNodeIdx; 
bdEdgeD = bdStruct.bdEdgeD;
g_D = pde.g_D;  Dw = pde.Du;

id = [bdNodeIdx; bdEdgeIdx+N];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
nodeD = node(bdNodeIdx,:); wD = g_D(nodeD);
% moments on the boundary
z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:); zc = (z1+z2)./2;
eb = z1-z2;  % e = z2-z1
neb = [-eb(:,2),eb(:,1)];  %scaled ne
wnD = sum(1/6*(Dw(z1)+4*Dw(zc)+Dw(z2)).*neb,2);
% rhs
w = zeros(NNdof,1); w(bdDof) = [wD; wnD];
ff = ff - kk*w;

%% Set solver
% solver = 'amg';
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
info.kk = kk; %info.DofI = freeDof;