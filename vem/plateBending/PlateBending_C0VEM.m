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
%   Vol 26. No 9., pp. 1671�C1687, 2016.
%
%  Note:
%  - The problem under consideration is in the general form.
%  - In the lowest order case k = 2: The right-hand side is approximated by
%  the L^2 projection, which is exactly the H^2 elliptic projection when we
%  consider the enhancement technique.
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
node2elem = auxT.node2elem;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+2*NE;   Nm = 6;
% characteristic lengths at nodes
hxi = cellfun(@(id) mean(diameter(id)), node2elem);

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  
elemLen = cellfun('length',elem); 
nnz = sum((3*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(3*elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
for iel = 1:NT
    % --------- element information ---------
    index = elem{iel};  indexEdge = elem2edge{iel};
    Nv = length(index);     Ndof = 3*Nv;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    hxiK = hxi(index); % characteristic length 
    
    % ------- scaled monomials ---------
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
        (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    Gradm = @(x,y) [[0*x,0*x]; [1+0*x, 0*x]/hK; [0*x, 1+0*x]/hK; ...
        [2*(x-xK), 0*x]/hK^2;[(y-yK), (x-xK)]/hK^2; [0*x, 2*(y-yK)]/hK^2]; % grad m
    
    % ------ transition matrix ---------
    D = zeros(Ndof,Nm);
    D(1:2*Nv,:) = [m(x,y); m(xe,ye)];  % vertices and mid-edge points
    for i = 1:Nv  % moments of normal derivatives on edges
        gradi = 0.5*(Gradm(x(v1(i)),y(v1(i)))+Gradm(x(v2(i)),y(v2(i))));
        D(2*Nv+i,:) = Ne(i,:)*gradi';
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
    % B, Bs, G, Gs
    B = zeros(Nm,Ndof);
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    for i = 1:Nv % loop of edges or vertices
        % Mnn(m)
        n1 = ne(i,1);  n2 = ne(i,2);
        Mnn = M11*(n1*n1) + M12*(n1*n2+n2*n1) + M22*(n2*n2); % (Nm,1)
        % nphi on ei
        nphi = zeros(1,Ndof); nphi(2*Nv+i) = 1;
        % Jump at zi
        tn11 = te(p2(i),1)*ne(p2(i),1) - te(p1(i),1)*ne(p1(i),1); % jump
        tn22 = te(p2(i),2)*ne(p2(i),2) - te(p1(i),2)*ne(p1(i),2);
        tn12 = (te(p2(i),1)*ne(p2(i),2) + te(p2(i),2)*ne(p2(i),1)) ...
            - (te(p1(i),1)*ne(p1(i),2) + te(p1(i),2)*ne(p1(i),1));
        Jump = M11*tn11 + M12*tn12 + M22*tn22;
        % phi at zi
        phi = zeros(1,Ndof); phi(i) = 1;
        % B1
        B = B - Mnn*nphi + Jump*phi;
    end
    Bs = B;
    % first constraint
    Bs(1,1:Nv) = 1/Nv; 
    % second constraint
    Bs(2:3,1:Nv) = te(p1,:)' - te(p2,:)';
    Bs(2:3,2*Nv+1:end) = ne';
    % consistency relation
    G = B*D;     Gs = Bs*D;
    
    % -------- sgnBase ----------------
    sgnedge = sign(diff(index([1:Nv,1])));
    id = indexEdge; sgnbd = E(id); sgnedge(sgnbd) = 1;
    sgnBase = ones(Ndof,1); sgnBase(2*Nv+1:end) = sgnedge;
    sgnK = sgnBase*sgnBase';  sgnF = sgnBase;
    
    % --------- local stiffness matrix --------- 
    hh = diag([hxiK; he; he].^(-2));
    Pis = Gs\Bs;   Pi = D*Pis;   I = eye(Ndof);
    AK  = Pis'*G*Pis + (hh*(I-Pi))'*(I-Pi);
    AK = AK.*sgnK;  A = reshape(AK',1,[]); 
    
    % --------- load vector -----------
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
        (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,5,nodeT,elemT); rhs = rhs';
    fK = Pis'*rhs;  fK = fK.*sgnF;
    
    % --------- assembly index for ellptic projection -----------
    indexDof = [index, indexEdge+N, indexEdge+N+NE];  
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
% boundary information
bdNodeIdx = bdStruct.bdNodeIdx; 
bdEdge = bdStruct.bdEdge;
g_D = pde.g_D;  Dw = pde.Du;
id = [bdNodeIdx; bdEdgeIdx+N; bdEdgeIdx+N+NE];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
nodeD = node(bdNodeIdx,:); wD = g_D(nodeD);
% mid-edge on the boundary
z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:); zc = (z1+z2)./2;
wDc = g_D(zc);
% moments on the boundary
eb = z1-z2;  % e = z2-z1
neb = [-eb(:,2),eb(:,1)];  % scaled ne
wnD = sum(1/6*(Dw(z1)+4*Dw(zc)+Dw(z2)).*neb,2);
% rhs
w = zeros(NNdof,1); w(bdDof) = [wD; wDc; wnD];
ff = ff - kk*w;

%% Set solver
if isfield(bdStruct,'solver')
    solver = bdStruct.solver;
else
    solver = 'direct';
end
if NNdof < 1e4, solver = 'direct'; end
switch solver
    case 'direct'
        w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
    case 'amg'
        option.solver = 'CG';
        w(freeDof) = amg(kk(freeDof,freeDof),ff(freeDof),option);    
    case 'cg'
        w(freeDof) = cgs(kk(freeDof,freeDof),ff(freeDof),1e-12,NNdof);
end

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; info.DofI = freeDof;