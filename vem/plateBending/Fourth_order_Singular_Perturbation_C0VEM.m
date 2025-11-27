function [w,info] = Fourth_order_Singular_Perturbation_C0VEM(node,elem,pde,bdStruct)
%This function solves fourth-order singular perturbation using C0-VEM
%
%     \epsilon^2 \Delta^2 u -\Delta u = f,  in Omega
%     Dirichlet boundary condition: 
%                    u = g1,  grad(u)n = g2 on \Gamma, 
%
%  Note:
%  - The problem under consideration is in the general form.
%  - In the lowest order case k = 2: The right-hand side is approximated by
%  the L^2 projection, which is exactly the H^2 elliptic projection when we
%  consider the enhancement technique.
%
% Copyright (C)  Terence Yu. 
%

%% Get auxiliary data
para = pde.para;
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
node2elem = auxT.node2elem;
elem2sgn = auxT.elem2sgn;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+2*NE;  Nm = 6;
% characteristic lengths at nodes
hxi = cellfun(@(id) mean(diameter(id)), node2elem);

%% Compute projection matrices
% Bbs, Gb, Gbs % b: biharmonic
% Bps, Gp, Gps % p: Poisson
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  
elemLen = cellfun('length',elem); 
nnz = sum((3*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(3*elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % ------------- element information ----------------
    index = elem{iel};   indexEdge = elem2edge{iel};
    Nv = length(index);  Ndof = 3*Nv; 
    xK = centroid(iel,1); yK = centroid(iel,2);  hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv; % jump index
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    hxiK = hxi(index); % characteristic length 

    % ------------ scaled monomials ---------------------
    m1 = @(x,y) 1+0*x;                    gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;               gradm2 = @(x,y) [1+0*x, 0+0*x]./hK;
    m3 = @(x,y) (y-yK)./hK;               gradm3 = @(x,y) [0+0*x, 1+0*x]./hK;
    m4 = @(x,y) (x-xK).^2/hK^2;           gradm4 = @(x,y) [2*(x-xK), 0+0*x]./hK^2;
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2;     gradm5 = @(x,y) [(y-yK), (x-xK)]./hK^2;
    m6 = @(x,y) (y-yK).^2./hK^2;          gradm6 = @(x,y) [0+0*x, 2*(y-yK)]./hK^2;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)]; 
    Gradm = @(x,y) [gradm1(x,y); gradm2(x,y); gradm3(x,y); ...
        gradm4(x,y); gradm5(x,y); gradm6(x,y)]; 
    Gradmc = {gradm1, gradm2, gradm3, gradm4, gradm5, gradm6};
    
    % --------------- transition matrices --------------
    D = zeros(Ndof,Nm);
    D(1:2*Nv,:) = [m(x,y); m(xe,ye)];  % vertices and mid-edge points
    for i = 1:Nv
        gradi = 0.5*(Gradm(x(v1(i)),y(v1(i)))+Gradm(x(v2(i)),y(v2(i))));
        D(2*Nv+i,:) = Ne(i,:)*gradi';
    end
    
    % 2. ---------- H2 elliptic projection ----------
    % \partial_ij (m)
    D11 = zeros(6,1); D11(4) = 2/hK^2;
    D12 = zeros(6,1); D12(5) = 1/hK^2;
    D22 = zeros(6,1); D22(6) = 2/hK^2;
    % Mij(m)
    nu = 0;  % D = 1, nu = 0 for biharmonic
    M11 = -((1-nu)*D11 + nu*(D11+D22));
    M12 = -(1-nu)*D12;
    M22 = -((1-nu)*D22 + nu*(D11+D22));
    % B, Bs, G, Gs
    B2 = zeros(Nm,Ndof);
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
        B2 = B2 - Mnn*nphi + Jump*phi;
    end
    B2s = B2;
    % first constraint
    B2s(1,1:Nv) = 1/Nv; 
    % second constraint
    B2s(2:3,1:Nv) = te(p1,:)' - te(p2,:)';
    B2s(2:3,2*Nv+1:end) = ne';
    % consistency relation
    G2 = B2*D;     G2s = B2s*D;

    % ----------------- L2 projection -----------------------
    % transition matrix: lifting   
    DL = zeros(Ndof+1,Nm);
    DL(1:Ndof,:) = D;    
    DL(end,:) = 1/area(iel)*integralTri(m,5,nodeT,elemT); %  interior d.o.f
    % B,Bs,G,Gs
    BLs = zeros(Nm,Ndof+1);  BLs(:,1:Ndof) = B2s;
    % G, Gs:  G = G2, GLs = G2s
    % Pi, Pis
    PLs = G2s\BLs;   PiL = DL*PLs;
    Dof = PiL(end,1:end-1); 
    % L2 projection
    C = area(iel)*Dof;
    H = C*D;
    Pi0s = H\C;
    
    % --------- H1 elliptic projection ----------
    % first term
    Lapm = zeros(Nm,1); Lapm([4,6]) = 2/hK^2;
    I1 = Lapm*area(iel)*Dof;
    % second term
    I2 = zeros(Nm,Ndof);
    elem1 = [v1(:), v2(:), v1(:)+Nv];
    for im = 1:Nm
        gradmc = Gradmc{im};
        F1 = 1/6*sum(gradmc(x(v1), y(v1)).*Ne, 2);
        F2 = 1/6*sum(gradmc(x(v2), y(v2)).*Ne, 2);
        F3 = 4/6*sum(gradmc(xe, ye).*Ne, 2);
        F = [F1, F2, F3];
        I2(im, 1:2*Nv) = accumarray(elem1(:), F(:), [2*Nv 1]);
    end
    B1 = -I1 + I2;
    % constraint
    B1s = B1; B1s(1,1:Nv) = 1/Nv;
    G1 = B1*D;   G1s = B1s*D;

    % --------- sign matrix and sign vector ----------
    sgnedge = elem2sgn{iel}; sgnedge(sgnedge==0) = 1;
    sgnBase = ones(Ndof,1); sgnBase(2*Nv+1:3*Nv) = sgnedge;
    sgnK = sgnBase*sgnBase';  sgnF = sgnBase;

    % --------- stiffness matrix ----------
    % Projection
    Pi2s = G2s\B2s;   Pi2 = D*Pi2s;   I = eye(size(Pi2)); % biharmonic term
    Pi1s = G1s\B1s;   Pi1 = D*Pi1s;                       % Poisson term     
    % Stiffness matrix
    hh = diag([hxiK; he; he].^(-2));
    AK  = Pi2s'*G2*Pi2s + (hh*(I-Pi2))'*(I-Pi2);   % biharmonic
    BK  = Pi1s'*G1*Pi1s + (I-Pi1)'*(I-Pi1);  % Poisson     
    ABK = (para.epsilon^2*AK+BK).*sgnK;
    ABK = reshape(ABK',1,[]); % straighten

    % --------- load vector ----------
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y);  
    rhs = integralTri(fm,5,nodeT,elemT);
    fK = Pi0s'*rhs(:);
    fK = fK.*sgnF;

    % --------- assembly index ----------
    % matrix
    indexDof = [index, indexEdge+N, indexEdge+N+NE];  
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = ABK(:);
    ia = ia + Ndof^2;
    % vector
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    sgnPis = repmat(sgnBase',size(Pi2s,1),1);
    Ph{iel} = sgnPis.*Pi2s; 
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,NNdof,NNdof);
ff = accumarray(elemb,Fb,[NNdof 1]);

%% Apply Dirichlet boundary conditions
% boundary information
bdEdgeIdx = bdStruct.bdEdgeIdx; 
bdNodeIdx = bdStruct.bdNodeIdx; bdEdgeD = bdStruct.bdEdgeD;
g_D = pde.g_D;  Dw = pde.Du;
id = [bdNodeIdx; bdEdgeIdx+N; bdEdgeIdx+N+NE];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
pD = node(bdNodeIdx,:); wD = g_D(pD);
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
% if N < 2e3, solver = 'direct'; end
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
info.kk = kk; info.DofI = freeDof;