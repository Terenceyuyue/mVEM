function [uh,info] = Fourth_order_Singular_Perturbation_IPVEMTri(node,elem,pde,bdStruct)

para = pde.para;
para.D = 1;  para.nu = 0;

%% Preparations
% auxgeometry
isTri = true;
aux = auxgeometry(node,elem,isTri);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% triangulation
nodeTri = aux.nodeTri;
elemTri = aux.elemTri;
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge;   bdEdgeIdx = auxT.bdEdgeIdx;
node2elem = auxT.node2elem;
elem2edge = auxT.elem2edge;
edge2elem = auxT.edge2elem;
edge2elemLocalIndex = auxT.edge2elemLocalIndex;
elemLen = cellfun('length',elem); % length of each elem
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+NE+NT;   Nm = 6;
% characteristic lengths of derivatives at nodes
hxi = cellfun(@(id) mean(diameter(id)), node2elem);
% Gauss-Lobatto: Simpson for k = 2
ng = 3;
w = [1/6, 4/6, 1/6];

%% Interior biliner form
% store the projections
Ph = cell(NT,2); % H1,H2
% store the connectivity list
elem2dof = cell(NT,1);
nnz = sum((2*elemLen+1).^2);
ii = zeros(nnz,1);   jj = zeros(nnz,1);
ssA = zeros(nnz,1);  ssB = zeros(nnz,1);
nnz = sum(2*elemLen+1);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
% sgnelem
sgnelem = cell(NT,1);
E = false(NE,1); E(bdEdgeIdx) = 1;
% projection of basis functions at quadrature points
[Base,Basex,Basey,Basexx,Basexy,Baseyy] = deal(cell(NT,1));
for iel = 1:NT
    % --------- element information ---------
    index = elem{iel};  indexEdge = elem2edge{iel};  Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    xe = (x(v1)+x(v2))/2; ye = (y(v1)+y(v2))/2; % mid-edge
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    hxiK = hxi(index); % characteristic length
    nodeT = nodeTri{iel};
    elemT = elemTri{iel};
    % sgnelem
    sgnL = sign(index([2:Nv,1]) - index(1:Nv));
    id = elem2edge{iel}; sgnbd = E(id); sgnL(sgnbd) = 0;  % on the domain boundary
    sgnelem{iel} = sgnL;

    % ------- scaled monomials ---------
    m1 = @(x,y) 1+0*x;                     gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)/hK;                 gradm2 = @(x,y) [1/hK+0*x, 0+0*x];
    m3 = @(x,y) (y-yK)/hK;                 gradm3 = @(x,y) [0+0*x, 1/hK+0*x];
    m4 = @(x,y) (x-xK).^2/hK^2;            gradm4 = @(x,y) [2*(x-xK)/hK^2, 0+0*x];
    m5 = @(x,y) (x-xK).*(y-yK)/hK^2;       gradm5 = @(x,y) [(y-yK)/hK^2,  (x-xK)/hK^2];
    m6 = @(x,y) (y-yK).^2/hK^2;            gradm6 = @(x,y) [0+0*x, 2*(y-yK)/hK^2];
    mc = {m1,m2,m3,m4,m5,m6};
    Gradmc = {gradm1,gradm2,gradm3,gradm4,gradm5,gradm6};
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)];
    mx = @(x,y) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x];
    my = @(x,y) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2];
    Gradm = @(x,y) [mx(x,y); my(x,y)]';
    mxx = @(x,y) [0*x, 0*x, 0*x, 2/hK^2+0*x, 0*x, 0*x];
    mxy = @(x,y) [0*x, 0*x, 0*x, 0*x, 1/hK^2+0*x, 0*x];
    myy = @(x,y) [0*x, 0*x, 0*x, 0*x, 0*x, 2/hK^2+0*x];
    Lapm = [0,0,0,2/hK^2,0,2/hK^2];

    % ------------ H1-elliptic projection on the lifting space ------------
    % dof number
    NdofL = 6*Nv+6; % lifting
    Ndof = 2*Nv+1;  % Poisson or IPVEM
    % transition matrix （lifting)
    DL = zeros(NdofL,Nm);
    DL(1:Nv,:) = m(x,y);
    DL(Nv+1:2*Nv,:) = repmat(hxiK,1,Nm).*mx(x,y);
    DL(2*Nv+1:3*Nv,:) = repmat(hxiK,1,Nm).*my(x,y);
    DL(3*Nv+1:4*Nv,:) = m(xe,ye);
    for i = 1:Nv  % moments of normal derivatives on edges
        ga = Gradm(x(v1(i)),y(v1(i)));   ca = -1/2;  % (s-se)/he
        ge = Gradm(xe(i),ye(i));         ce = 0;
        gb = Gradm(x(v2(i)),y(v2(i)));   cb = 1/2;
        DL(4*Nv+i,:) = 1/2*Ne(i,:)*(ga+gb)';
        DL(5*Nv+i,:) = 1/6*Ne(i,:)*(ca*ga+4*ce*ge+cb*gb)';
    end
    for j = 1:6
        mjv = @(x,y) repmat(mc{j}(x,y),1,Nm).*m(x,y);
        DL(6*Nv+j,:) = 1/area(iel)*integralTri(mjv,5,nodeT,elemT);
    end
    % B, Bs, G, Gs
    idof = [1:Nv, 3*Nv+1:4*Nv, 6*Nv+1];
    I1 = zeros(Nm,Ndof);
    I1(:,end) = -area(iel)*Lapm';
    I2 = zeros(Nm,Ndof);
    elem1 = [v1(:), v1(:)+Nv, v2(:)];
    for im = 1:Nm
        Gradma = Gradmc{im};
        F1 = 1/6*sum(Gradma(x(v1),y(v1)).*Ne,2);
        F2 = 4/6*sum(Gradma(xe,ye).*Ne,2);
        F3 = 1/6*sum(Gradma(x(v2),y(v2)).*Ne,2);
        I2(im,:) = accumarray(elem1(:), [F1; F2; F3], [Ndof, 1]);
    end
    BL = zeros(Nm,NdofL);  BL(:,idof) = I1 + I2;
    BsL = BL;
    BsL(1,1:Nv) = 1;
    GL = BL*DL; GsL = BsL*DL;
    % H1 elliptic
    PisL = GsL\BsL;   PiL = DL*PisL;

    % ------------------- H2-elliptic projection ---------------------
    % transition matrix
    D = DL(idof,:);
    % use code for plate bending problems
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
    % H1-projection on Vk(K)
    Pi1s = PisL(:,idof);
    % B, Bs, G, Gs
    Ndof = 2*Nv+1;
    B = zeros(Nm,Ndof);
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    for i = 1:Nv % loop of edges or vertices
        % int[\partial_n (phi')] on ei
        Dma = Ne(i,:)*Gradm(x(v1(i)),y(v1(i)))';
        Dmb = Ne(i,:)*Gradm(x(v2(i)),y(v2(i)))';
        g = 1/2*(Dma+Dmb);
        nphi = g*Pi1s;
        % Jump(m) at zi
        Jump = Mtn(:,p2(i))-Mtn(:,p1(i));
        % phi' at zi
        phi = zeros(1,Ndof);  phi(i) = 1;
        % B1 on e and at zi
        B = B - Mnn(:,i)*nphi + Jump*phi;
    end
    Bs = B;
    % first constraint
    Bs(1,1:Nv) = 1;
    % second constraint
    Bs(2:3,1:Nv) = te(p1,:)' - te(p2,:)';
    for i = 1:Nv % loop of edges
        % int[\partial_n (phi')] on ei
        Dma = Ne(i,:)*Gradm(x(v1(i)),y(v1(i)))';
        Dmb = Ne(i,:)*Gradm(x(v2(i)),y(v2(i)))';
        g = 1/2*(Dma+Dmb);
        Nphi = ne(i,:)'*(g*Pi1s);
        Bs(2:3,:) = Bs(2:3,:) + Nphi;
    end
    % consistency relation
    G = B*D;     Gs = Bs*D;

    % ------------------- L2 projection ---------------------
    C = zeros(Nm,Ndof);
    C(1,2*Nv+1) = area(iel);
    C(2:end,:) = area(iel)*PiL(6*Nv+2:end,idof);
    H = C*D;
    Pi0s = H\C;

    % ------------- H2 interior stiffness matrix ------------
    Pis = Gs\Bs;   Pi = D*Pis;   I = eye(Ndof);
    kA  = Pis'*G*Pis + hK^(-2)*(I-Pi)'*(I-Pi);

    % ------------- H1 interior stiffness matrix ------------
    Pi1 = PiL(idof,idof);
    kB  = Pi1s'*GL*Pi1s + (I-Pi1)'*(I-Pi1); % G1 = GL = (\nabla m, \nabla m')

    % --------------- load vector ----------------------
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,5,nodeT,elemT);
    fK = Pi0s'*rhs(:);

    % --------- assembly index for interior stiffness matrix -----------
    indexDof = [index, indexEdge+N, iel+N+NE];
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ssA(ia+1:ia+Ndof^2) = reshape(kA',1,[]);
    ssB(ia+1:ia+Ndof^2) = reshape(kB',1,[]);
    ia = ia + Ndof^2;

    % --------- assembly index for right hand side -----------
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;

    % --------- matrix for error evaluation  ---------
    Ph{iel,1} = Pi1s; % H1
    Ph{iel,2} = Pis; % H2
    elem2dof{iel} = indexDof;

    % --- Elementwise evaluations of H1 projections of basis ---
    % quadrature points
    xq = zeros(ng*Nv,1);   yq = zeros(ng*Nv,1);   % ng = 3
    xq(1:3:end) = x(v1);   yq(1:3:end) = y(v1);
    xq(2:3:end) = xe;      yq(2:3:end) = ye;
    xq(3:3:end) = x(v2);   yq(3:3:end) = y(v2);

    % evaluations at quadrature points
    Base{iel} = (m(xq,yq)*Pi1s)';   % (Ndof, Nv*ng)
    Basex{iel} = (mx(xq,yq)*Pi1s)';
    Basey{iel} = (my(xq,yq)*Pi1s)';
    Basexx{iel} = (mxx(xq,yq)*Pi1s)';
    Basexy{iel} = (mxy(xq,yq)*Pi1s)';
    Baseyy{iel} = (myy(xq,yq)*Pi1s)';
end
A = sparse(ii,jj,ssA,NNdof,NNdof);
B = sparse(ii,jj,ssB,NNdof,NNdof);
ff = accumarray(elemb,Fb,[NNdof 1]);

%% Edgewise jumps and averages
% ne,he
v12 = node(edge(:,1),:)-node(edge(:,2),:);
Ne = [-v12(:,2),v12(:,1)];
he = vecnorm(v12,2,2);
ne = Ne./he;
% penalty parameter
if isfield(pde,'cp')
    cp = pde.cp*ones(NE,1);
else
    a = 10;
    ce = 1/4*ones(NE,1);  ce(bdEdgeIdx) = 1;
    k1 = edge2elem(:,1);  k2 = edge2elem(:,2);
    Nv1 = elemLen(k1);    Nv2 = elemLen(k2);
    NK = max([Nv1,Nv2],2);
    T1 = area(k1)./Nv1;   T2 = area(k2)./Nv2;
    k = 2;
    cp = a*k*(k-1)*ce.*NK.*he.^2.*(1./T1+1./T2);
end
% adjustment for average: for domain boundary {u} = u = 2*(u+0)/2
Cav = ones(NE,1);  Cav(bdEdgeIdx) = 2;
% number of macro basis functions
edgeDof = sum(2*elemLen(edge2elem)+1,2)-3;
edgeDof(bdEdgeIdx) = (edgeDof(bdEdgeIdx) + 3)/2;
nnz = sum(edgeDof.^2);
[ii,jj,ssJ,ssC] = deal(zeros(nnz,1));
% loop of edges
ia = 0;
for ie = 1:NE
    % --------- info of left and right elements -----------
    k1 = edge2elem(ie,1);  k2 = edge2elem(ie,2);
    indexDofk1 = elem2dof{k1};   Ndofk1 = length(indexDofk1);
    indexDofk2 = elem2dof{k2};   Ndofk2 = length(indexDofk2);
    % local index of e on the left or right element
    e1 = edge2elemLocalIndex(ie,1);  e2 = edge2elemLocalIndex(ie,2);
    % sign of the edge on the left element
    sgnk1 = sgnelem{k1};  sgne1 = sgnk1(e1);
    % macro indices of left and right basis functions
    [indexDofe,~, idMacro] = unique([indexDofk1,indexDofk2]);
    Ndofe = length(indexDofe);

    % ---------  macro base matrix ---------
    phix1 = Basex{k1};     phix2 = Basex{k2};
    phiy1 = Basey{k1};     phiy2 = Basey{k2};
    phixx1 = Basexx{k1};   phixx2 = Basexx{k2};
    phixy1 = Basexy{k1};   phixy2 = Basexy{k2};
    phiyy1 = Baseyy{k1};   phiyy2 = Baseyy{k2};
    [Px,Py,Pxx,Pxy,Pyy] = deal(zeros(Ndofe, 2*ng));
    % left
    for k = 1:Ndofk1
        row = idMacro(k); % r
        Px(row, 1:ng) = phix1(k, (1:ng)+(e1-1)*ng);
        Py(row, 1:ng) = phiy1(k, (1:ng)+(e1-1)*ng);
        Pxx(row, 1:ng) = phixx1(k, (1:ng)+(e1-1)*ng);
        Pxy(row, 1:ng) = phixy1(k, (1:ng)+(e1-1)*ng);
        Pyy(row, 1:ng) = phiyy1(k, (1:ng)+(e1-1)*ng);
    end
    % right: k1~=k2  % zero for the exterior of the boundary edges (k1=k2)
    if k1~=k2
        for k = 1:Ndofk2
            row = idMacro(k+Ndofk1); % s
            Px(row, (1:ng)+ng) = phix2(k, (1:ng)+(e2-1)*ng);
            Py(row, (1:ng)+ng) = phiy2(k, (1:ng)+(e2-1)*ng);
            Pxx(row, (1:ng)+ng) = phixx2(k, (1:ng)+(e2-1)*ng);
            Pxy(row, (1:ng)+ng) = phixy2(k, (1:ng)+(e2-1)*ng);
            Pyy(row, (1:ng)+ng) = phiyy2(k, (1:ng)+(e2-1)*ng);
        end
    end
    % reorder w.r.t the positive sign
    ReverseCols = [(1:ng), (ng:-1:1)+ng];
    if sgne1<0
        ReverseCols = [(1:ng)+ng, (ng:-1:1)];
    end
    Px = Px(:, ReverseCols);
    Py = Py(:, ReverseCols);
    Pxx = Pxx(:, ReverseCols);
    Pxy = Pxy(:, ReverseCols);
    Pyy = Pyy(:, ReverseCols);

    % ---------  jumps and averages of Projections ---------
    Jx = Px(:,1:ng) - Px(:,(1:ng)+ng);
    Jy = Py(:,1:ng) - Py(:,(1:ng)+ng);
    Axx = Cav(ie)*0.5*(Pxx(:,1:ng) + Pxx(:,(1:ng)+ng));
    Axy = Cav(ie)*0.5*(Pxy(:,1:ng) + Pxy(:,(1:ng)+ng));
    Ayy = Cav(ie)*0.5*(Pyy(:,1:ng) + Pyy(:,(1:ng)+ng));

    % -------- Local stiffness matrices of Jump and Penalty terms -------------
    [KJ,KC] = deal(zeros(Ndofe,Ndofe));
    nx = ne(ie,1);  ny = ne(ie,2);
    for i = 1:Ndofe
        for j = 1:Ndofe
            % J
            vi = Jx(i,:)*nx + Jy(i,:)*ny;
            uj = Axx(j,:)*nx*nx + 2*Axy(j,:)*nx*ny + Ayy(j,:)*ny*ny;
            KJ(i,j) = he(ie)*sum(w.*vi.*uj);
            % C
            vi = Jx(i,:)*nx + Jy(i,:)*ny;
            uj = Jx(j,:)*nx + Jy(j,:)*ny;
            KC(i,j) = cp(ie)*sum(w.*vi.*uj);
        end
    end

    % ---------  macro connectivity -------
    ii(ia+1:ia+Ndofe^2) = reshape(repmat(indexDofe, Ndofe, 1), [], 1);
    jj(ia+1:ia+Ndofe^2) = repmat(indexDofe(:), Ndofe, 1);
    ssJ(ia+1:ia+Ndofe^2) = reshape(KJ',1,[]);
    ssC(ia+1:ia+Ndofe^2) = reshape(KC',1,[]);
    ia = ia + Ndofe^2;
end

%% Assemble the jump and penalty terms
J = sparse(ii,jj,ssJ,NNdof,NNdof);
C = sparse(ii,jj,ssC,NNdof,NNdof);
kk = para.epsilon^2*(A - J - J' + C) + B;

%% Apply Dirichlet boundary conditions
% boundary information
g_D = pde.g_D;
bdEdge = bdStruct.bdEdge;
bdNodeIdx = bdStruct.bdNodeIdx;
z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:); zc = (z1+z2)./2;
id = [bdNodeIdx; bdEdgeIdx+N];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
pD = node(bdNodeIdx,:); uD = g_D(pD);
% mid-edge on the boundary
uDc = g_D(zc);
% rhs
uh = zeros(NNdof,1); uh(bdDof) = [uD; uDc];
ff = ff - kk*uh;

%% Set solver
uh(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm