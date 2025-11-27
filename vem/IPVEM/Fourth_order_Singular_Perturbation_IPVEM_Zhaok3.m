function [uh,info] = Fourth_order_Singular_Perturbation_IPVEM_Zhaok3(node,elem,pde,bdStruct)

para = pde.para;

%% Preparations
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge;   
node2elem = auxT.node2elem;
elem2edge = auxT.elem2edge;
edge2elem = auxT.edge2elem;
edge2elemLocalIndex = auxT.edge2elemLocalIndex;
elem2sgn = auxT.elem2sgn;
% number
elemLen = cellfun('length',elem); % length of each elem
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+2*NE+3*NT;   Nm = 10;
% characteristic lengths of derivatives at nodes
hxi = cellfun(@(id) mean(diameter(id)), node2elem);
% Gauss-Lobatto weights and points: k = 3
ng = 4;
r = [-1, -1/sqrt(5), 1/sqrt(5), 1]; % [-1,1]
r = (r+1)/2; % [0,1]: gives the ratios of interior nodes
w = [1/12, 5/12, 5/12, 1/12];
cr = -1/2+r;

%% Interior biliner form
% store the projections
Ph = cell(NT,2); % H1,H2
% store the connectivity list
elem2dof = cell(NT,1);
nnz = sum((3*elemLen+3).^2);
ii = zeros(nnz,1);   jj = zeros(nnz,1);
ssA = zeros(nnz,1);  ssB = zeros(nnz,1);
nnz = sum(3*elemLen+3);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
% projection of basis functions at quadrature points
[Base,Basex,Basey,Basexx,Basexy,Baseyy] = deal(cell(NT,1));
for iel = 1:NT
    % --------- element information ---------
    index = elem{iel};  indexEdge = elem2edge{iel};  Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);  z = [x,y]; % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    hxiK = hxi(index); % characteristic length
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    za = z(v1,:)+r(2)*(z(v2,:)-z(v1,:)); % Gauss-Lobatto
    zb = z(v1,:)+r(3)*(z(v2,:)-z(v1,:));

    % ------- scaled monomials ---------
    m1 = @(x,y) 1+0*x;                     gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)/hK;                 gradm2 = @(x,y) [1/hK+0*x, 0+0*x];
    m3 = @(x,y) (y-yK)/hK;                 gradm3 = @(x,y) [0+0*x, 1/hK+0*x];
    m4 = @(x,y) (x-xK).^2/hK^2;            gradm4 = @(x,y) [2*(x-xK)/hK^2, 0+0*x];
    m5 = @(x,y) (x-xK).*(y-yK)/hK^2;       gradm5 = @(x,y) [(y-yK)/hK^2,  (x-xK)/hK^2];
    m6 = @(x,y) (y-yK).^2/hK^2;            gradm6 = @(x,y) [0+0*x, 2*(y-yK)/hK^2];
    m7 = @(x,y) (x-xK).^3/hK^3;            gradm7 = @(x,y) [3*(x-xK).^2/hK^3, 0+0*x];
    m8 = @(x,y) (x-xK).^2.*(y-yK)/hK^3;    gradm8 = @(x,y) [2*(x-xK).*(y-yK)/hK^3, (x-xK).^2/hK^3];
    m9 = @(x,y) (x-xK).*(y-yK).^2/hK^3;    gradm9 = @(x,y) [(y-yK).^2/hK^3, 2*(x-xK).*(y-yK)/hK^3];
    m10 = @(x,y) (y-yK).^3/hK^3;           gradm10 = @(x,y) [0+0*x, 3*(y-yK).^2/hK^3];
    mc = {m1,m2,m3,m4,m5,m6,m7,m8,m9,m10};
    Gradmc = {gradm1,gradm2,gradm3,gradm4,gradm5,gradm6,gradm7,...
        gradm8,gradm9,gradm10};   
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y),...
        m6(x,y), m7(x,y), m8(x,y), m9(x,y), m10(x,y)];
    mx = @(x,y) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x,...
        3*(x-xK).^2/hK^3,2*(x-xK).*(y-yK)/hK^3,(y-yK).^2/hK^3,0*x];
    my = @(x,y) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2,...
        0*x,(x-xK).^2/hK^3,2*(x-xK).*(y-yK)/hK^3,3*(y-yK).^2/hK^3];
    Gradm = @(x,y) [mx(x,y); my(x,y)]';
    mxx = @(x,y) [0*x, 0*x, 0*x, 2/hK^2+0*x, 0*x, 0*x,...
        6*(x-xK)/hK^3, 2*(y-yK)/hK^3,0*x,0*x];
    mxy = @(x,y) [0*x, 0*x, 0*x, 0*x, 1/hK^2+0*x, 0*x,...
        0*x,2*(x-xK)/hK^3,2*(y-yK)/hK^3,0*x];
    myy = @(x,y) [0*x, 0*x, 0*x, 0*x, 0*x, 2/hK^2+0*x,...
        0*x,0*x,2*(x-xK)/hK^3,6*(y-yK)/hK^3];

    % ------------ H1-elliptic projection on the lifting space ------------
    % dof number
    NdofL = 8*Nv+10; % lifting
    Ndof = 3*Nv+3;  % Poisson or IPVEM
    % transition matrix （lifting)
    DL = zeros(NdofL,Nm);
    DL(1:Nv,:) = m(x,y);
    DL(Nv+1:2*Nv,:) = repmat(hxiK,1,Nm).*mx(x,y);
    DL(2*Nv+1:3*Nv,:) = repmat(hxiK,1,Nm).*my(x,y);
    DL(3*Nv+1:4*Nv,:) = m(za(:,1), za(:,2)); % at ai
    DL(4*Nv+1:5*Nv,:) = m(zb(:,1), zb(:,2)); % at bi
    for i = 1:Nv  % moments of normal derivatives on edges
        g1 = Gradm(x(v1(i)),y(v1(i)));  
        ga = Gradm(za(i,1), za(i,2));   
        gb = Gradm(zb(i,1), zb(i,2));   
        g2 = Gradm(x(v2(i)),y(v2(i)));  
        DL(5*Nv+i,:) = Ne(i,:)*(w(1)*g1+w(2)*ga+w(3)*gb+w(4)*g2)';
        DL(6*Nv+i,:) = Ne(i,:)*(w(1)*cr(1)*g1+w(2)*cr(2)*ga+w(3)*cr(3)*gb+w(4)*cr(4)*g2)';
        DL(7*Nv+i,:) = Ne(i,:)*(w(1)*cr(1)^2*g1+w(2)*cr(2)^2*ga+w(3)*cr(3)^2*gb+w(4)*cr(4)^2*g2)';
    end
    for j = 1:Nm
        mjv = @(x,y) repmat(mc{j}(x,y),1,Nm).*m(x,y);
        DL(8*Nv+j,:) = 1/area(iel)*integralTri(mjv,5,nodeT,elemT);
    end
    % B, Bs, G, Gs
    idof = [1:Nv, 3*Nv+1:5*Nv, 8*Nv+(1:3)];
    I1 = zeros(Nm,Ndof);
    I1([4,6],3*Nv+1) = 2*area(iel)/hK^2;
    I1(7,3*Nv+2) = 6*area(iel)/hK^2;
    I1(8,3*Nv+3) = 2*area(iel)/hK^2;
    I1(9,3*Nv+2) = 2*area(iel)/hK^2;
    I1(10,3*Nv+3) = 6*area(iel)/hK^2;
    I2 = zeros(Nm,Ndof);
    elem1 = [v1(:), v1(:)+Nv, v1(:)+2*Nv, v2(:)];
    for im = 1:Nm
        Gradma = Gradmc{im};
        F1 = w(1)*sum(Gradma(x(v1),y(v1)).*Ne,2);
        F2 = w(2)*sum(Gradma(za(:,1), za(:,2)).*Ne,2);
        F3 = w(3)*sum(Gradma(zb(:,1), zb(:,2)).*Ne,2);
        F4 = w(4)*sum(Gradma(x(v2),y(v2)).*Ne,2);
        I2(im,1:3*Nv) = accumarray(elem1(:), [F1; F2; F3; F4], [3*Nv, 1]);
    end
    BL = zeros(Nm,NdofL);  BL(:,idof) = -I1 + I2;
    BsL = BL;
    BsL(1,1:Nv) = 1;
    GL = BL*DL; GsL = BsL*DL;
    % H1 elliptic
    PisL = GsL\BsL;   PiL = DL*PisL;
    % 
    % % ------------------- H2-elliptic projection ---------------------
    % transition matrix
    D = DL(idof,:);
    % use code for plate bending problems
    % Mij(m')
    M11 = @(x,y) -mxx(x,y);
    M12 = @(x,y) -mxy(x,y);
    M22 = @(x,y) -myy(x,y);
    % M_{ij,k}
    M111 = zeros(Nm,1);  M111(7) = -6/hK^3;
    M112 = zeros(Nm,1);  M112(8) = -2/hK^3;
    M121 = zeros(Nm,1);  M121(8) = -2/hK^3;
    M122 = zeros(Nm,1);  M122(9) = -2/hK^3;
    M221 = zeros(Nm,1);  M221(9) = -2/hK^3;
    M222 = zeros(Nm,1);  M222(10) = -6/hK^3;
    % H1-projection on Vk(K)
    Pi1s = PisL(:,idof);
    % B, Bs, G, Gs
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    [J0,J1,J2] = deal(zeros(Nm,Ndof));
    for i = 1:Nv % loop of edges or vertices
        % t,n
        t1 = te(i,1);  t2 = te(i,2);
        n1 = ne(i,1);  n2 = ne(i,2);
        % J0
        Q3t = (M111+M122)*n1+(M121+M222)*n2 ...
            + (M111*t1+M112*t2)*t1*n1 + (M121*t1+M122*t2)*t1*n2 ...
            + (M121*t1+M122*t2)*t2*n1 + (M221*t1+M222*t2)*t2*n2;
        phi1 = zeros(1,Ndof); phi1(i) = 1;  % at zi
        phia = zeros(1,Ndof); phia(Nv+i) = 1; % za
        phib = zeros(1,Ndof); phib(2*Nv+i) = 1; % zb
        phi2 = zeros(1,Ndof); phi2(v2(i)) = 1; % z_{i+1}
        J0 = J0 + he(i)*Q3t*(w(1)*phi1+w(2)*phia+w(3)*phib+w(4)*phi2);        
        % Gauss-Lobatto points on ei
        xg = [x(v1(i)), za(i,1), zb(i,1), x(v2(i))];
        yg = [y(v1(i)), za(i,2), zb(i,2), y(v2(i))];
        % J1
        for ig = 1:ng  % loop of G-L points on ei
            Mnn = M11(xg(ig),yg(ig))*n1*n1 + 2*M12(xg(ig),yg(ig))*n1*n2 ...
                + M22(xg(ig),yg(ig))*n2*n2; % (1,Nm)
            nphi = Ne(i,:)*(Gradm(xg(ig),yg(ig)))'*Pi1s;  % (1,Ndof), scaled
            J1 = J1 + w(ig)*Mnn'*nphi; % (Nm,Ndof)
        end
        % J2
        tn11 = te(p2(i),1)*ne(p2(i),1) - te(p1(i),1)*ne(p1(i),1); % jump
        tn22 = te(p2(i),2)*ne(p2(i),2) - te(p1(i),2)*ne(p1(i),2);
        tn12 = (te(p2(i),1)*ne(p2(i),2) + te(p2(i),2)*ne(p2(i),1)) ... % 12 and 21
            - (te(p1(i),1)*ne(p1(i),2) + te(p1(i),2)*ne(p1(i),1));
        Jump = M11(x(i),y(i))*tn11 + M12(x(i),y(i))*tn12 ...
            + M22(x(i),y(i))*tn22; % (1,Nm)
        phi = zeros(1,Ndof);  phi(i) = 1;  % at zi
        J2 = J2 + Jump'*phi; % (Nm,Ndof)
    end
    B = J0 -J1 + J2;
    Bs = B;
    % first constraint
    Bs(1,1:Nv) = 1;
    % second constraint
    Bs(2:3,1:Nv) = te(p1,:)' - te(p2,:)';
    for i = 1:Nv % loop of edges
        % Gauss-Lobatto points on ei
        xg = [x(v1(i)), za(i,1), zb(i,1), x(v2(i))];
        yg = [y(v1(i)), za(i,2), zb(i,2), y(v2(i))];
        for ig = 1:ng  % loop of G-L points
            nphi = Ne(i,:)*(Gradm(xg(ig),yg(ig)))'*Pi1s;  % (1,Ndof)
            Bs(2:3,:) = Bs(2:3,:) + w(ig)*ne(i,:)'*nphi;
        end
    end
    % consistency relation
    G = B*D;     Gs = Bs*D;

    % ------------------- L2 projection ---------------------
    C = zeros(Nm,Ndof);
    C(1,3*Nv+1) = area(iel);
    C(2,3*Nv+2) = area(iel);
    C(3,3*Nv+3) = area(iel);
    C(4:end,:) = area(iel)*PiL(8*Nv+4:end,idof);
    H = C*D;
    Pi0s = H\C;

    % ------------- H2 interior stiffness matrix ------------
    hh = diag([hxiK; he; he; hK; hK; hK].^(-2));
    Pis = Gs\Bs;   Pi = D*Pis;   I = eye(Ndof);
    kA  = Pis'*G*Pis + (hh*(I-Pi))'*(I-Pi);

    % ------------- H1 interior stiffness matrix ------------
    Pi1 = PiL(idof,idof);
    kB  = Pi1s'*GL*Pi1s + (I-Pi1)'*(I-Pi1); % G1 = GL = (\nabla m, \nabla m')
     
    % --------------- load vector ----------------------
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,5,nodeT,elemT);
    fK = Pi0s'*rhs(:);

    % ---------local to global index -----------
    % sgnBase
    sgnBase = elem2sgn{iel};
    sgnBase(sgnBase==0) = 1; % set as positive on the domain boundary
    sgnBase(sgnBase<0) = 0;
    elem1a = indexEdge+N*sgnBase+(N+NE)*(~sgnBase);
    elem1b = indexEdge+(N+NE)*sgnBase+N*(~sgnBase);
    indexDof = [index, elem1a, elem1b,...
        iel+N+2*NE, iel+N+2*NE+NT, iel+N+2*NE+2*NT];

    % --------- assembly index for interior stiffness matrix -----------
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
    xq = zeros(ng*Nv,1);     yq = zeros(ng*Nv,1);   % ng = 4
    xq(1:ng:end) = x(v1);    yq(1:ng:end) = y(v1);
    xq(2:ng:end) = za(:,1);  yq(2:ng:end) = za(:,2);
    xq(3:ng:end) = zb(:,1);  yq(3:ng:end) = zb(:,2);
    xq(4:ng:end) = x(v2);    yq(4:ng:end) = y(v2);

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
bdEdgeIdx = bdStruct.bdEdgeIdx;
bdEdge = bdStruct.bdEdge;
edge(bdEdgeIdx,:) = bdEdge; % adjustment on the boundary
v12 = node(edge(:,1),:)-node(edge(:,2),:);
Ne = [-v12(:,2),v12(:,1)];
he = vecnorm(v12,2,2);
ne = Ne./he;
% penalty parameter
if isfield(pde,'a')
    a = pde.a;
else
    a = 10;
end
if isfield(pde,'cp')
    cp = pde.cp*ones(NE,1);
else
    ce = 1/4*ones(NE,1);  ce(bdEdgeIdx) = 1;
    k1 = edge2elem(:,1);  k2 = edge2elem(:,2);
    Nv1 = elemLen(k1);    Nv2 = elemLen(k2);
    NK = max([Nv1,Nv2],2);
    T1 = area(k1)./Nv1;   T2 = area(k2)./Nv2;
    k = 3;
    cp = a*k*(k-1)*ce.*NK.*he.^2.*(1./T1+1./T2);
end
% adjustment for average: for domain boundary {u} = u = 2*(u+0)/2
Cav = ones(NE,1);  Cav(bdEdgeIdx) = 2;
% number of macro basis functions
edgeDof = sum(3*elemLen(edge2elem)+3,2)-ng;  % 3Nv+3
edgeDof(bdEdgeIdx) = (edgeDof(bdEdgeIdx) + ng)/2;
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
    sgnk1 = elem2sgn{k1};  sgne1 = sgnk1(e1);
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

 % -------- Local stiffness matrices of Jump terms -------------
    KJ = zeros(Ndofe,Ndofe);
    nx = ne(ie,1);  ny = ne(ie,2);
    for i = 1:Ndofe
        for j = 1:Ndofe
            % J2,J3
            vi = Jx(i,:)*nx + Jy(i,:)*ny;
            uj = Axx(j,:)*nx*nx + 2*Axy(j,:)*nx*ny + Ayy(j,:)*ny*ny;
            KJ(i,j) = he(ie)*sum(w.*vi.*uj);
        end
    end

    % -------- Local stiffness matrices of Penalty terms -------------
    KC = zeros(Ndofe,Ndofe);
    % expansion coefficients
    mm = zeros(5,5);
    for i = 0:4
        for j = 0:4
            cij = i+j+1;
            mm(i+1,j+1) = he(ie)/cij*(1/2^cij - 1/(-2)^cij);
        end
    end
    Aa = zeros(5,5);
    Aa(1,:) = (-1/2).^(0:4);
    Aa(2,:) = (1/2).^(0:4);
    Aa(3:5,:) = mm(1:3,:);
    a = zeros(5,Ndofe);  % (k+1)-th order polynomials: a0,...,a_{k+1}
    cr = r-1/2;
    for i = 1:Ndofe
        b = zeros(5,1);
        vi = Jx(i,:)*nx + Jy(i,:)*ny;
        b(1) = vi(1);  % zi
        b(2) = vi(ng); % z_{i+1}
        b(3) = he(ie)*dot(w,vi);
        b(4) = he(ie)*dot(w,cr.*vi);
        b(5) = he(ie)*dot(w,cr.^2.*vi);
        a(:,i) = Aa\b;
    end
    for i = 1:Ndofe
        for j = 1:Ndofe
            KC(i,j) = cp(ie)/he(ie)*a(:,j)'*mm*a(:,i);
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

%% Apply "Neumann" boundary conditions
% boundary information
NEb = size(bdEdge,1);  
idElemb = edge2elem(bdEdgeIdx,1);
z1 = node(bdEdge(:,1),:); 
z2 = node(bdEdge(:,2),:); 
za = z1+r(2)*(z2-z1); % Gauss-Lobatto
zb = z1+r(3)*(z2-z1);
% loop of boundary edges
nnz = sum(3*elemLen(idElemb)+3);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ib = 0;
for s = 1:NEb
    % element information    
    ie = bdEdgeIdx(s); % index of boundary edge   
    iel = idElemb(s); % index of the element    
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);      
    Dofe = elem2dof{iel};
    ndofe = length(Dofe);
    Pi1s = Ph{iel,1};  % H1 projection on e

    % scaled monomials
    mx = @(x,y) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x,...
        3*(x-xK).^2/hK^3,2*(x-xK).*(y-yK)/hK^3,(y-yK).^2/hK^3,0*x];
    my = @(x,y) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2,...
        0*x,(x-xK).^2/hK^3,2*(x-xK).*(y-yK)/hK^3,3*(y-yK).^2/hK^3];
    mxx = @(x,y) [0*x, 0*x, 0*x, 2/hK^2+0*x, 0*x, 0*x,...
        6*(x-xK)/hK^3, 2*(y-yK)/hK^3,0*x,0*x];
    mxy = @(x,y) [0*x, 0*x, 0*x, 0*x, 1/hK^2+0*x, 0*x,...
        0*x,2*(x-xK)/hK^3,2*(y-yK)/hK^3,0*x];
    myy = @(x,y) [0*x, 0*x, 0*x, 0*x, 0*x, 2/hK^2+0*x,...
        0*x,0*x,2*(x-xK)/hK^3,6*(y-yK)/hK^3];

    % basis matrix
    nx = ne(ie,1);  ny = ne(ie,2);  
    xq = [z1(s,1);za(s,1);zb(s,1);z2(s,1)];  % Gauss-Lobatto
    yq = [z1(s,2);za(s,2);zb(s,2);z2(s,2)];  
    zq = [xq, yq];
    phix = (mx(xq,yq)*Pi1s)';  % (ndofe,ng)
    phiy = (my(xq,yq)*Pi1s)'; 
    phixx = (mxx(xq,yq)*Pi1s)'; 
    phixy = (mxy(xq,yq)*Pi1s)'; 
    phiyy = (myy(xq,yq)*Pi1s)'; 

     % expansion coefficients for jump term
    mm = zeros(5,5); % (k+1)-th order polynomials: a0,...,a_{k+1}
    for i = 0:4
        for j = i:4
            cij = i+j+1;
            mm(i+1,j+1) = he(ie)/cij*(1/2^cij - 1/(-2)^cij);
            mm(j+1,i+1) = mm(i+1,j+1);
        end
    end
    Aa = zeros(5,5);    
    Aa(1,:) = cr(1).^(0:4);
    Aa(2,:) = cr(end).^(0:4);
    Aa(3:5,:) = mm(1:3,:);
    a = zeros(5,ndofe);      
    for i = 1:ndofe
        vi = phix(i,:)*nx + phiy(i,:)*ny;  
        b = zeros(4,1);        
        b(1) = vi(1);  % zi
        b(2) = vi(ng); % z_{i+1}
        b(3) = he(ie)*dot(w,vi);
        b(4) = he(ie)*dot(w,cr.*vi);
        b(5) = he(ie)*dot(w,cr.^2.*vi);
        a(:,i) = Aa\b;
    end

    % jump and average    
    ave = phixx*nx*nx + 2*phixy*nx*ny + phiyy*ny*ny; % average: (ndofe,ng)
    jump = zeros(ndofe,ng);
    for i = 1:ndofe
        for p = 0:4 % 0:k+1
            jump(i,:) = jump(i,:) + a(p+1,i)*cr.^p;  % a_p*m_p^e on Gauss-Lobatto
        end
    end
    ajump = -ave + cp(ie)/he(ie)*jump;  % (ndofe,ng)

    % local load vector
    gN = pde.Du(zq);  
    wgN = w(:).*sum(repmat(Ne(ie,:),ng,1).*gN,2);  % (ng,1)
    fe = ajump*wgN; 

    % assembly index for right hand side 
    ndofe = length(Dofe);
    elemb(ib+1:ib+ndofe) = Dofe(:);
    Fb(ib+1:ib+ndofe) = fe(:);
    ib = ib + ndofe;   
end
ff = ff + para.epsilon^2*accumarray(elemb,Fb,[NNdof 1]);

%% Apply Dirichlet boundary conditions
% boundary information
g_D = pde.g_D;
bdEdge = bdStruct.bdEdge;
bdNodeIdx = bdStruct.bdNodeIdx;
bdEdgeIdxa = bdEdgeIdx + N;
bdEdgeIdxb = bdEdgeIdx + N+NE;
id = [bdNodeIdx; bdEdgeIdxa; bdEdgeIdxb];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
pD = node(bdNodeIdx,:); uD = g_D(pD);
% ai, bi on the boundary
z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:);
za = z1+r(2)*(z2-z1); zb = z1+r(3)*(z2-z1); % Gauss-Lobatto
uDa = g_D(za); uDb = g_D(zb);
% rhs
uh = zeros(NNdof,1); uh(bdDof) = [uD; uDa; uDb];
ff = ff - kk*uh;

%% Set solver
uh(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm