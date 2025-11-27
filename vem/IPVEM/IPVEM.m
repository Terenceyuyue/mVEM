function [uh,info] = IPVEM(node,elem,pde,bdStruct)
%IPVEM solves biharmonic equation using the interior penalty virtual
%element method in the lowest-order case:
%
%     \Delta^2 u = f,  in Omega
%     Dirichlet boundary condition u=g_D on \Gamma_D, 
%     Neumann boundary condition   grad(u)*n=g_N on \Gamma_N
%
% Note: - H1-elliptic projector is included in J1(v,w), J2(v,w), J3(v,w)
%       - The handling of inhomogeneous boundary conditions is included
%
%   References
%   [1] J. Zhao, S. Mao, B. Zhang and F. Wang, "The interior penalty virtual element method 
%   for the biharmonic problem", Math. Comp., Vol 92. No 342., pp. 1543-1574, 2023.
%   [2] F. Feng and Y. Yu, "A modified interior penalty virtual element method for fourth-order
%   singular perturbation problems", J. Sci. Comput., Vol 101., p. 21, 2024.
%
% Copyright (C)  Terence Yu. 

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
NNdof = N+NE+NT;   Nm = 6;
% characteristic lengths of derivatives at nodes
hxi = cellfun(@(id) mean(diameter(id)), node2elem);
% Gauss-Lobatto: Simpson for k = 2
ng = 3;
r = [-1, 0, 1]; % [-1,1]
r = (r+1)/2; % [0,1]: gives the ratios
w = [1/6, 4/6, 1/6];
cr = -1/2+r;

%% Interior biliner form
% store the projections
Ph = cell(NT,2); % H1,H2
Ph0 = cell(NT,1); % L2
[Dh,Hh] = deal(cell(NT,1));
% store the connectivity list
elem2dof = cell(NT,1);
nnz = sum((2*elemLen+1).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
nnz = sum(2*elemLen+1);
elemf = zeros(nnz,1); Ff = zeros(nnz,1);
ia = 0; ib = 0;
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
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];

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
        g1 = Gradm(x(v1(i)),y(v1(i)));   
        ge = Gradm(xe(i),ye(i));         
        g2 = Gradm(x(v2(i)),y(v2(i)));   
        DL(4*Nv+i,:) = Ne(i,:)*(w(1)*g1+w(2)*ge+w(3)*g2)';
        DL(5*Nv+i,:) = Ne(i,:)*(w(1)*cr(1)*g1+w(2)*cr(2)*ge+w(3)*cr(3)*g2)';
    end
    for j = 1:Nm
        mjv = @(x,y) repmat(mc{j}(x,y),1,Nm).*m(x,y);
        DL(6*Nv+j,:) = 1/area(iel)*integralTri(mjv,5,nodeT,elemT);
    end
    % B, Bs, G, Gs
    idof = [1:Nv, 3*Nv+1:4*Nv, 6*Nv+1]; % index for Poisson
    I1 = zeros(Nm,Ndof);
    I1(:,end) = area(iel)*Lapm';
    I2 = zeros(Nm,Ndof);
    elem1 = [v1(:), v1(:)+Nv, v2(:)];
    for im = 1:Nm
        Gradma = Gradmc{im};
        F1 = w(1)*sum(Gradma(x(v1),y(v1)).*Ne,2);
        F2 = w(2)*sum(Gradma(xe,ye).*Ne,2);
        F3 = w(3)*sum(Gradma(x(v2),y(v2)).*Ne,2);
        I2(im,:) = accumarray(elem1(:), [F1; F2; F3], [Ndof, 1]);
    end
    BL = zeros(Nm,NdofL);  BL(:,idof) = -I1 + I2;
    BsL = BL;
    BsL(1,1:Nv) = 1;
    GsL = BsL*DL;
    % H1 elliptic
    PisL = GsL\BsL;   PiL = DL*PisL;

    % ------------------- H2-elliptic projection ---------------------
    % transition matrix
    D = DL(idof,:);
    Dh{iel} = D;
    % use code for plate bending problems
    % Mij(m')
    M11 = @(x,y) -mxx(x,y);
    M12 = @(x,y) -mxy(x,y);
    M22 = @(x,y) -myy(x,y);
    % H1-projection on Vk(K)
    Pi1s = PisL(:,idof);
    % B, Bs, G, Gs
    Ndof = 2*Nv+1;
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    J1 = zeros(Nm,Ndof); 
    J2 = zeros(Nm,Ndof);
    for i = 1:Nv % loop of edges or vertices
        % n
        n1 = ne(i,1);  n2 = ne(i,2);
        % Gauss-Lobatto points on ei
        xg = [x(v1(i)), xe(i), x(v2(i))];
        yg = [y(v1(i)), ye(i), y(v2(i))];
        % J1
        for ig = 1:ng
            Mnn = M11(xg(ig),yg(ig))*n1*n1 + 2*M12(xg(ig),yg(ig))*n1*n2 ...
                + M22(xg(ig),yg(ig))*n2*n2; % (1,Nm)
            nphi = Ne(i,:)*(Gradm(xg(ig),yg(ig)))'*Pi1s;  % (1,Ndof), scaled
            J1 = J1 + w(ig)*Mnn'*nphi; % (Nm,Ndof)
        end
        % J2
        tn11 = te(p2(i),1)*ne(p2(i),1) - te(p1(i),1)*ne(p1(i),1); % jump
        tn22 = te(p2(i),2)*ne(p2(i),2) - te(p1(i),2)*ne(p1(i),2);
        tn12 = (te(p2(i),1)*ne(p2(i),2) + te(p2(i),2)*ne(p2(i),1)) ...
            - (te(p1(i),1)*ne(p1(i),2) + te(p1(i),2)*ne(p1(i),1));
        Jump = M11(x(i),y(i))*tn11 + M12(x(i),y(i))*tn12 ...
            + M22(x(i),y(i))*tn22; % (1,Nm)
        % phi' at zi
        phi = zeros(1,Ndof);  phi(i) = 1;
        % B1 on e and at zi
        J2 = J2 + Jump'*phi;
    end
    B = -J1 + J2;
    Bs = B;
    % first constraint
    Bs(1,1:Nv) = 1;
    % second constraint
    Bs(2:3,1:Nv) = te(p1,:)' - te(p2,:)';
    for i = 1:Nv % loop of edges
        % Gauss-Lobatto points on ei
        xg = [x(v1(i)), xe(i), x(v2(i))];
        yg = [y(v1(i)), ye(i), y(v2(i))];
        for ig = 1:ng
            nphi = Ne(i,:)*(Gradm(xg(ig),yg(ig)))'*Pi1s;  % (1,Ndof), scaled
            Bs(2:3,:) = Bs(2:3,:) + w(ig)*ne(i,:)'*nphi;
        end        
    end
    % consistency relation
    G = B*D;     Gs = Bs*D;

    % ------------------- L2 projection ---------------------
    C = zeros(Nm,Ndof);
    C(1,2*Nv+1) = area(iel);
    C(2:end,:) = area(iel)*PiL(6*Nv+2:end,idof);
    H = C*D;
    Pi0s = H\C;
    Hh{iel} = H;

    % ------------- interior stiffness matrix ------------
    hh = diag([hxiK; he; hK].^(-2));
    Pis = Gs\Bs;   Pi = D*Pis;   I = eye(Ndof);
    kA  = Pis'*G*Pis + (hh*(I-Pi))'*(I-Pi);

    % --------------- load vector ----------------------
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,5,nodeT,elemT);
    fK = Pi0s'*rhs(:);

    % --------- assembly index for interior stiffness matrix -----------
    indexDof = [index, indexEdge+N, iel+N+NE];
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = reshape(kA',1,[]);
    ia = ia + Ndof^2;

    % --------- assembly index for right hand side -----------
    elemf(ib+1:ib+Ndof) = indexDof(:);
    Ff(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;

    % --------- matrix for error evaluation  ---------
    Ph{iel,1} = Pi1s; % H1
    Ph{iel,2} = Pis; % H2
    elem2dof{iel} = indexDof;

    % ---- Elementwise evaluations of H1 projections of basis functions ---
    % quadrature points: Nv edges
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
A = sparse(ii,jj,ss,NNdof,NNdof);
ff = accumarray(elemf,Ff,[NNdof 1]);

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
    a = 1;
end
if isfield(pde,'cp')
    cp = pde.cp*ones(NE,1);
else
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
    [Px,Py,Pxx,Pxy,Pyy] = deal(zeros(Ndofe, 2*ng)); % initialized as zero matrix
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
            % J2, J3 = J2'
            vi = Jx(i,:)*nx + Jy(i,:)*ny;
            uj = Axx(j,:)*nx*nx + 2*Axy(j,:)*nx*ny + Ayy(j,:)*ny*ny;
            KJ(i,j) = he(ie)*sum(w.*vi.*uj);
            % J1
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
J = sparse(ii,jj,ssJ,NNdof,NNdof); % J2, J3 = J2'
C = sparse(ii,jj,ssC,NNdof,NNdof); % J1
kk = A - J - J' + C;

%% Apply "Neumann" boundary conditions
% boundary information
NEb = size(bdEdge,1);  
idElemb = edge2elem(bdEdgeIdx,1);
z1 = node(bdEdge(:,1),:); 
z2 = node(bdEdge(:,2),:); 
zc = (z1+z2)/2;
% loop of boundary edges
nnz = sum(2*elemLen(idElemb)+1);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ib = 0;
for s = 1:NEb
    % element information    
    ie = bdEdgeIdx(s); % index of boundary edge   
    iel = idElemb(s); % index of the element
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);      
    Dofe = elem2dof{iel};
    Pi1s = Ph{iel,1};  % H1 projection on e

    % scaled monomials
    mx = @(x,y) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x];
    my = @(x,y) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2];
    mxx = @(x,y) [0*x, 0*x, 0*x, 2/hK^2+0*x, 0*x, 0*x];
    mxy = @(x,y) [0*x, 0*x, 0*x, 0*x, 1/hK^2+0*x, 0*x];
    myy = @(x,y) [0*x, 0*x, 0*x, 0*x, 0*x, 2/hK^2+0*x];
	
	% base matrix
	nx = ne(ie,1);  ny = ne(ie,2);  
    xq = [z1(s,1);zc(s,1);z2(s,1)];  % Gauss-Lobatto
    yq = [z1(s,2);zc(s,2);z2(s,2)]; 
    zq = [xq, yq];
    phix = (mx(xq,yq)*Pi1s)';  % (ndofe,ng)
    phiy = (my(xq,yq)*Pi1s)'; 
    phixx = (mxx(xq,yq)*Pi1s)'; 
    phixy = (mxy(xq,yq)*Pi1s)'; 
    phiyy = (myy(xq,yq)*Pi1s)'; 

    % jump and average    
    ave = phixx*nx*nx + 2*phixy*nx*ny + phiyy*ny*ny; % average 
    jump = phix*nx + phiy*ny; % jump
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
ff = ff + accumarray(elemb,Fb,[NNdof 1]);

%% Apply Dirichlet boundary conditions
% boundary information
g_D = pde.g_D;
bdNodeIdx = bdStruct.bdNodeIdx;
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
info.Ph = Ph;  info.Ph0 = Ph0;
info.elem2dof = elem2dof;
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm
info.Dh = Dh; info.Hh = Hh;  info.cp = cp;
info.diameter = diameter;