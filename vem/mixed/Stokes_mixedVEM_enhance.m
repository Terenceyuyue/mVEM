function [uh,ph,info] = Stokes_mixedVEM_enhance(node,elem,pde,bdStruct)
%Stokes_mixedVEM_enhance solves Stokes problem using the mixed virtual element method
% in the lowest order k=2 with enhancement technique used:  
% 
%     u = [u1,u2], f = [f1,f2]     
%     - nu*\Delta(u) - \nabla(p) = f   in Omega
%     div(u) = 0  in Omega
%     Dirichlet boundary condition   u = g on \Gamma
%
%  The mixed formulation is: Find (u,p) \in (V,Q) such that
%
%      \nu*a(u,v) + b(v,p) = (f,v),   v in V,
%      b(u,q) = 0, q in Q,
%
% where,
% - V = (H_0^1(Omega))^2,   Q = L_0^2(Omega)
% - u is approximated by the lowest order divergence-free virtual element
% with enhancement technique used (k = 2)
% - p by piecewise linear element (k = 2)
%      Q_{k-1}(K) = P_{k-1}(K).
%
%   References
%   L. Beirao da Veiga, C. Lovadina and G.Vacca, "Virtual elements for the
%   Navier-Stokes problem on polygonal meshes", SIAM J. Numer. Anal., Vol 56. 
%   No 3., pp. 1210--1242, 2018.
%
% Copyright (C)  Terence Yu. 


%% Gauss-Lobatto weights and points
% sp, wp
sp = [-1, -1/sqrt(5), 1/sqrt(5), 1]; % [-1,1]
sp = (sp+1)/2; % [0,1]: gives the ratios
wp = [1/12, 5/12, 5/12, 1/12];
% rhs
f1 = @(x,y) pde.f([x,y])*[1;0]; % f(p) = f([x,y]), f = [f1, f2]
f2 = @(x,y) pde.f([x,y])*[0;1];

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdofA = (N+NE)*2 + 2*NT;  NNdofB = 3*NT;  NNdof = NNdofA + NNdofB;
Nm = 6; Nmm = 2*Nm;

%% Compute and assemble the linear system
elemLen = cellfun('length',elem); 
nnzA = sum((4*elemLen+2).^2);  nnzB = sum((4*elemLen+2)*3); 
iiA = zeros(nnzA,1); jjA = zeros(nnzA,1); ssA = zeros(nnzA,1);
iiB = zeros(nnzB,1); jjB = zeros(nnzB,1); ssB = zeros(nnzB,1);
elemb = zeros(NNdofA,1);  Fb = zeros(NNdofA,1); 
idA = 0; idB = 0;  ib = 0;
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
dLag = zeros(NT, 3); % for Lagrange multiplier
for iel = 1:NT
    % ------- element information ----------
    index = elem{iel};   indexEdge = elem2edge{iel};   Nv = length(index);
    xK = centroid(iel,1);  yK = centroid(iel,2); 
    hK = diameter(iel);    areaK = area(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % --------------- scaled monomials -----------------
    m1 = @(x,y) 1+0*x;                  gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;             gradm2 = @(x,y) [1/hK+0*x, 0+0*x];
    m3 = @(x,y) (y-yK)./hK;             gradm3 = @(x,y) [0+0*x, 1/hK+0*x];
    m4 = @(x,y) (x-xK).^2/hK^2;         gradm4 = @(x,y) [2*(x-xK)./hK^2, 0+0*x];
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2;   gradm5 = @(x,y) [(y-yK)./hK^2, (x-xK)./hK^2];
    m6 = @(x,y) (y-yK).^2./hK^2;        gradm6 = @(x,y) [0+0*x, 2*(y-yK)./hK^2];
    % we need to add m7, ..., m10 for L2 projection
    m7 = @(x,y) (x-xK).^3/hK^3;         gradm7 = @(x,y) [3*(x-xK).^2/hK^3, 0+0*x];
    m8 = @(x,y) (x-xK).^2.*(y-yK)/hK^3; gradm8 = @(x,y) [2*(x-xK).*(y-yK)/hK^3, (x-xK).^2/hK^3];
    m9 = @(x,y) (x-xK).*(y-yK).^2/hK^3; gradm9 = @(x,y) [(y-yK).^2/hK^3, 2*(x-xK).*(y-yK)/hK^3];
    m10 = @(x,y) (y-yK).^3/hK^3;        gradm10 = @(x,y) [0+0*x, 3*(y-yK).^2/hK^3];
    
    mk1 = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    mk2 = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)];
    mc = {m1,m2,m3,m4,m5,m6,m7,m8,m9,m10};
    Gradmc = {gradm1,gradm2,gradm3,gradm4,gradm5,gradm6,gradm7,gradm8,gradm9,gradm10};   
    divmm = @(x,y) [0+0*x, 1/hK+0*x, 0+0*x, 2*(x-xK)./hK^2, (y-yK)./hK^2, 0+0*x, ...
        0+0*x, 0+0*x, 1/hK+0*x, 0+0*x, (x-xK)./hK^2, 2*(y-yK)./hK^2];
    
    % -------- transition matrix ----------
    NdofBd = 2*Nv; NdofA = 4*Nv+2;
    % D
    divmm2 = @(x,y) divmm(x,y).*repmat(m2(x,y),1,Nmm);
    divmm3 = @(x,y) divmm(x,y).*repmat(m3(x,y),1,Nmm);
    D = zeros(NdofA, Nmm);
    Dbd = [mk2(x,y); mk2(xe,ye)];
    D(1:4*Nv, :) = blkdiag(Dbd, Dbd);
    D(end-1,:) = hK/areaK*integralTri(divmm2,5,nodeT,elemT);
    D(end,:) = hK/areaK*integralTri(divmm3,5,nodeT,elemT);
    % DL: in lifting space
    DL = zeros(NdofA+3,Nmm);
    DL(1:NdofA,:) = D;
    for j = 1:3
        gjv = @(x,y) [repmat(m3(x,y).*mc{j}(x,y),1,Nm).*mk2(x,y), ...
            repmat(-m2(x,y).*mc{j}(x,y),1,Nm).*mk2(x,y)];
        DL(NdofA+j,:) = 1/area(iel)*integralTri(gjv,5,nodeT,elemT);
    end
    
    % --------- elliptic projection -----------
    % B
    % first term
    Lapm = [0, 0, 0, 2/hK^2, 0, 2/hK^2];
    I1 = zeros(Nmm, NdofA);
    I1(1:Nm, end-1) = areaK*Lapm;
    I1(Nm+1:end, end) = areaK*Lapm;
    % second term
    I2 = zeros(Nmm, NdofA);   
    elem1 = [v1(:), v2(:), v1(:)+Nv]; % elem2dof for [ae, be, me]
    Gradmm = cell(2,Nmm);
    Gradmm(1,1:Nm) = Gradmc(1:Nm); 
    Gradmm(2,Nm+1:end) = Gradmc(1:Nm);
    for im = 1:Nm
        Gradmm{1,im+Nm} = @(x,y) [0+0*x, 0+0*x];
        Gradmm{2,im} = @(x,y) [0+0*x, 0+0*x];
    end
    qmm = cell(1,Nmm);  % q = hK(c2*m2 + c3*m3)
    c2 = [Lapm, zeros(1,Nm)];
    c3 = [zeros(1,Nm), Lapm];    
    for im = 1:Nmm
        qmm{im} = @(x,y) hK*(c2(im)*m2(x,y) + c3(im)*m3(x,y));
    end
    for s = 1:2
        id = (1:NdofBd) + (s-1)*NdofBd;
        for im = 1:Nmm
            pm = @(x,y) Gradmm{s,im}(x,y);
            qa = @(x,y) qmm{im}(x,y);
            F1 = 1/6*(sum(pm(x(v1),y(v1)).*Ne, 2) - qa(x(v1),y(v1)).*Ne(:,s));
            F2 = 1/6*(sum(pm(x(v2),y(v2)).*Ne, 2) - qa(x(v2),y(v2)).*Ne(:,s));
            F3 = 4/6*(sum(pm(xe,ye).*Ne, 2) - qa(xe,ye).*Ne(:,s));            
            I2(im, id) = accumarray(elem1(:), [F1; F2; F3], [NdofBd, 1]);
        end
    end
    B = I1 + I2;
    % constraint 
    P0K = zeros(2,NdofA);
    P0K(1, end-1) = -1; P0K(2, end) = -1;
    for s = 1:2
        ms = mc{s+1}; % m2, m3
        F1 = 1/6*(ms(x(v1),y(v1)).*Ne);  % [n1, n2]
        F2 = 1/6*(ms(x(v2),y(v2)).*Ne);
        F3 = 4/6*(ms(xe,ye).*Ne);
        F = [F1; F2; F3];
        P0K(s, 1:NdofBd) = hK/areaK*accumarray(elem1(:), F(:,1), [NdofBd 1]);
        P0K(s, NdofBd+1:2*NdofBd) = hK/areaK*accumarray(elem1(:), F(:,2), [NdofBd 1]);
    end
    % Bs, G, Gs
    Bs = B;  Bs([1,7], :) = P0K;
    G = B*D;  Gs = Bs*D;
    % elliptic projection in lifting space
    %BL = zeros(Nmm,NdofA+3);
    BsL = zeros(Nmm,NdofA+3);
    %BL(:,1:NdofA) = B;
    BsL(:,1:NdofA) = Bs;
    GsL = BsL*DL;
    PisL = GsL\BsL;  PiL = DL*PisL;

    % --------- L2 projection: Pi_k^0 -----------
    % expansion of div(phi)
    H1 = zeros(10,3); % add more rows for later use
    for i = 1:10
        fun = @(x,y) repmat(mc{i}(x,y),1,3).*mk1(x,y);
        H1(i,:) = integralTri(fun,5,nodeT,elemT);
    end
    r = zeros(3,NdofA); % each column for a basis function
    elemr = [elem1(:); elem1(:)+NdofBd];
    Fr = 1/6*[Ne; Ne; 4*Ne];
    r(1,:) = accumarray(elemr, Fr(:), [NdofA, 1]);
    r(2,end-1) = areaK/hK;
    r(3,end) = areaK/hK;
    a = H1(1:3, :)\r;
    % expansion of phi|_e
    A = [1 -1/2 1/4; 1 1/2 1/4; 1 0 0];
    ao = cell(1,Nv);  % \overline{a}
    au = cell(1,Nv);  % \underline{a}
    for j = 1:Nv % loop of edges
        Mato = zeros(3,2*Nv); % vanish for remaining columns
        Mato(1,v1(j)) = 1;
        Mato(2,v2(j)) = 1;
        Mato(3,Nv+j) = 1;
        aoj = zeros(3,NdofA);  auj = zeros(3,NdofA);
        aoj(:,1:2*Nv) = A\Mato;  
        auj(:,2*Nv+(1:2*Nv)) = aoj(:,1:2*Nv);
        ao{j} = aoj;  au{j} = auj;
    end
    % expansion of m_\alpha \in (P_k(K))^2
    Hoplus = zeros(12,12);
    gm1 = @(x,y) [m3(x,y), -m2(x,y)].*repmat(m1(x,y),1,2);
    gm2 = @(x,y) [m3(x,y), -m2(x,y)].*repmat(m2(x,y),1,2);
    gm3 = @(x,y) [m3(x,y), -m2(x,y)].*repmat(m3(x,y),1,2);
    dcmpBase = horzcat(Gradmc(2:10), {gm1,gm2,gm3});
    for i = 1:12
        for j = i:12
            fun = @(x,y) dot(dcmpBase{i}(x,y), dcmpBase{j}(x,y));
            Hoplus(j,i) = integralTri(fun,5,nodeT,elemT);
            Hoplus(i,j) = Hoplus(j,i);
        end
    end
    roplus = zeros(12,12);  % each column for \alpha
    for im = 1:Nm  % loop of alpha
        mo = @(x,y) [mc{im}(x,y), 0*x];  % [m,0]
        mu = @(x,y) [0*x, mc{im}(x,y)];  % [0,m]
        for i = 1:12
            funou = @(x,y) [dot(dcmpBase{i}(x,y), mo(x,y)),...
                dot(dcmpBase{i}(x,y), mu(x,y))];
            roplus(i,[im, im+Nm]) = integralTri(funou,5,nodeT,elemT);
        end
    end
    aoplus = Hoplus\roplus;
    aNabla = aoplus(1:9,:);  aBot = aoplus(10:12,:);

    % C(\alpha,i):
    % first term
    Inabla = -H1(2:end,:)*a;
    for i = 1:Nv % loop of edges
        xp = x(v1(i)) + sp*(x(v2(i))-x(v1(i))); % 4 nodes
        yp = y(v1(i)) + sp*(y(v2(i))-y(v1(i)));
        for s = 2:10
            Hes = zeros(1,3);
            msp = mc{s}(xp,yp);
            for j = 1:3                
                mjep = (sp-0.5).^(j-1);
                Hes(j) = sum(wp.*msp.*mjep);
            end
            Inabla(s-1,:) = Inabla(s-1,:) + Hes*(ao{i}*Ne(i,1)+au{i}*Ne(i,2));
        end
    end
    % second term
    dofPiL = PiL(end-2:end, 1:NdofA);
    Ibot = areaK*dofPiL;
    % combination
    C = zeros(12,NdofA);
    for al = 1:12 % for alpha
        C(al,:) =  C(al,:) + sum(repmat(aNabla(:,al),1,NdofA).*Inabla) ...
            + sum(repmat(aBot(:,al),1,NdofA).*Ibot);
    end
    % H
    H = C*D;

    % ---------- projection matrices -----------------------
    % elliptic projection
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));
    % L2 projection: Pi_k^0
    Pi0s = H\C;    %Pi0 = D*Pi0s;
        
    % ------------- stiffness matrix -------------
    % stiffness matrix A
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    AK = pde.nu*reshape(AK',1,[]); % straighten as row vector for easy assembly
    % stiffness matrix B for Qh
    % BK = zeros(NdofA,3);
    % BK(end-1, 2) = areaK/hK;  BK(end, 3) = areaK/hK;
    % F = 1/6*[(1*Ne); (1*Ne); (4*Ne)]; % [n1, n2]
    % BK(1:end-2,1) = accumarray([elem1(:);elem1(:)+NdofBd], F(:), [2*NdofBd 1]);
    BK = r';
    BK = reshape(BK',1,[]); % straighten as row vector for easy assembly 
    
    
    % -------- load vector ----------
    fun = @(x,y) [repmat(f1(x,y),1,Nm).*mk2(x,y), ...
            repmat(f2(x,y),1,Nm).*mk2(x,y)];
    Pf = integralTri(fun,5,nodeT,elemT);
    fK = Pi0s'*Pf(:);
    
    % ------ assembly index for bilinear forms --------
    NdofA = 4*Nv+2; NdofB = 3;
    indexDofA = [index, indexEdge+N, ...
                 index+N+NE, indexEdge+2*N+NE, ...
                 iel+2*N+2*NE, iel+2*N+2*NE+NT];
    indexDofB = [iel, iel+NT, iel+2*NT];
    iiA(idA+1:idA+NdofA^2) = reshape(repmat(indexDofA, NdofA,1), [], 1);
    jjA(idA+1:idA+NdofA^2) = repmat(indexDofA(:), NdofA, 1);
    ssA(idA+1:idA+NdofA^2) = AK(:);
    idA = idA + NdofA^2;
    iiB(idB+1:idB+NdofA*NdofB) = reshape(repmat(indexDofA, NdofB,1), [], 1);
    jjB(idB+1:idB+NdofA*NdofB) = repmat(indexDofB(:), NdofA, 1);
    ssB(idB+1:idB+NdofA*NdofB) = BK(:);
    idB = idB + NdofA*NdofB;
    
    % ------- assembly index for rhs -------
    elemb(ib+1:ib+NdofA) = indexDofA(:);
    Fb(ib+1:ib+NdofA) = fK(:);
    ib = ib + NdofA;

    % ------------ Lagrange multiplier -------------
    dLag(iel,:) = integralTri(mk1,5,nodeT,elemT);
    
    % ------- matrix for error evaluation -------
    Ph{iel} = Pis; 
    elem2dof{iel} = indexDofA;
end

%% Get block linear system
% Lagrange multiplier 
iid = (1:NNdofB)'+NNdofA;  
jjd = (NNdof+1)*ones(NNdofB,1);
ssd = dLag(:);  
% kk: A, B, B', d, d'
ii = [iiA; iiB;        jjB+NNdofA; iid;  jjd];
jj = [jjA; jjB+NNdofA; iiB;        jjd;  iid];
ss = [ssA; ssB;        ssB;        ssd;  ssd];
kk = sparse(ii,jj,ss,NNdof+1,NNdof+1);
% ff
ff = accumarray(elemb,Fb,[NNdof+1 1]);

%% Apply Dirichlet boundary conditions 
% bdDof, freeDof
bdNodeIdx = bdStruct.bdNodeIdx;
bdEdge = bdStruct.bdEdge;  bdEdgeIdx = bdStruct.bdEdgeIdx;
id = [bdNodeIdx; bdEdgeIdx+N; bdNodeIdx+N+NE; bdEdgeIdx+2*N+NE];
isBdDof = false(NNdof+1,1); isBdDof(id) = true;
bdDof = (isBdDof); freeDof = (~isBdDof);
% bdval
g_D = pde.g_D;
z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:); ze = (z1+z2)./2;
gv = g_D(node(bdNodeIdx,:));
ge = g_D(ze);
bdval = [gv(:,1); ge(:,1); gv(:,2); ge(:,2)]; 
% sol
sol = zeros(NNdof+1,1); 
sol(bdDof) = bdval;
ff = ff - kk*sol;

%% Set solver
sol(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
uh = sol(1:NNdofA); % u = [u1(zv); u1(ze); u2(zv); u2(ze); divu*m2; divu*m3]
ph = sol(NNdofA+1:end-1);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
% info.kk = kk; %info.freeDof = freeDof;
info.D = D;