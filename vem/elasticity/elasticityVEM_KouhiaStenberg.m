function [u,info] = elasticityVEM_KouhiaStenberg(node,elem,pde,bdStruct)
%elasticityVEM_FNC solves linear elasticity equation of tensor form using
% the Kouhia-Stenberg type virtual element: V_nc \times V_c
%
% The problem is
%
%       u = [u1, u2]
%       -div sigma(u) = f in \Omega,
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D,
%       Neumann boundary condition  sigma(u)n = [g1_N, g2_N] on \Gamma_N
%
%   References
%   D. Y. Kwak and H. Park, "Lowest-order virtual element methods for linear 
%   elasticity problems", Math. Comp., Vol 59. No 200., pp. 321¨C338, 1992.
%
% Copyright (C)  Terence Yu.

%% PDE data
f1xy = @(x,y) pde.f([x,y])*[1;0];
f2xy = @(x,y) pde.f([x,y])*[0;1];
para = pde.para;

%% Get auxiliary data
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge;  edge = auxT.edge;
% auxgeometry
aux = auxgeometry(node,elem);
centroid = aux.centroid;  diameter = aux.diameter;  area = aux.area;
% numbers
N = size(node,1); NT = size(elem,1); NE = size(edge,1);

%% Compute projection matrices
%% Also derive load vector and local-global index
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(2*elemLen);
elemf = zeros(nnz,1); Ff = zeros(nnz,1);
ia = 0; ib = 0;

for iel = 1:NT    
    % ------------- element information -------------
    % commonly used information
    index = elem{iel};  indexEdge = elem2edge{iel};
    Nv = length(index); Ndof = 2*Nv;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1];  % starting and ending
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    Ne = [y(v2)-y(v1), x(v1)-x(v2)];
    he = sqrt(Ne(:,1).^2 + Ne(:,2).^2);
    Te = [-Ne(:,2), Ne(:,1)];
    % for integration over auxiliary triangles
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    % scaled monomials
    m1 = @(x,y) 1 + 0*x; m2 = @(x,y) (x-xK)./hK; m3 = @(x,y) (y-yK)./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
    
    % ---------------- transition matrix D ------------------
    D = zeros(Ndof,6);
    D(1:Nv,1:3) = 0.5*( m(x(v1),y(v1)) + m(x(v2),y(v2)) );
    D(Nv+1:end,4:end) = m(x,y);
    
    % -------------- integration over edges ------------------
    % C0
    C01xe = zeros(1,Nv);  C01ye = zeros(1,Nv);  
    C01xv = zeros(1,Nv);  C01yv = zeros(1,Nv);
    % B0 for constraints
    B0x = zeros(3,Nv);  B0y = zeros(3,Nv);
    % fK
    f1K = zeros(Nv,1);
    f1int = integralTri(f1xy,3,nodeT,elemT);
    f2int = integralTri(f2xy,3,nodeT,elemT);
    for i = 1:Nv   % integrating basis functions over edges
        % basis for edges
        phie = zeros(1,Nv); phie(i) = 1; 
        phinxe = Ne(i,1)*phie;  phinye = Ne(i,2)*phie;
        phitxe = Te(i,1)*phie;  
        % C0
        C01xe = C01xe + phinxe;
        C01ye = C01ye + phinye;
        C01xv(i) = 0.5*(Ne(p1(i),1)+Ne(p2(i),1)); % loop of basis
        C01yv(i) = 0.5*(Ne(p1(i),2)+Ne(p2(i),2)); 
        % B0
        B0x(1,:) = B0x(1,:) + phitxe;
        B0x(2,:) = B0x(2,:) + he(i)*phie;
        B0y(1,i) = 0.5*(Te(p1(i),2)+Te(p2(i),2));        
        B0y(3,i) = 0.5*(he(p1(i))+he(p2(i)));
        % right-hand side
        f1K = f1K + phie';
    end
    
    % ------------------ H0,C0 ---------------------
    H0 = area(iel);
    C0 = [C01xe, C01yv];
    
    % ---------------- B,Bs,G,Gs -------------------
    % B0
    B0 = [B0x, B0y];
    % B
    E = zeros(6,4);
    E(2,1) = 1/hK;  E([3,5],[2,3]) = 1/(2*hK); E(6,4) = 1/hK;
    B = [E(:,1)*C01xe+E(:,2)*C01ye, E(:,3)*C01xv+E(:,4)*C01yv];
    % Bs
    Bs = B;
    Bs([1,3,4],:) = B0;   
    % G,Gs
    G = B*D;     Gs = Bs*D;
    
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));
    Pi0s = H0\C0;
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi); 
    AK = 2*para.mu*AK;
    BK = para.lambda*Pi0s'*H0*Pi0s;
    AB = reshape(AK'+BK',1,[]);
    
    % ------------- local load vector ------------------
    f1K = f1int/Nv*f1K;
    f2K = f2int*ones(Nv,1)/Nv;
    fK = [f1K; f2K];
    
    % --------- assembly index for ellptic projection -----------
    indexDof = [indexEdge+N, index];  Ndof = length(indexDof);
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemf(ib+1:ib+Ndof) = indexDof(:);
    Ff(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,N+NE,N+NE);
ff = accumarray(elemf,Ff,[N+NE 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN; bdEdgeIdxN = bdStruct.bdEdgeIdxN;
if ~isempty(bdEdgeN)
    g_N = pde.g_N;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    Sig1 = g_N(z1); Sig2 = g_N(z2);  % Sig = [Sig11, Sig22, Sig12]
    F11 = sum(Ne.*Sig1(:,[1,3]),2); % g1 at z1 
    F12 = sum(Ne.*Sig2(:,[1,3]),2); % g1 at z2
    F21 = sum(Ne.*Sig1(:,[3,2]),2)./2;    
    F22 = sum(Ne.*Sig2(:,[3,2]),2)./2;
    Fe = (F11+F12)/2;  Fv = [F21,F22];
    FN = [Fe,Fv];
    ff = ff + accumarray([bdEdgeIdxN(:)+N; bdEdgeN(:)], FN(:),[N+NE 1]);
end

%% Apply Dirichlet boundary conditions
bdEdgeD = bdStruct.bdEdgeD;  bdEdgeIdxD = bdStruct.bdEdgeIdxD;
bdNodeIdxD = bdStruct.bdNodeIdxD;
g_D = pde.g_D;
isBdNode = false(N+NE,1);
isBdNode([bdEdgeIdxD+N; bdNodeIdxD]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
z1 = node(bdEdgeD(:,1),:);  z2 = node(bdEdgeD(:,2),:);
uDe = 0.5*(g_D(z1) + g_D(z2)); uDe = uDe(:,1); % u1
nodeD = node(bdNodeIdxD,:);
uDv = g_D(nodeD); uDv = uDv(:,2); % u2
uD = [uDe; uDv];
u = zeros(N+NE,1); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; %info.DofI = freeDof;
