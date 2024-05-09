function [u,info] = elasticityVEM_NavierNC(node,elem,pde,bdStruct)
%elasticityVEM_NavierNC solves linear elasticity equation of Navier form 
% using the nonconforming virtual element method in lowest order case 
%
%       u = [u1, u2]
%       -mu \Delta u - (lambda + mu)*grad(div(u)) = f in \Omega
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D.
%
% Copyright (C)  Terence Yu.

%% PDE data
fxy = @(x,y) pde.f([x,y]); % [f1,f2]
para = pde.para;

%% Get auxiliary data
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge;  edge = auxT.edge;
% auxgeometry
aux = auxgeometry(node,elem);
centroid = aux.centroid;  diameter = aux.diameter; area = aux.area;
% numbers
NT = size(elem,1); NE = size(edge,1); 
Nm = 3;

%% Compute and assemble linear system 
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(2*elemLen);
elemf = zeros(nnz,1); Ff = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % -------- element information ----------
    index = elem{iel};  indexEdge = elem2edge{iel};  Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1];
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    he = sqrt(Ne(:,1).^2 + Ne(:,2).^2);
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % -------- scaled monomials ----------
    m1 = @(x,y) 1+0*x;  m2 = @(x,y) (x-xK)./hK;  m3 = @(x,y) (y-yK)./hK;          
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];   
    Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]];
    
    % ---------- transition matrix ------------
    D = 0.5*(m(x(v1),y(v1)) + m(x(v2),y(v2)));
    
    % ------ integration over edges ---------
    % B
    B = zeros(Nm,Nv); Bs = zeros(Nm,Nv);
    % C0
    C01x = zeros(1,Nv);  C01y = zeros(1,Nv);
    % f
    fK = zeros(2*Nv,1);
    fint = integralTri(fxy,3,nodeT,elemT); % [f1,f2]
    for i = 1:Nv % loop of edges
        % basis for edges
        phie = zeros(1,Nv); phie(i) = 1; % moment values of all basis functions on ei
        phinxe = Ne(i,1)*phie;  phinye = Ne(i,2)*phie;
        % gradm*Ne
        gi = sum(Gradm.*repmat(Ne(i,:),3,1), 2);
        Gradmphie =  gi*phie;
        B = B + Gradmphie;
        % add constraint
        Bs(1,:) = Bs(1,:) + he(i)*phie;
        % C0
        C01x = C01x + phinxe;
        C01y = C01y + phinye;
        % right-hand side
        fK = fK + [fint(1)*phie, fint(2)*phie]'; % [f1,f2]
    end  
    
    % -------- elliptic projection ------------
    % B,Bs
    Bs(2:end,:) = B(2:end,:);
    % consistency relation
    G = B*D;  Gs = Bs*D;  
    
    % -------- L2 projection ------------
    H0 = area(iel);
    C0 = [C01x, C01y];
    
    % --------- local load vector -----------
    fK = 1/Nv*fK;  
       
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi)); 
    Pi0s = H0\C0;
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    AK = para.mu*blkdiag(AK,AK);
    BK = (para.mu+para.lambda)*Pi0s'*H0*Pi0s;
    AB = reshape(AK'+BK',1,[]); % straighten   

    % --------- assembly index for ellptic projection -----------
    indexDof = [indexEdge,indexEdge+NE];  Ndof = length(indexDof);  % local to global index
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemf(ib+1:ib+Ndof) = indexDof(:);
    Ff(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = blkdiag(Pis,Pis);
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,2*NE,2*NE);
ff = accumarray(elemf,Ff,[2*NE 1]);

%% Apply Dirichlet boundary conditions
bdEdgeD = bdStruct.bdEdgeD;  bdEdgeIdxD = bdStruct.bdEdgeIdxD;
g_D = pde.g_D;
isBdNode = false(2*NE,1);
isBdNode([bdEdgeIdxD, bdEdgeIdxD+NE]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
z1 = node(bdEdgeD(:,1),:);  z2 = node(bdEdgeD(:,2),:);
uD = 0.5*(g_D(z1) + g_D(z2));
u = zeros(2*NE,1); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk; %info.DofI = freeDof;