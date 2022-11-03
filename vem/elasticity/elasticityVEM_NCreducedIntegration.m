function [u,info] = elasticityVEM_NCreducedIntegration(node0,elem0,varargin)
%function [u,info] = elasticityVEM_NCreducedIntegration(node0,elem0,pde,bdNeumann,refineType)
%function [u,info] = elasticityVEM_NCreducedIntegration(node0,elem0,node,elem,pde,bdStruct,refineType)

if isnumeric(varargin{1}) % node
    [node,elem,pde,bdStruct,refineType] = deal(varargin{:});
else
    [pde,bdNeumann,refineType] = deal(varargin{:});
end

%elasticityVEM_NCreducedIntegration solves linear elasticity equation of tensor
%form using nonconforming virtual element method in V1 with reduced integration applied.
%
%       u = [u1, u2]
%       -div sigma(u) = f in \Omega,
%       Dirichlet boundary condition u = [g1_D, g2_D] on \Gamma_D,
%       Neumann boundary condition  sigma(u)n = [g1_N, g2_N] on \Gamma_N
%
%   References
%   R. S. Falk, "Nonconforming finite element methods for the equations of
%   linear elasticity", Math. Comp., Vol 57. No 196., pp. 529¨C550, 1991.
%
% Copyright (C)  Terence Yu.

%% PDE data
fxy = @(x,y) pde.f([x,y]); % [f1,f2]
para = pde.para;

%% Information of coarse mesh
% auxgeometry
aux0 = auxgeometry(node0,elem0);
node0 = aux0.node; elem0 = aux0.elem;
% numbers
NT0 = size(elem0,1);  

%% Information of the refined mesh
if ~isnumeric(varargin{1}) 
    % refine the mesh
    [node,elem] = PolyMeshUniformRefine(node0,elem0,refineType);
    % boundary information
    bdStruct = setboundary(node,elem,bdNeumann);
end
% auxstructure
aux = auxstructure(node,elem);
edge = aux.edge; 
elem2edge = aux.elem2edge;  
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
% numbers
NT = size(elem,1);  NE = size(edge,1);
Nm = 3;

%% Compute and assemble L2 rot projection
elemLen0 = cellfun('length',elem0); 
nnz = sum((2*(2*elemLen0)).^2); % e -- e1+e2
ii0 = zeros(nnz,1); jj0 = zeros(nnz,1); ss0 = zeros(nnz,1);
ia = 0;  s = 0;
for iK = 1:NT0
    % the structure of the refined mesh
    if refineType<=2
        Ns = elemLen0(iK); % number of the needed subcells in K 
        Nc = 2;  % number of the edges that has contribution
    else % refineType==3
        Ns = 1;  Nc = 2*elemLen0(iK);
    end
    
    % loop of the sub-elements in K
    Crotx = zeros(Ns,Nc);  Croty = zeros(Ns,Nc);  
    indexrot = zeros(Ns,Nc);
    for i = 1:Ns
        iel = s+i;
        index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
        x = node(index,1); y = node(index,2);
        v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
        Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
        Te = [-Ne(:,2), Ne(:,1)];
        % Crot
        Crotx(i,:) = Te(1:Nc,1);  % The first two edges of T are the bisected edges
        Croty(i,:) = Te(1:Nc,2); 
        indexrot(i,:) = indexEdge(1:Nc);
    end
    s = s+Ns;
    
    % L2 rot projection
    H0K = aux0.area(iK);
    C0x = reshape(Crotx',1,[]);
    C0y = reshape(Croty',1,[]);
    C0K = [C0x, C0y];     
    PirotK = H0K\C0K; 
    BrotK = PirotK'*H0K*PirotK;
    
    % local to global index
    Ndof = Nc*Ns; Ndof2 = 2*Ndof; 
    indexDofx = reshape(indexrot',1,[]);
    indexDof = [indexDofx, indexDofx+NE];

    % assembly index
    ii0(ia+1:ia+Ndof2^2) = reshape(repmat(indexDof, Ndof2, 1), [], 1);
    jj0(ia+1:ia+Ndof2^2) = repmat(indexDof(:), Ndof2, 1);
    ss0(ia+1:ia+Ndof2^2) = reshape(BrotK',[],1);
    ia = ia + Ndof2^2;
end
ss0 = -para.mu*ss0; 

%% Compute and assemble elliptic projection and L2 div projection
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
    
    % -------- elliptic projection ------------
    % B
    B = zeros(Nm,Nv);
    for i = 1:Nv % loop of edges
        % gradm*Ne
        gi = sum(Gradm.*repmat(Ne(i,:),3,1), 2);
        B(:,i) = gi;
    end  
    % Bs
    Bs = zeros(Nm,Nv);
    Bs(2:end,:) = B(2:end,:);
    Bs(1,:) = he;
    % consistency relation
    G = B*D;  Gs = Bs*D;  
    
    % -------- L2 projection ------------
    H0 = area(iel);
    C01x = Ne(:,1)';  C01y = Ne(:,2)';
    C0 = [C01x, C01y];
    
    % --------- local load vector -----------
    fint = integralTri(fxy,3,nodeT,elemT); % [f1,f2]
    fK = repmat(fint,Nv,1); 
    fK = 1/Nv*fK(:);  
       
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi)); 
    Pi0s = H0\C0;
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    AK = 2*para.mu*blkdiag(AK,AK);
    BK = para.lambda*Pi0s'*H0*Pi0s;
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
kk = sparse([ii;ii0],[jj;jj0],[ss;ss0],2*NE,2*NE);
ff = accumarray(elemf,Ff,[2*NE 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN; bdEdgeIdxN = bdStruct.bdEdgeIdxN;
if ~isempty(bdEdgeN)
    g_N = pde.g_N;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    Sig1 = g_N(z1); Sig2 = g_N(z2);
    F11 = sum(Ne.*Sig1(:,[1,3]),2);    F12 = sum(Ne.*Sig2(:,[1,3]),2); % g1 at endpoints
    F21 = sum(Ne.*Sig1(:,[3,2]),2);    F22 = sum(Ne.*Sig2(:,[3,2]),2);
    F1 = (F11+F12)/2;  F2 = (F21+F22)/2;
    FN = [F1, F2];
    ff = ff + accumarray([bdEdgeIdxN(:); bdEdgeIdxN(:)+NE], FN(:),[2*NE 1]);
end

%% Lagrange multiplier for pure traction problem
bdEdgeD = bdStruct.bdEdgeD; 
nd = 0; % initialized as zero
if isempty(bdEdgeD)
    % information of boundary
    nd = 3; 
    bdEdge = bdStruct.bdEdge;
    bdEdgeIdx = bdStruct.bdEdgeIdx;
    z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:);
    Te = z2-z1;
    % d1,d2
    d1 = zeros(2*NE,1); d1(bdEdgeIdx) = 1;
    d2 = zeros(2*NE,1); d2(bdEdgeIdx+NE) = 1;
    % d3
    d3 = zeros(2*NE,1); 
    d3([bdEdgeIdx;bdEdgeIdx+NE]) = Te(:);
    % kkd, ffd
    kkd = sparse(2*NE+nd, 2*NE+nd);
    kkd(1:2*NE,1:2*NE) = kk;
    kkd(1:2*NE, 2*NE+(1:nd)) = [d1,d2,d3];
    kkd(2*NE+(1:nd), 1:2*NE) = [d1';d2';d3'];
    ffd = zeros(2*NE+nd,1);
    ffd(1:2*NE) = ff;
    % kk, ff
    kk = kkd;  ff = ffd;
end

%% Apply Dirichlet boundary conditions
bdEdgeD = bdStruct.bdEdgeD;  bdEdgeIdxD = bdStruct.bdEdgeIdxD;
g_D = pde.g_D;
isBdNode = false(2*NE+nd,1);
isBdNode([bdEdgeIdxD, bdEdgeIdxD+NE]) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
z1 = node(bdEdgeD(:,1),:);  z2 = node(bdEdgeD(:,2),:);
uD = 0.5*(g_D(z1) + g_D(z2));
u = zeros(2*NE+nd,1); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
u = u(1:2*NE);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof;
info.kk = kk(1:2*NE,1:2*NE); %info.DofI = freeDof;
info.node = node; info.elem = elem;