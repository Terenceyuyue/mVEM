function [u,info] = PoissonVEM_NC(node,elem,pde,bdStruct)
%PoissonVEM_NC solves Poisson equation using the standard nonconforming 
% virtual element method in the lowest order case.
% 
% The problem is 
%
%     -\Delta u = f,  in Omega, 
%     Dirichlet boundary condition u=g_D on \Gamma_D, 
%     Neumann boundary condition   grad(u)*n=g_N on \Gamma_N
%
% Copyright (C)  Terence Yu. 

fxy = @(x,y) pde.f([x,y]);

%% Get auxiliary data
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge;  edge = auxT.edge;
% auxgeometry
aux = auxgeometry(node,elem);
centroid = aux.centroid;  diameter = aux.diameter;
% numbers
NT = size(elem,1); NE = size(edge,1); 
Nm = 3;
NNdof = NE;

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  
elemLen = cellfun('length',elem); 
nnz = sum(elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % ------- element information --------
    index = elem{iel};  indexEdge = elem2edge{iel}; Nv = length(index);    
    xK = centroid(iel,1); yK = centroid(iel,2); 
    hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1]; 
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    he = sqrt(Ne(:,1).^2 + Ne(:,2).^2);
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % ------- scaled monomials ----------
    m1 = @(x,y) 1+0*x;  m2 = @(x,y) (x-xK)./hK;  m3 = @(x,y) (y-yK)./hK;          
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];   
    Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]]; 
    
    % -------- transition matrix ----------
    D = 0.5*( m(x(v1),y(v1)) + m(x(v2),y(v2)) );    
    
    % ---------------- B,Bs,G,Gs -------------------
    % B
    B = zeros(Nm,Nv);     
    for i = 1:Nv % loop of edges        
        gi = sum(Gradm.*repmat(Ne(i,:),3,1), 2); % gradm*Ne
        B(:,i) = gi;  
    end  
    % Bs
    Bs = zeros(Nm,Nv);  
    Bs(2:end,:) = B(2:end,:);
    Bs(1,:) = he;
    % consistency relation
    G = B*D;  Gs = Bs*D;  
       
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi)); 
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    A = reshape(AK',1,[]); 
    
    % --------- local load vector -----------
    fint = integralTri(fxy,3,nodeT,elemT); 
    fK = fint/Nv*ones(Nv,1);       

    % --------- assembly index for ellptic projection -----------
    indexDof = indexEdge;  Ndof = Nv;  % local to global index
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = A(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,NNdof,NNdof);
ff = accumarray(elemb,Fb,[NNdof 1]);

%% Assemble Neumann boundary conditions
bdEdgeN = bdStruct.bdEdgeN;  bdEdgeIdxN = bdStruct.bdEdgeIdxN;
if ~isempty(bdEdgeN)
    g_N = pde.Du;
    z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:);
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)]; % scaled ne
    F1 = sum(Ne.*g_N(z1),2);
    F2 = sum(Ne.*g_N(z2),2);
    FN = (F1+F2)/2;
    ff = ff + accumarray(bdEdgeIdxN(:), FN(:),[NNdof 1]);
end

%% Apply Dirichlet boundary conditions
bdEdgeD = bdStruct.bdEdgeD;  bdEdgeIdxD = bdStruct.bdEdgeIdxD;
g_D = pde.g_D;  
isBdNode = false(NNdof,1); isBdNode(bdEdgeIdxD) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
z1 = node(bdEdgeD(:,1),:);  z2 = node(bdEdgeD(:,2),:);
uD = 0.5*(g_D(z1) + g_D(z2));
u = zeros(NNdof,1); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm