function [u,info] = PoissonVEMk2_NC(node,elem,pde,bdStruct)
%PoissonVEMk2_NC solves Poisson equation using the standard nonconforming 
% virtual element method in the case of k=2.
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
centroid = aux.centroid;  diameter = aux.diameter; area = aux.area;
% numbers
NT = size(elem,1); NE = size(edge,1); 
Nm = 6;
NNdof = 2*NE + NT;

%% Get elementwise signs of basis functions
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
sgnBase = cell(NT,1);
for iel = 1:NT
    index = elem{iel};   Nv = length(index);   Ndof = 2*Nv+1;
    sgnedge = sign(diff(index([1:Nv,1])));
    id = elem2edge{iel}; sgnbd = E(id); sgnedge(sgnbd) = 1;
    sgnelem = ones(Ndof,1); sgnelem(Nv+1:2*Nv) = sgnedge; 
    sgnBase{iel} = sgnelem;
end

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  
elemLen = cellfun('length',elem); 
nnz = sum((2*elemLen+1).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(2*elemLen+1);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % ------- element information --------
    index = elem{iel};  indexEdge = elem2edge{iel}; Nv = length(index);    
    xK = centroid(iel,1); yK = centroid(iel,2); 
    hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1]; 
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % ------- scaled monomials ----------
    m1 = @(x,y) 1+0*x;                    gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;               gradm2 = @(x,y) [1+0*x, 0+0*x]./hK;
    m3 = @(x,y) (y-yK)./hK;               gradm3 = @(x,y) [0+0*x, 1+0*x]./hK;
    m4 = @(x,y) (x-xK).^2/hK^2;           gradm4 = @(x,y) [2*(x-xK), 0+0*x]./hK^2;
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2;     gradm5 = @(x,y) [(y-yK), (x-xK)]./hK^2;
    m6 = @(x,y) (y-yK).^2./hK^2;          gradm6 = @(x,y) [0+0*x, 2*(y-yK)]./hK^2;     
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)];
    Gradmc = {gradm1, gradm2, gradm3, gradm4, gradm5, gradm6};
    
    % -------- transition matrix ----------
    Ndof = 2*Nv+1;
    D = zeros(Ndof,Nm);
    for i = 1: Nv % loop of edges
        % chi_i , i = 1 ,... , Nv
        D(i,:) = 1/6*( m(x(v1(i)),y(v1(i))) + 4*m(xe(i),ye(i)) +  m(x(v2(i)),y(v2(i))) );
        % chi_{N_v+i} , i = 1 ,... , Nv
        D(Nv+i,:) = 1/12*( -m(x(v1(i)),y(v1(i))) + m(x(v2(i)),y(v2(i))) );
    end
    % chi_{2*Nv+1}
    D(end,:) = 1/area(iel)*integralTri(m,5,nodeT,elemT);
    
    % --------- elliptic projection -----------
    % first term
    Lapm = zeros(6,1); Lapm([4,6]) = 2/hK^2;
    Dof = [zeros(1,2*Nv), area(iel)];
    I1 = Lapm*Dof;
    % second term
    I2 = zeros(Nm, Ndof);
    for im = 1:Nm
        % he*ce0, he*ce1
        ce0 = zeros(1,Nv);  ce1 = zeros(1,Nv);
        for i = 1:Nv
            ce0(i) = dot(Gradmc{im}(xe(i),ye(i)), Ne(i,:));
            ce1(i) = 2*( ce0(i) - dot(Gradmc{im}(x(v1(i)),y(v1(i))), Ne(i,:)) );
        end
        I2(im,:) = [ce0, ce1, 0];
    end
    B = -I1 + I2;
    % constraint
    Bs = B;  Bs(1,:) = Dof;  
    % consistency relation
    G = B*D;  Gs = Bs*D;  

    % -------- sign matrix and sign vector -------
    sgnelem = sgnBase{iel}; 
    sgnK = sgnelem*sgnelem'; 
       
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi)); 
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    AK = AK.*sgnK;
    A = reshape(AK',1,[]); 
    
    % --------- local load vector -----------
    fm = @(x,y) repmat(pde.f([x,y]),1,6).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,5,nodeT,elemT); rhs = rhs';
    fK = Pis'*rhs;  
    fK = fK.*sgnelem;

    % --------- assembly index for ellptic projection -----------
    indexDof = [indexEdge, indexEdge+NE, iel+2*NE];  % local to global index
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = A(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    sgnPis = repmat(sgnelem',size(Pis,1),1);
    Ph{iel} = sgnPis.*Pis; 
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,NNdof,NNdof);
ff = accumarray(elemb,Fb,[NNdof 1]);

% %% Assemble Neumann boundary conditions
% bdEdgeN = bdStruct.bdEdgeN;  bdEdgeIdxN = bdStruct.bdEdgeIdxN;
% if ~isempty(bdEdgeN)
%     g_N = pde.Du;
%     z1 = node(bdEdgeN(:,1),:); z2 = node(bdEdgeN(:,2),:); zc = (z1+z2)/2;
%     e = z1-z2;  % e = z2-z1
%     Ne = [-e(:,2),e(:,1)]; % scaled ne
%     F1 = sum(Ne.*g_N(z1),2);
%     Fc = sum(Ne.*g_N(zc),2);
%     F2 = sum(Ne.*g_N(z2),2);
%     FN = (F1+4*Fc+F2)/6;
%     ff = ff + accumarray(bdEdgeIdxN(:), FN(:),[NNdof 1]);
% end

%% Apply Dirichlet boundary conditions
% boundary information
g_D = pde.g_D;  
bdEdgeD = bdStruct.bdEdgeD;
bdEdgeIdxD = bdStruct.bdEdgeIdxD; 
id = [bdEdgeIdxD; bdEdgeIdxD+NE];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% moments on the boundaries
z1 = node(bdEdgeD(:,1),:);  z2 = node(bdEdgeD(:,2),:); zc = (z1+z2)/2;
dofD0 = 1/6*(g_D(z1) + 4*g_D(zc) + g_D(z2));
dofD1 = 1/12*(-g_D(z1) + g_D(z2));
% rhs
u = zeros(NNdof,1); u(bdDof) = [dofD0(:); dofD1(:)];
ff = ff - kk*u;

%% Set solver
u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm