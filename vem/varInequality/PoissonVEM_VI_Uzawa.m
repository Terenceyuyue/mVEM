function [uh,info] = PoissonVEM_VI_Uzawa(node,elem,pde,bdStruct)
%PoissonVEM_VI_Uzawa solves the simplified friction problem using conforming
% virtual element method in the lowest order case. The problem is iteratively 
% solved by using the Uzawa algorithm, hence the implementation is reduced to
% that of Poisson equation in each iteration.
%
% The problem is
%
%  D.E.
%     -\Delta u + cu = f,  in Omega,
%  Frictional boundary conditions
%     |grad(u)*n| <= g,  u*grad(u)*n + g*|u| = 0 on \Gamma_C
%  Dirichlet boundary conditions
%     u = g_D    on \Gamma_D
%  \partial(\Omega) = \Gamma_C + \Gamma_D
%
%   References
%   [1] F. Wang and H. Wei. Virtual element method for simplified friction problem. 
%       Appl. Math. Lett., 85: 125-131, 2018.
%   [2] B. Wu, F. Wang, and W. Han. Virtual element method for a frictional contact 
%       problem with normal compliance. Commun. Nonlinear Sci. Numer. Simul., 107: 
%       Paper No. 106125, 13 pp., 2022.
%
% Copyright (C)  Terence Yu.

%% Get auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid;  diameter = aux.diameter;  area = aux.area;
N = size(node,1);  NT = size(elem,1);
Nm = 3;

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  Dm = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum(elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % ------- element information --------
    index = elem{iel};  Nv = length(index);    
    xK = centroid(iel,1); yK = centroid(iel,2); 
    hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    v1 = 1:Nv;  v2 = [2:Nv,1]; 
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    nodeT = [node(index,:); centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % ------- scaled monomials ----------
    m1 = @(x,y) 1+0*x;                gradm1 = @(x,y) [0+0*x, 0+0*x];
    m2 = @(x,y) (x-xK)./hK;           gradm2 = @(x,y) [1+0*x, 0+0*x]./hK;
    m3 = @(x,y) (y-yK)./hK;           gradm3 = @(x,y) [0+0*x, 1+0*x]./hK;
    m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];   
    mc = {m1,m2,m3};  % cell 
    Gradmc = {gradm1, gradm2, gradm3};    
    
    % -------- transition matrix ----------
    D = m(x,y);   Dm{iel} = D;
    
    % --------- elliptic projection -----------
    % first term  = 0
    B = zeros(Nm, Nv);
    % second term   
    elem1 = [v1(:), v2(:)];
    for im = 1:Nm
        gradmc = Gradmc{im};
        F1 = 0.5*sum(gradmc(x(v1), y(v1)).*Ne, 2);
        F2 = 0.5*sum(gradmc(x(v2), y(v2)).*Ne, 2);
        F = [F1, F2];
        B(im, :) = accumarray(elem1(:), F(:), [Nv 1]);
    end  
    % constraint
    Bs = B;  Bs(1,:) = 1/Nv; 
    % consistency relation
    G = B*D;  Gs = Bs*D;  
    
    % --------- L2 projection -----------   
    H = zeros(Nm,Nm);
    for i = 1:Nm
        fun = @(x,y) repmat(mc{i}(x,y),1,Nm).*m(x,y);
        H(i,:) = integralTri(fun,3,nodeT,elemT);
    end
       
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi)); 
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    BK = pde.c*Pis'*H*Pis + pde.c*hK^2*(I-Pi)'*(I-Pi);  % Pis = Pi0s;
    AB = reshape(AK'+BK',1,[]); 
    
    % --------- load vector -----------
    fK = Pis'*[pde.f(centroid(iel,:))*area(iel);0;0];        

    % --------- assembly index for ellptic projection -----------
    indexDof = index;  Ndof = Nv;  % local to global index
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AB(:);
    ia = ia + Ndof^2;
    
    % --------- assembly index for right hand side -----------
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,N,N);
ff = accumarray(elemb,Fb,[N 1]);

%% Initialization of Uzawa iteration or Gamma_C is empty 
g_D = pde.g_D;  bdNodeIdxD = bdStruct.bdNodeIdxD;
isBdNode = false(N,1); isBdNode(bdNodeIdxD) = true;
bdDof = find(isBdNode); freeDof = find(~isBdNode);
nodeD = node(bdDof,:); 
uD = g_D(nodeD); 
uh0 = zeros(N,1); 
uh0(bdDof) = uD(:);
uh0(freeDof) = kk(freeDof,freeDof)\ff(freeDof);

%% Uzawa iteration
bdEdgeFri = bdStruct.bdEdgeN; 
if ~isempty(bdEdgeFri)
    g = pde.g_C;  
    z1 = node(bdEdgeFri(:,1),:); z2 = node(bdEdgeFri(:,2),:);
    gC = max([g(z1); g(z2)]);
    e = z2-z1;
    he = sqrt(e(:,1).^2+e(:,2).^2);       
    lambdah0 = ones(N,1); % Lagrange multiplizer 
    tol = 1e-8;  Err = 1;   
    iter = 0;  maxIt = 500;
    while Err>tol || iter<=maxIt
        % ----- assemble frictional boundary conditions -------        
        F1 = 0.5*gC*lambdah0(bdEdgeFri(:,1)).*he;
        F2 = 0.5*gC.*lambdah0(bdEdgeFri(:,2)).*he;
        Fri = [F1,F2];
        ffNew = ff - accumarray(bdEdgeFri(:), Fri(:),[N 1]);
        
        % ------ apply Dirichlet boundary conditions -------
        uh = zeros(N,1); uh(bdDof) = uD(:);
        ffNew = ffNew - kk*uh;
        
        % ------- set solver -------
        uh(freeDof) = kk(freeDof,freeDof)\ffNew(freeDof);
        
        % update the Lagrangian multiplier
        rho = 10;
        lambdah = lambdah0 + rho*gC*uh;
        lambdah = max(-ones(N,1), min(ones(N,1),lambdah));
        
        % compute errors between two steps
        Err = norm(uh-uh0);
        
        % update uh0 and lambdah0
        uh0 = uh; lambdah0 = lambdah;
        iter = iter + 1;
    end
end

%% The discrete solution
uh = uh0;

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm