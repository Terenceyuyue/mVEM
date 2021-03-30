function [uh,ph,info] = Darcy_mixedVEM(node,elem,pde,bdStruct)
%Darcy_mixedVEM solves Darcy problem using the mixed virtual element method
% in the lowest order
%
%     -div(K*grad(p)) = f   in Omega
%     Neumann boundary condition   K*grad(u)*n=g on \Gamma
%
%  The mixed formulatoin is: 
%      u = K*grad(p), div u = -f  in \Omega
%      u*n = g                    on \Gamma
%  The velocity u = K*grad(p) is approximated by the lowest order
%  H(div)-conforming virtual element and p by piecewise constant element.
%
%   References
%   F. Brezzi, R.S. Falk and L.D. Marini, "Basic principles of mixed
%   virtual element methods", ESAIM Math. Model. Numer. Anal.,
%   Vol 48. No 4., pp. 1227--1240, 2014.
%
% Copyright (C)  Terence Yu. 

K = pde.K;  % coefficient matrix

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
NT = size(elem,1); NE = size(edge,1);
NNdofA = 2*NE+NT;  NNdofB = NT;  NNdof = NNdofA + NNdofB;
Nm = 6;   Nmh = Nm-1;

%% Compute projection matrices
D = cell(NT,1);  Bs = cell(NT,1); 
Gs = cell(NT,1);  
for iel = 1:NT
    % ------- element information ----------
    index = elem{iel};     Nv = length(index);  
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); 
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne
    Te = [-Ne(:,2), Ne(:,1)]; % he*te
    
    % ------- scaled monomials --------     
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;
    m4 = @(x,y) (x-xK).^2./hK^2; 
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2; 
    m6 = @(x,y) (y-yK).^2./hK^2;
    % \hat{m}_a = K*grad(hK*m_{a+1})
    mh1 = @(x,y) [K(1,1)+0*x; K(2,1)+0*x];
    mh2 = @(x,y) [K(1,2)+0*x; K(2,2)+0*x];
    mh3 = @(x,y) [2*K(1,1)*m2(x,y); 2*K(2,1)*m2(x,y)];
    mh4 = @(x,y) [K(1,1)*m3(x,y)+K(1,2)*m2(x,y); K(2,1)*m3(x,y)+K(2,2)*m2(x,y)];
    mh5 = @(x,y) [2*K(1,2)*m3(x,y); 2*K(2,2)*m3(x,y)];
    mh = @(x,y) [mh1(x,y), mh2(x,y), mh3(x,y), mh4(x,y), mh5(x,y)];    

    % -------- transition matrix ----------
    NdofA = 2*Nv+1; 
    D1 = zeros(NdofA,Nmh); 
    for i = 1:Nv % loop of edges
        % v at z_i, z_{i+1} for v = [mK1,...,mK5]
        va = mh(x(v1(i)),y(v1(i))); vb = mh(x(v2(i)),y(v2(i)));
        % chi_i, i = 1,...,Nv
        D1(i,:) =  1/2*Ne(i,:)*(va+vb);
        % chi_{Nv+i}, i = 1,...,Nv
        D1(Nv+i,:) =  1/12*Ne(i,:)*(vb-va);
        % chi_{2Nv+1}
        D1(end,:) = D1(end,:) + 1/2*Te(i,:)*(va+vb);
    end
    D{iel} = D1;
    
    % --------- elliptic projection -----------
    % first term
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    m = @(x,y) [m2(x,y), m3(x,y), m4(x,y), m5(x,y), m6(x,y)]; % m_{a+1},...
    ci = zeros(1,NdofA); ci(1:Nv) = 1/area(iel);
    Intm = integralTri(m,3,nodeT,elemT);
    I1 = Intm'*ci;
    % second term
    rij = zeros(NdofA,Nv); rij(1:Nv,:) = eye(Nv);
    sij = zeros(NdofA,Nv); sij(Nv+1:2*Nv,:) = eye(Nv);
    rija = rij - 6*sij;
    rije = rij;
    rijb = rij + 6*sij;
    I2 = zeros(Nmh,NdofA);
    for j = 1:Nv  
        ma = m(x(v1(j)),y(v1(j)));
        me = m(xe(j),ye(j));
        mb = m(x(v2(j)),y(v2(j)));
        rja = rija(:,j)';  rje = rije(:,j)';  rjb = rijb(:,j)';
        I2 = I2 + 1/6*(ma'*rja + 4*me'*rje + mb'*rjb);
    end
    % Bs=B, Gs=G
    B1 = hK*(-I1+I2);
    Bs{iel} = B1;  Gs{iel} = B1*D1;    
end

%% Get elementwise signs of basis functions
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
sgnBase = cell(NT,1);
for iel = 1:NT
    index = elem{iel};   Nv = length(index);   NdofA = 2*Nv+1;
    sgnedge = sign(diff(index([1:Nv,1])));
    id = elem2edge{iel}; sgnbd = E(id); sgnedge(sgnbd) = 1;
    sgnelem = ones(NdofA,1); sgnelem(1:Nv) = sgnedge; 
    sgnBase{iel} = sgnelem;
end

%% Get elementwise stiffness matrix and load vector
Aelem = cell(NT,1); Belem = cell(NT,1); belem = cell(NT,1);
Ph = cell(NT,1); % matrix for error evaluation
for iel = 1:NT
    % element information
    index = elem{iel};   Nv = length(index);   NdofA = 2*Nv+1;
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    % sign matrix and sign vector
    sgnelem = sgnBase{iel}; 
    sgnK = sgnelem*sgnelem';  
    % Projection
    Pis = Gs{iel}\Bs{iel};   Pi  = D{iel}*Pis;   I = eye(size(Pi));    
    % Stiffness matrix A
    AK  = Pis'*Gs{iel}*Pis + norm(inv(K),'fro')*(I-Pi)'*(I-Pi);  % G = Gs
    AK = AK.*sgnK;
    Aelem{iel} = reshape(AK',1,[]); % straighten as row vector for easy assembly 
    % Stiffness matrix B
    BK = zeros(NdofA,1); BK(1:Nv) = 1;
    BK = BK.*sgnelem;
    Belem{iel} = reshape(BK',1,[]); % straighten as row vector for easy assembly 
    % Load vector f
    fxy = @(x,y) pde.f([x,y]); % f(p) = f([x,y])
    rhs = integralTri(fxy,3,nodeT,elemT); rhs = rhs';
    fK = -rhs;
    belem{iel} = fK'; % straighten as row vector for easy assembly
    % matrix for error evaluation
    sgnPis = repmat(sgnelem',size(Pis,1),1);
    Ph{iel} = sgnPis.*Pis; 
end

%% Assemble matrices A,B, and vector f,d
elemLen = cellfun('length',elem); vertNum = unique(elemLen);
nnzA = sum((2*elemLen+1).^2);  nnzB = sum((2*elemLen+1)*1); 
iiA = zeros(nnzA,1); jjA = zeros(nnzA,1); ssA = zeros(nnzA,1);
iiB = zeros(nnzB,1); jjB = zeros(nnzB,1); ssB = zeros(nnzB,1);
ffB = zeros(NNdofB,1); 
idA = 0; idB = 0; 
elem2dof = cell(NT,1);
for Nv = vertNum(:)' % only valid for row vector
    
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices
    
    % elem2dof    
    elem2edgeNv = cell2mat(elem2edge(idNv)); 
    elem2dofA = [elem2edgeNv,elem2edgeNv+NE,idNv+2*NE];
    elem2dofB = idNv;
    NdofA = 2*Nv+1; NdofB = 1;
    
    % assemble the matrix A    
    KA = cell2mat(Aelem(idNv));    
    iiA(idA+1:idA+NTv*NdofA^2) = reshape(repmat(elem2dofA, NdofA,1), [], 1);
    jjA(idA+1:idA+NTv*NdofA^2) = repmat(elem2dofA(:), NdofA, 1);
    ssA(idA+1:idA+NTv*NdofA^2) = KA(:);
    idA = idA + NTv*NdofA^2;
    
    % assemble the matrix B    
    KB = cell2mat(Belem(idNv));    
    iiB(idB+1:idB+NTv*NdofA*NdofB) = reshape(repmat(elem2dofA, NdofB,1), [], 1);
    jjB(idB+1:idB+NTv*NdofA*NdofB) = repmat(elem2dofB(:), NdofA, 1);
    ssB(idB+1:idB+NTv*NdofA*NdofB) = KB(:);
    idB = idB + NTv*NdofA*NdofB;

    % assemble the vector
    FB = cell2mat(belem(idNv));
    ffB = ffB +  accumarray(elem2dofB(:),FB(:),[NNdofB 1]);
    
    % elementwise global indices
    elem2dof(idNv) = mat2cell(elem2dofA, ones(NTv,1), NdofA);
end
A = sparse(iiA,jjA,ssA,NNdofA,NNdofA);
B = sparse(iiB,jjB,ssB,NNdofA,NNdofB);
d = area;  % for Lagrange multiplier 

%% Get block linear system
kk = sparse(NNdof+1,NNdof+1);  ff = zeros(NNdof+1,1);
kk(1:NNdofA,1:NNdofA) = A;
kk(1:NNdofA, (1:NNdofB)+NNdofA) = B;
kk((1:NNdofB)+NNdofA, 1:NNdofA) = B';
kk((1:NNdofB)+NNdofA, end) = d;
kk(end, (1:NNdofB)+NNdofA) = d';
ff((1:NNdofB)+NNdofA) = ffB;

%% Apply Dirichlet boundary conditions 
% bdDof, freeDof
bdEdge = bdStruct.bdEdge;  bdEdgeIdx = bdStruct.bdEdgeIdx;
id = [bdEdgeIdx; bdEdgeIdx+NE];
isBdDof = false(NNdof+1,1); isBdDof(id) = true;
bdDof = (isBdDof); freeDof = (~isBdDof);
% bdval
u = pde.uexact;
z1 = node(bdEdge(:,1),:); z2 = node(bdEdge(:,2),:); ze = (z1+z2)./2;
e = z1-z2;  % e = z2-z1
Ne = [-e(:,2),e(:,1)];  
chi1 = 1/6*sum((u(z1)+4*u(ze)+u(z2)).*Ne,2); % u*n = g
chi2 = 1/12*sum((u(z2)-u(z1)).*Ne,2);
bdval = [chi1;chi2]; 
% sol
sol = zeros(NNdof+1,1); 
sol(bdDof) = bdval;
ff = ff - kk*sol;

%% Set solver
sol(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
uh = sol(1:NNdofA); % u = [u1,u2]
ph = sol(NNdofA+1:end-1);

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; %info.freeDof = freeDof;

