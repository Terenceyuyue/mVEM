function Eigen = PoissonEigenVEM(node,elem,bdStruct)

%% Get auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid;  diameter = aux.diameter;  
N = size(node,1);  NT = size(elem,1);
Nm = 3;

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  Dm = cell(NT,1);
elemLen = cellfun('length',elem); 
nnz = sum(elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); 
ssA = zeros(nnz,1);  ssB = zeros(nnz,1);
ia = 0; 
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
    D = m(x,y);   
    
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
        H(i,:) = integralTri(fun,5,nodeT,elemT);
    end
       
    % --------- local stiffness matrix --------- 
    Pis = Gs\Bs;   Pi = D*Pis;   I = eye(size(Pi)); 
    AK  = Pis'*G*Pis + (I-Pi)'*(I-Pi);
    BK = Pis'*H*Pis + hK^2*(I-Pi)'*(I-Pi);  % Pis = Pi0s;
    A = reshape(AK',1,[]); 
    B = reshape(BK',1,[]);      

    % --------- assembly index for ellptic projection -----------
    indexDof = index;  Ndof = Nv;  % local to global index
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ssA(ia+1:ia+Ndof^2) = A(:);
    ssB(ia+1:ia+Ndof^2) = B(:);
    ia = ia + Ndof^2;
    
    
    % --------- matrix for L2 and H1 error evaluation  ---------
    Ph{iel} = Pis;
    elem2dof{iel} = indexDof;
end
kkA = sparse(ii,jj,ssA,N,N);
kkB = sparse(ii,jj,ssB,N,N);

%% Apply Dirichlet boundary conditions
bdNodeIdx = bdStruct.bdNodeIdx;
isBdNode = false(N,1); isBdNode(bdNodeIdx) = true;
freeDof = (~isBdNode);

%% Set solver
A = kkA(freeDof,freeDof);
B = kkB(freeDof,freeDof);


%Xi= zeros(N,1);
[~,lam] = eigs(A,B,1,'smallestabs');
lam = diag(lam);

lamExact = pi^2*2;

Eigen = [lamExact, lam(:)'];

