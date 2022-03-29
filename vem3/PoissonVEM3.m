function [u,info] = PoissonVEM3(node3,elem3,pde,bdStruct)
%PoissonVEM3 solves Poisson equation in 3-D
%
%     -\Delta u + cu = f,  in Omega (c=0)
%     Dirichlet boundary condition u=g_D on \Gamma_D,
%     Neumann boundary condition   grad(u)*n=g_N on \Gamma_N
%
%   References
%   L. Beirao da Veiga, F. Dassi and A. Russo, "High-order virtual element method on 
%   polyhedral meshes", Comput. Math. Appl., Vol 74. No. 5, pp. 1110¨C1122, 2017.
%
% Copyright (C)  Terence Yu.

%% Get auxiliary data
% auxgeometry3
aux = auxgeometry3(node3,elem3);
node3 = aux.node3; elem3 = aux.elem3;
centroid3 = aux.centroid3;  diameter3 = aux.diameter3; 
% auxstructure3
auxT = auxstructure3(node3,elem3);
face = auxT.face;  elem2face = auxT.elem2face;
% numbers
N = size(node3,1); NT = size(elem3,1); NF = length(face);

%% Derive elliptic projections of all faces
faceProj = cell(NF,1);
for s = 1:NF
    % Pifs
    faces = face{s};  P = node3(faces,:);
    Pifs = faceEllipticProjection(P);
    % sort the columns
    [~,idx] = sort(faces);
    faceProj{s} = Pifs(:,idx);
end

%% Compute and assemble the linear system
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);
elemLen = cellfun(@(elemf) length(unique(horzcat(elemf{:}))), elem3);
nnz = sum(elemLen.^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1);
nnz = sum(elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % ------- element information --------
    % faces
    elemf = elem3{iel};   indexFace = elem2face{iel};
    % global index of vertices and local index of elemf
    [~,index3,~,elemfLocal] = faceTriangulation(elemf);
    % centroid and diameter
    Nv = length(index3);  Ndof = Nv;
    V = node3(index3,:);
    xK = centroid3(iel,1); yK = centroid3(iel,2); zK = centroid3(iel,3);
    hK = diameter3(iel);
    x = V(:,1);  y = V(:,2);  z = V(:,3);
    
    % ------- scaled monomials ----------
    m1 = @(x,y,z) 1+0*x;
    m2 = @(x,y,z) (x-xK)/hK;
    m3 = @(x,y,z) (y-yK)/hK;
    m4 = @(x,y,z) (z-zK)/hK;
    m = @(x,y,z) [m1(x,y,z),m2(x,y,z),m3(x,y,z),m4(x,y,z)]; % m1,m2,m3,m4
    mc = {m1,m2,m3,m4};
    gradmMat = [0 0 0; 1/hK 0 0; 0 1/hK 0; 0 0 1/hK];
    
    % -------- transition matrix ----------
    D = m(x,y,z);
    
    % ----------- elliptic projection -------------
    B = zeros(4,Ndof);
    for s = 1:size(elemf,1)
        % --- information of current face
        % vertices of face
        faces = elemf{s};  P = node3(faces,:);
        % elliptic projection on the face
        idFace = indexFace(s);   
        Pifs = faceProj{idFace}; % the order may be not correct
        [~,~,idx] = unique(faces);
        Pifs = Pifs(:,idx);
        % normal vector           
        e1 = P(2,:)-P(1,:);  en = P(1,:)-P(end,:);
        nf = cross(e1,en); nf = nf./norm(nf);
        % area
        areaf = polyarea3(P);
        % --- integral of Pifs
        intFace = [areaf,0,0]*Pifs;  % local
        intProj = zeros(1,Ndof);  % global adjustment
        faceLocal = elemfLocal{s};        
        intProj(faceLocal) = intFace;
        % add grad(m)*nf
        Bf = dot(gradmMat, repmat(nf,4,1), 2)*intProj;
        B = B + Bf;
    end
    % constraint
    Bs = B;  Bs(1,:) = 1/Ndof;
    % consistency relation
    G = B*D;  Gs = Bs*D;
    
    % --------- L2 projection -----------   
    H = zeros(4,4);
    if abs(pde.c)>1e-8
        for i = 1:4
            fun = @(x,y,z) repmat(mc{i}(x,y,z),1,4).*m(x,y,z);
            H(i,:) = integralPolyhedron(fun,3,node3,elemf);
        end
    end
    
    % --------- local stiffness matrix ---------
    Pis = Gs\Bs;   Pi  = D*Pis;   I = eye(size(Pi));
    AK  = Pis'*G*Pis + pde.c*Pis'*H*Pis ...
        + hK*(1+pde.c*hK^2)*(I-Pi)'*(I-Pi); 
    AK = reshape(AK,1,[]); % straighten
    
    % --------- load vector -----------
    %fK = Pis'*[pde.f(centroid3(iel,:))*aux.volume(iel);0;0;0];
    % or use the Gaussian rule
    fun = @(x,y,z) repmat(pde.f([x,y,z]),1,4).*m(x,y,z);
    fK = integralPolyhedron(fun,3,node3,elemf);
    fK = Pis'*fK(:);
    
    % --------- assembly index for ellptic projection -----------
    indexDof = index3;
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = AK(:);
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

%% Assemble Neumann boundary conditions
bdFaceN = bdStruct.bdFaceN;  bdFaceIdxN = bdStruct.bdFaceIdxN;
if ~isempty(bdFaceN)
    Du = pde.Du;
    faceLen = cellfun('length',bdFaceN);
    nnz = sum(faceLen);
    elemb = zeros(nnz,1); FN = zeros(nnz,1);
    ib = 0;
    for s = 1:size(bdFaceN,1)
        % vertices of face
        faces = bdFaceN{s};   nv = length(faces);
        P = node3(faces,:);
        % elliptic projection on the face
        idFace = bdFaceIdxN(s);   
        Pifs = faceProj{idFace}; % the order may be not correct
        [~,idx] = sort(faces);
        Pifs = Pifs(:,idx);
        % 3-D polygon -> 2-D polygon
        poly = localPolygon3(P);        
        nodef = poly.nodef;  % local coordinates
        nf = poly.nf;   % outer normal vector of face
        centroidf = poly.centroidf;
        sc = centroidf(1); tc = centroidf(2);
        hf = poly.diameterf;
        Coord = poly.Coord; % (s,t) ---> (x,y,z)
        % g_N
        g_N = @(s,t) dot(Du(Coord(s,t)), nf);
        fun = @(s,t) g_N(s,t).*[1+0*s, (s-sc)/hf, (t-tc)/hf];
        Ff = integralPolygon(fun,3,nodef);
        Ff = Pifs'*Ff(:);
        % assembly index
        elemb(ib+1:ib+nv) = faces(:);
        FN(ib+1:ib+nv) = Ff(:);
        ib = ib + nv;
    end
    ff = ff + accumarray(elemb(:), FN(:),[N 1]);
end

%% Apply Dirichlet boundary conditions
g_D = pde.g_D;  bdNodeIdxD = bdStruct.bdNodeIdxD;
isBdNode = false(N,1); isBdNode(bdNodeIdxD) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
nodeD = node3(bdDof,:);
u = zeros(N,1); uD = g_D(nodeD); u(bdDof) = uD(:);
ff = ff - kk*u;

%% Set solver
solver = 'amg';
if N < 2e3, solver = 'direct'; end
% solve
switch solver
    case 'direct'
        u(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
    case 'amg'
        option.solver = 'CG';
        u(freeDof) = amg(kk(freeDof,freeDof),ff(freeDof),option);
end

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; info.DofI = freeDof; % d.o.f.s for computing errors in the energy norm
info.aux = aux; % auxgeometry3