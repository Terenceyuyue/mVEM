function [w,info] = PlateBending_C1VEM(node,elem,pde,bdStruct,option)
%PlateBending_C1VEM solves plate bending problem using C1-VEM in the lowest
% order case:
%
%       -D_{ij} M_{ij}(w) = f in \Omega,
%       Dirichlet boundary condition:
%               w = g1, grad(w)n = g2    on \Gamma.
%   References
%   [1] F. Brezzi and L.D. Marini, "Virtual Element Methods for plate being
%   problems", Comput. Methods Appl. Mech. Engrg., Vol 253., pp. 455¨C462, 2013.
%   [2] C. Chinosi and L.D. Marini, "Virtual Element Methods for fourth
%   order problems: L^2-estimates", Computer and Mathematices with
%   Applications, Vol 72., pp. 1959-1967, 2016.
%
%  Note: 
%  - The problem under consideration is in the general form.
%  - In the lowest order case k = 2: The right-hand side is approximated by
%  the L^2 projection, which is exactly the H^2 elliptic projection when we
%  consider the enhancement technique (see (3.23) in Ref. [2]).
% 
%
% Copyright (C)  Terence Yu. 


%% Get auxiliary data
para = pde.para; 
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; 
% auxstructure
auxT = auxstructure(node,elem);
node2elem = auxT.node2elem;
% number
N = size(node,1); NT = size(elem,1); NNdof = 3*N;   Nm = 6;
% characteristic lengths at nodes
hxi = cellfun(@(id) mean(diameter(id)), node2elem);

%% Compute projection matrices
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  
elemLen = cellfun('length',elem); 
nnz = sum((3*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(3*elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
for iel = 1:NT
    % --------- element information ---------
    index = elem{iel};    Nv = length(index);   Ndof = 3*Nv;
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % edge index
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    hxiK = hxi(index); % characteristic length 
    nodeT = [node(index,:);centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % ------- scaled monomials ---------
    % m'
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
                (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    % Dx(m'), Dy(m')
    mx = @(x,y) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x];
    my = @(x,y) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2];
    
    % ------ transition matrix ---------
    D = zeros(Ndof,Nm);
    D(1:Nv,:) = m(x,y);
    D(Nv+1:2*Nv,:) = repmat(hxiK,1,Nm).*mx(x,y); 
    D(2*Nv+1:3*Nv,:) = repmat(hxiK,1,Nm).*my(x,y);
    
    % --------- elliptic projection -----------
    % \partial_ij (m)
    D11 = zeros(Nm,1); D11(4) = 2/hK^2;
    D12 = zeros(Nm,1); D12(5) = 1/hK^2;
    D22 = zeros(Nm,1); D22(6) = 2/hK^2;
    % Mij(m)
    M11 = -para.D*((1-para.nu)*D11 + para.nu*(D11+D22));
    M12 = -para.D*(1-para.nu)*D12;
    M22 = -para.D*((1-para.nu)*D22 + para.nu*(D11+D22));
    % Dx(phi'), Dy(phi') at z1,...,zNv (each row)
    phix = zeros(Nv,Ndof); phiy = zeros(Nv,Ndof);
    phix(:,Nv+1:2*Nv) = diag(1./hxiK);  
    phiy(:,2*Nv+1:3*Nv) = diag(1./hxiK);
    % B, Bs, G, Gs
    B = zeros(Nm,Ndof);
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    for i = 1:Nv % loop of edges or vertices
        % Mnn(m)
        n1 = ne(i,1);  n2 = ne(i,2);
        Mnn = M11*(n1*n1) + M12*(n1*n2+n2*n1) + M22*(n2*n2); % (Nm,1)
        % int[\partial_n (phi')] on ei
        nphi1 = phix(v1(i),:)*n1 + phiy(v1(i),:)*n2; % zi
        nphi2 = phix(v2(i),:)*n1 + phiy(v2(i),:)*n2; % z_{i+1}
        nphi = 0.5*he(i)*(nphi1+nphi2);
        % Jump(m) at zi
        tn11 = te(p2(i),1)*ne(p2(i),1) - te(p1(i),1)*ne(p1(i),1); % jump
        tn22 = te(p2(i),2)*ne(p2(i),2) - te(p1(i),2)*ne(p1(i),2);
        tn12 = (te(p2(i),1)*ne(p2(i),2) + te(p2(i),2)*ne(p2(i),1)) ...
            - (te(p1(i),1)*ne(p1(i),2) + te(p1(i),2)*ne(p1(i),1));
        Jump = M11*tn11 + M12*tn12 + M22*tn22;
        % phi' at zi
        phi = zeros(1,Ndof);  phi(i) = 1;
        % B1 on e and at zj
        B = B - Mnn*nphi + Jump*phi;
    end
    Bs = B;
    % first constraint
    Bs(1,1:Nv) = 1/Nv;
    % second constraint
    Bs(2:3,1:Nv) = te(p1,:)' - te(p2,:)';
    for i = 1:Nv % loop of edges 
        n1 = ne(i,1);  n2 = ne(i,2);
        nphi1 = phix(v1(i),Nv+1:3*Nv)*n1 + phiy(v1(i),Nv+1:3*Nv)*n2; % zi
        nphi2 = phix(v2(i),Nv+1:3*Nv)*n1 + phiy(v2(i),Nv+1:3*Nv)*n2; % z_{i+1}
        Nphi = 0.5*Ne(i,:)'*(nphi1+nphi2); % scaled        
        Bs(2:3,Nv+1:3*Nv) = Bs(2:3,Nv+1:3*Nv) + Nphi;
    end
    % consistency relation
    G = B*D;     Gs = Bs*D;
    
    % --------- local stiffness matrix --------- 
    hh = diag(repmat(hxiK.^(-2), 3, 1));
    Pis = Gs\Bs;   Pi = D*Pis;   I = eye(Ndof);
    AK  = Pis'*G*Pis + (hh*(I-Pi))'*(I-Pi);
    A = reshape(AK',1,[]); 
    
    % --------- load vector -----------
    fm = @(x,y) repmat(pde.f([x,y]),1,Nm).*m(x,y); % f(p) = f([x,y])
    rhs = integralTri(fm,5,nodeT,elemT); rhs = rhs';
    fK = Pis'*rhs;
    
    % --------- assembly index for ellptic projection -----------
    indexDof = [index, index+N, index+2*N];  
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

%% Apply Dirichlet boundary conditions
% boundary information 
bdNodeIdx = bdStruct.bdNodeIdx; 
hxib = hxi(bdNodeIdx);% characteristic length on the boundary
% bdDof, freeDof
id = [bdNodeIdx; bdNodeIdx+N; bdNodeIdx+2*N];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% values on the boundary
nodeD = node(bdNodeIdx,:); g_D = pde.g_D;  wD = g_D(nodeD);  
% derivatives on the boundary with characteristic length
Dw = pde.Du;  
hDw = Dw(nodeD).*repmat(hxib,1,2);   % wxD = Dw(:,1); wyD = Dw(:,2);
w = zeros(NNdof,1); w(bdDof) = [wD; hDw(:)];
ff = ff - kk*w;

%% Set solver
% solver = option.solver;
% if NNdof < 2e3, solver = 'direct'; end
% solve
solver = 'direct';
switch solver
    case 'direct'
        w(freeDof) = kk(freeDof,freeDof)\ff(freeDof);
    case 'amg'
        option.solver = 'CG';
        w(freeDof) = amg(kk(freeDof,freeDof),ff(freeDof),option);                 
end

%% Store information for computing errors
info.Ph = Ph; info.elem2dof = elem2dof; 
info.kk = kk; %info.freeDof = freeDof;