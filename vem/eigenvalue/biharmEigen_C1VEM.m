function lam = biharmEigen_C1VEM(node,elem,bdStruct)

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
node2elem = auxT.node2elem;
% number
N = size(node,1); NT = size(elem,1); NNdof = 3*N;   Nm = 6;
% characteristic lengths of derivatives at nodes
hxi = cellfun(@(id) mean(diameter(id)), node2elem);

%% Compute projection matrices
% Bbs, Gb, Gbs % b: biharmonic
% Bps, Gp, Gps % p: Poisson
D = cell(NT,1);
Bbs = cell(NT,1); Gb = cell(NT,1); Gbs = cell(NT,1); % b: biharmonic
Bps = cell(NT,1); Gp = cell(NT,1); Gps = cell(NT,1); % p: Poisson
for iel = 1:NT
    % 0. ------------- element information ----------------
    index = elem{iel};  Nv = length(index);  Ndof = 3*Nv;
    xK = centroid(iel,1); yK = centroid(iel,2);  hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    hxiK = hxi(index); % characteristic length
    % scaled monomials
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
        (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    Gradm = @(x,y) [[0,0]; [1, 0]/hK; [0, 1]/hK; [2*(x-xK), 0]/hK^2;
                    [(y-yK), (x-xK)]/hK^2; [0, 2*(y-yK)]/hK^2]; % grad m
    % Dx(m'), Dy(m')
    Dxm = @(x,y) [0*x, 1/hK+0*x, 0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0*x];
    Dym = @(x,y) [0*x, 0*x, 1/hK+0*x, 0*x, (x-xK)/hK^2, 2*(y-yK)/hK^2];

    % 1. --------------- transition matrices --------------
    D1 = zeros(Ndof,Nm);
    D1(1:Nv,:) = m(x,y);
    D1(Nv+1:2*Nv,:) = repmat(hxiK,1,Nm).*Dxm(x,y);
    D1(2*Nv+1:end,:) = repmat(hxiK,1,Nm).*Dym(x,y);
    D{iel} = D1;

    % 2. ------ elliptic projection of plate bending term --------
    % \partial_ij (m)
    D11 = zeros(Nm,1); D11(4) = 2/hK^2;
    D12 = zeros(Nm,1); D12(5) = 1/hK^2;
    D22 = zeros(Nm,1); D22(6) = 2/hK^2;
    % Mij(m)
    nu = 0.3;
    M11 = -((1-nu)*D11 + nu*(D11+D22));
    M12 = -(1-nu)*D12;
    M22 = -((1-nu)*D22 + nu*(D11+D22));

    n1 = ne(:,1);  n2 = ne(:,2);
    Mnn = M11*(n1.*n1)' + M12*(n1.*n2+n2.*n1)' + M22*(n2.*n2)';
    % Mtn(m) on e1,...,eNv
    t1 = te(:,1); t2 = te(:,2);
    Mtn = M11*(t1.*n1)' + M12*(t1.*n2+t2.*n1)' + M22*(t2.*n2)';

    Dxphi = zeros(Nv,Ndof); Dyphi = zeros(Nv,Ndof);
    Dxphi(:,Nv+1:2*Nv) = eye(Nv)./repmat(hxiK,1,Nv);
    Dyphi(:,2*Nv+1:end) = eye(Nv)./repmat(hxiK,1,Nv);
    B1 = zeros(Nm,Ndof);
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv;
    % B, Bs, G, Gs
    for j = 1:Nv % loop of edges
        Dnphi1 = Dxphi(v1(j),:)*n1(j) + Dyphi(v1(j),:)*n2(j); % zj
        Dnphi2 = Dxphi(v2(j),:)*n1(j) + Dyphi(v2(j),:)*n2(j); % z_{j+1}
        nphi = 0.5*he(j)*(Dnphi1+Dnphi2);
        % Jump(m) at zj
        Jump = Mtn(:,p2(j))-Mtn(:,p1(j));
        % phi' at zj
        phi = zeros(1,Ndof);  phi(j) = 1;
        B1 = B1 - Mnn(:,j)*nphi + Jump*phi;
    end
    B1s = B1;
    % first constraint
    B1s(1,1:Nv) = 1/Nv;
    % second constraint
    for j = 1:Nv % loop of edges
        B1s(2:3,1:Nv) = te([Nv,1:Nv-1],:)' - te';
        Dnphi1 = Dxphi(v1(j),Nv+1:end)*n1(j) + Dyphi(v1(j),Nv+1:end)*n2(j); % zj
        Dnphi2 = Dxphi(v2(j),Nv+1:end)*n1(j) + Dyphi(v2(j),Nv+1:end)*n2(j); % z_{j+1}
        Nphi = 0.5*Ne(j,:)'*(Dnphi1+Dnphi2); % scaled
        B1s(2:3,Nv+1:end) = B1s(2:3,Nv+1:end) + Nphi;
    end

    Bbs{iel} = B1s;
    Gb{iel} = B1*D1;     Gbs{iel} = B1s*D1;

    % 3. --- Interior d.o.fs of elliptic projection of basis functions ---
    Di = zeros(Ndof+1,Nm);
    Di(1:Ndof,:) = D1;
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    Di(end,:) = 1/area(iel)*integralTri(m,4,nodeT,elemT); %  interior d.o.f
    Bis = zeros(6,Ndof+1); Bis(:,1:Ndof) = B1s;
    % G, Gs
    Gis = Bis*Di;
    % Pi, Pis
    Pis = Gis\Bis;   Pi = Di*Pis;
    Dof = Pi(end,1:end-1);

    % 4. --------- elliptic projection for Poisson term ----------
    % first term
    Lapm = zeros(6,1); Lapm([4,6]) = 2/hK^2;
    Dof1 = area(iel)*Dof;
    I1 = Lapm*Dof1;
    % second term
    I2 = 0;
    for j = 1:Nv % loop of edges
        % he*\partial_n(m) at x_j, xe, x_j+1
        DnL = Gradm(x(j),y(j))*Ne(j,:)';
        Dne = Gradm(xe(j),ye(j))*Ne(j,:)';
        DnR = Gradm(x(v2(j)),y(v2(j)))*Ne(j,:)'; 
        % [ei,0,0]
        e1 = zeros(1,Ndof); e1(j) = 1;
        e2 = zeros(1,Ndof); e2(Nv+j) = 1;
        e3 = zeros(1,Ndof); e3(v2(j)) = 1;
        % f
        f1 = DnL*e1; f2 = Dne*e2; f3 = DnR*e3;
        I2 = I2 + 1/2*(f1+0*f2+f3);  % -------- 积分修改
    end
    B2 = -I1 + I2; 
    % constraint
    B2s = B2; B2s(1,1:Nv) = 1/Nv;
    Bps{iel} = B2s;
    Gp{iel} = B2*D{iel}; Gps{iel} = B2s*D{iel};
end

%% Get elementwise stiffness matrix and load vector
Aelem = cell(NT,1); Belem = cell(NT,1);
for iel = 1:NT
    % Projection
    hxiK = hxi(elem{iel});
    hh = diag(repmat(hxiK.^(-2), 3, 1));    
    Pibs = Gbs{iel}\Bbs{iel};   Pib = D{iel}*Pibs;   Ib = eye(size(Pib)); % biharmonic term
    Pips = Gps{iel}\Bps{iel};   Pip = D{iel}*Pips;   Ip = eye(size(Pip)); % Poisson term     
    % Stiffness matrix
    AK  = Pibs'*Gb{iel}*Pibs + (hh*(Ib-Pib))'*(Ib-Pib); % biharmonic
    BK  = Pips'*Gp{iel}*Pips + (Ip-Pip)'*(Ip-Pip);    % Poisson 
    Aelem{iel} = reshape(AK',1,[]); % straighten
    Belem{iel} = reshape(BK',1,[]);
end

%% Assemble stiffness matrix and load vector
elemLen = cellfun('length',elem);
nnz = sum((3*elemLen).^2);
vertNum = unique(elemLen);
ii = zeros(nnz,1); jj = zeros(nnz,1);
ssA = zeros(nnz,1);  ssB = zeros(nnz,1);
id = 0;
for Nv = vertNum(:)' % only valid for row vector

    % assemble the matrix
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv); % number of elements with Nv vertices

    elemNv = cell2mat(elem(idNv)); % elem

    elem2 = [elemNv, elemNv+N, elemNv+2*N];
    KA = cell2mat(Aelem(idNv));  KB = cell2mat(Belem(idNv));
    Ndof = 3*Nv;
    ii(id+1:id+NTv*Ndof^2) = reshape(repmat(elem2, Ndof,1), [], 1);
    jj(id+1:id+NTv*Ndof^2) = repmat(elem2(:), Ndof, 1);
    ssA(id+1:id+NTv*Ndof^2) = KA(:);
    ssB(id+1:id+NTv*Ndof^2) = KB(:);
    id = id + NTv*Ndof^2;
end
kkA = sparse(ii,jj,ssA,NNdof,NNdof);
kkB = sparse(ii,jj,ssB,NNdof,NNdof);

%% Apply Dirichlet boundary conditions
% boundary information 
bdNodeIdx = bdStruct.bdNodeIdx; 
% bdDof, freeDof
id = [bdNodeIdx; bdNodeIdx+N; bdNodeIdx+2*N];
isBdNode = false(NNdof,1); isBdNode(id) = true;
freeDof = (~isBdNode);

%% Set solver
A = kkA(freeDof,freeDof);
B = kkB(freeDof,freeDof);

[~, lam] = eigs(A,B,1,'smallestabs');