function [w,info] = Fourth_order_Singular_Perturbation_C0VEM(node,elem,pde,bdStruct)
%This function solves fourth-order singular perturbation using C0-VEM
%
%     \epsilon^2 \Delta^2 u -\Delta u = f,  in Omega
%     Dirichlet boundary condition: 
%                    u = g1,  grad(u)n = g2 on \Gamma, 
%
% Copyright (C)  Terence Yu. 

%% Get auxiliary data
para = pde.para;
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; area = aux.area;
% auxstructure
auxT = auxstructure(node,elem);
elem2edge = auxT.elem2edge; edge = auxT.edge;
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
NNdof = N+2*NE;  % total dof number

%% Compute projection matrices
% Bbs, Gb, Gbs % b: biharmonic
% Bps, Gp, Gps % p: Poisson
Ph = cell(NT,1); % matrix for error evaluation
elem2dof = cell(NT,1);  
elemLen = cellfun('length',elem); 
nnz = sum((3*elemLen).^2);
ii = zeros(nnz,1); jj = zeros(nnz,1); ss = zeros(nnz,1); 
nnz = sum(3*elemLen);
elemb = zeros(nnz,1); Fb = zeros(nnz,1);
ia = 0; ib = 0;
bdEdgeIdx = bdStruct.bdEdgeIdx;  E = false(NE,1); E(bdEdgeIdx) = 1;
for iel = 1:NT
    % 0. ------------- element information ----------------
    index = elem{iel};   indexEdge = elem2edge{iel};
    Nv = length(index);  Ndof = 3*Nv; 
    xK = centroid(iel,1); yK = centroid(iel,2);  hK = diameter(iel);
    x = node(index,1); y = node(index,2); % vertices
    v1 = 1:Nv; v2 = [2:Nv,1]; % loop index for vertices or edges
    p1 = [Nv,1:Nv-1]; p2 = 1:Nv; % jump index
    xe = (x(v1)+x(v2))/2;  ye = (y(v1)+y(v2))/2; % mid-edge points
    Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % scaled outer normal vectors
    he = sqrt(sum(Ne.^2,2));
    ne = Ne./repmat(he,1,2); % outer normal vectors
    te = [-ne(:,2), ne(:,1)];
    n1 = ne(:,1)';  n2 = ne(:,2)'; 
    t1 = te(:,1)';  t2 = te(:,2)';
    % scaled monomials
    m = @(x,y) [1+0*x, (x-xK)/hK, (y-yK)/hK, (x-xK).^2/hK^2, ...
               (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2]; % m1,...,m6
    Gradm = @(x,y) [[0,0]; [1, 0]/hK; [0, 1]/hK; [2*(x-xK), 0]/hK^2;
                    [(y-yK), (x-xK)]/hK^2; [0, 2*(y-yK)]/hK^2]; % grad m
    
    % 1. --------------- transition matrices --------------
    D = zeros(Ndof,6);
    D(1:2*Nv,:) = [m(x,y); m(xe,ye)];  % vertices and mid-edge points
    for i = 1:Nv
        gradi = 0.5*(Gradm(x(v1(i)),y(v1(i)))+Gradm(x(v2(i)),y(v2(i))));
        D(2*Nv+i,:) = Ne(i,:)*gradi';
    end
    
    % 2. ------ elliptic projection of biharmonic term --------
    % \partial_ij (m)
    D11 = zeros(6,1); D11(4) = 2/hK^2;
    D12 = zeros(6,1); D12(5) = 1/hK^2;
    D22 = zeros(6,1); D22(6) = 2/hK^2;
    % Mij(m)
    nu = 0;  % D = 1, nu = 0 for biharmonic
    M11 = -((1-nu)*D11 + nu*(D11+D22));
    M12 = -(1-nu)*D12;
    M22 = -((1-nu)*D22 + nu*(D11+D22));
    % Mnn(m) on e1,...,eNv
    Mnn = M11*(n1.*n1) + M12*(n1.*n2+n2.*n1) + M22*(n2.*n2);
    % Mtn(m) on e1,...,eNv
    Mtn = M11*(t1.*n1) + M12*(t1.*n2+t2.*n1) + M22*(t2.*n2);
    % B, Bs, G, Gs
    Bb = zeros(6,Ndof);
    for j = 1:Nv % loop of edges
        % nphi on ej
        nphi = zeros(1,Ndof); nphi(2*Nv+j) = 1;
        % Jump at zj
        Jump = Mtn(:,p2(j))-Mtn(:,p1(j));
        % phi at zj
        phi = zeros(1,Ndof); phi(j) = 1;
        % B1
        Bb = Bb - Mnn(:,j)*nphi + Jump*phi;
    end
    Bbs = Bb;
    % first constraint
    Bbs(1,1:Nv) = 1/Nv; 
    % second constraint
    Bbs(2:3,1:Nv) = te([Nv,1:Nv-1],:)' - te';
    Bbs(2:3,2*Nv+1:end) = ne';
    % consistency relation
    Gb = Bb*D;     Gbs = Bbs*D;
    
    % 3. --- Interior d.o.fs of elliptic projection of basis functions ---    
    Di = zeros(Ndof+1,6);
    Di(1:Ndof,:) = D;
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    Di(end,:) = 1/area(iel)*integralTri(m,2,nodeT,elemT); %  interior d.o.f
    Bis = zeros(6,Ndof+1); Bis(:,1:Ndof) = Bbs;
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
        I2 = I2 + (f1+4*f2+f3);
    end
    Bp = -I1 + 1/6*I2;
    % constraint
    Bps = Bp; Bps(1,1:Nv) = 1/Nv;
    Gp = Bp*D;   Gps = Bps*D;

    % 5. --------- sign matrix and sign vector ----------
    sgnedge = sign(diff(index([1:Nv,1])));
    id = elem2edge{iel}; sgnbd = E(id); sgnedge(sgnbd) = 1;
    sgnBase = ones(Ndof,1); sgnBase(2*Nv+1:3*Nv) = sgnedge;
    sgnK = sgnBase*sgnBase';  sgnF = sgnBase;

    % 6. --------- stiffness matrix ----------
    % Projection
    Pibs = Gbs\Bbs;   Pib = D*Pibs;   Ib = eye(size(Pib)); % biharmonic term
    Pips = Gps\Bps;   Pip = D*Pips;   Ip = eye(size(Pip)); % Poisson term     
    % Stiffness matrix
    AK  = Pibs'*Gb*Pibs + hK^(-2)*(Ib-Pib)'*(Ib-Pib);   % biharmonic
    BK  = Pips'*Gp*Pips + (Ip-Pip)'*(Ip-Pip);  % Poisson     
    ABK = (para.epsilon^2*AK+BK).*sgnK;
    ABK = reshape(ABK',1,[]); % straighten

    % 7. --------- load vector ----------
    fK = zeros(Ndof,1);
    P0f = pde.f(centroid(iel,:))*area(iel)/Nv;
    fK(1:Nv,:) = repmat(P0f,Nv,1);    
    fK = (fK.*sgnF)'; % straighten

    % 7. --------- assembly index ----------
    % matrix
    indexDof = [index, indexEdge+N, indexEdge+N+NE];  
    ii(ia+1:ia+Ndof^2) = reshape(repmat(indexDof, Ndof, 1), [], 1);
    jj(ia+1:ia+Ndof^2) = repmat(indexDof(:), Ndof, 1);
    ss(ia+1:ia+Ndof^2) = ABK(:);
    ia = ia + Ndof^2;
    % vector
    elemb(ib+1:ib+Ndof) = indexDof(:);
    Fb(ib+1:ib+Ndof) = fK(:);
    ib = ib + Ndof;
    
    % 8. --------- matrix for L2 and H1 error evaluation  ---------
    sgnPis = repmat(sgnBase',size(Pibs,1),1);
    Ph{iel} = sgnPis.*Pibs; 
    elem2dof{iel} = indexDof;
end
kk = sparse(ii,jj,ss,NNdof,NNdof);
ff = accumarray(elemb,Fb,[NNdof 1]);

%% Apply Dirichlet boundary conditions
% boundary information
bdNodeIdx = bdStruct.bdNodeIdx; bdEdgeD = bdStruct.bdEdgeD;
g_D = pde.g_D;  Dw = pde.Du;
id = [bdNodeIdx; bdEdgeIdx+N; bdEdgeIdx+N+NE];
isBdNode = false(NNdof,1); isBdNode(id) = true;
bdDof = (isBdNode); freeDof = (~isBdNode);
% vertices on the boundary
pD = node(bdNodeIdx,:); wD = g_D(pD);
% mid-edge on the boundary
z1 = node(bdEdgeD(:,1),:); z2 = node(bdEdgeD(:,2),:); zc = (z1+z2)./2;
wDc = g_D(zc);
% moments on the boundary
eb = z1-z2;  % e = z2-z1
neb = [-eb(:,2),eb(:,1)];  % scaled ne
wnD = sum(1/6*(Dw(z1)+4*Dw(zc)+Dw(z2)).*neb,2);
% rhs
w = zeros(NNdof,1); w(bdDof) = [wD; wDc; wnD];
ff = ff - kk*w;

%% Set solver
% solver = 'amg';
% if N < 2e3, solver = 'direct'; end
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
info.kk = kk; info.DofI = freeDof;