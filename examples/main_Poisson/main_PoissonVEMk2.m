clc; close all; 
clear variables;

%% Parameters
nameV = [100, 200, 300, 400, 500];
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);
ErrI = zeros(maxIt,1);

%% PDE data
c = 1;
pde = Poissondata1(c);
bdNeumann = 'abs(x-0)<1e-4 | abs(x-1)<1e-4'; % string for Neumann

%% Virtual element method
for k = 1:maxIt
    % load mesh
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    % get boundary information
    bdStruct = setboundary(node,elem,bdNeumann);
    % solve
    [uh,info] = PoissonVEMk2(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1); 
    showresult(node,elem,pde.uexact,uh);
    pause(0.1);
    % compute errors in discrete L2 and H1 norms
    kOrder = 2; 
    ErrL2(k) = getL2error(node,elem,uh,info,pde,kOrder);
    ErrH1(k) = getH1error(node,elem,uh,info,pde,kOrder);
    % compute errors in discrete energy norm
    % auxgeometry
    aux = auxgeometry(node,elem);
    % auxstructure
    auxT = auxstructure(node,elem);
    edge = auxT.edge;
    % chi1: evaluation at all vertices
    uexact = pde.uexact;
    chi1 = uexact(node);
    % chi2: evaluation at all mid-edge points
    z1 = node(edge(:,1),:); z2 = node(edge(:,2),:); zc = (z1+z2)/2;
    chi2 = uexact(zc);
    % chi3: moments on element K
    NT = size(elem,1);
    chi3 = zeros(NT,1);
    fun = @(x,y) uexact([x,y]);
    for iel = 1:NT
        index = elem{iel};     Nv = length(index);     Ndof = 2*Nv+1;
        nodeT = [node(index,:);aux.centroid(iel,:)];
        elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
        chi3(iel) = 1/aux.area(iel)*integralTri(fun,3,nodeT,elemT);
    end
    % errors in energy norm
    kk = info.kk; freeDof = info.freeDof;
    kk = kk(freeDof,freeDof);
    chie = [chi1; chi2; chi3]; chi = uh;
    chie = chie(freeDof); chi = chi(freeDof);
    ErrI(k) = sqrt(abs((chi-chie)'*kk*(chi-chie))/abs(chie'*kk*chie));
end

%% Plot convergence rates and display error table
figure(2);
showrateErr(h,ErrL2,ErrH1,ErrI);

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1','||u_I-u_h||_E'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e',ErrI,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (2nd order) and L2-norm
% (3rd order) is observed for k = 2. The order of ||uI-uh||_E is 2.7 order and
% thus superconvergence exists.

figure,spy(info.kk)