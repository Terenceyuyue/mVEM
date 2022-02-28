clc; close all; 
clear variables;

%% Parameters
nameV = [200, 400, 600, 800, 1000]./2;
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);
ErrH2 = zeros(maxIt,1);
ErrI = zeros(maxIt,1);

%% PDE data
pde = PlateBendingData();

%% Virtual element method
for k = 1:maxIt
    % load mesh
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    [w,info] = PlateBending_C0VEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(w);  h(k) = 1/sqrt(size(elem,1));
    figure(1); 
    showresult(node,elem,pde.uexact,w);
    drawnow; %pause(0.1);
    % compute errors in discrete norms    
    kOrder = 2;
    ErrL2(k) = getL2error(node,elem,w,info,pde,kOrder);
    ErrH1(k) = getH1error(node,elem,w,info,pde,kOrder);
    ErrH2(k) = getH2error(node,elem,w,info,pde,kOrder);
    % compute errors in discrete energy norm
    % auxstructure
    auxT = auxstructure(node,elem);
    edge = auxT.edge;
    % chi1: evaluation at all vertices
    uexact = pde.uexact;
    chi1 = uexact(node);
    % chi2: evaluation at all mid-edge points
    z1 = node(edge(:,1),:); z2 = node(edge(:,2),:); zc = (z1+z2)/2;
    chi2 = uexact(zc);
    % chi3: moments of \partial_n v on edges with given orientation
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)];  % scaled ne
    Dw = pde.Du;
    Dw1 = Dw(z1); Dwc = Dw(zc); Dw2 = Dw(z2);
    wnD = sum(1/6*(Dw1+4*Dwc+Dw2).*Ne,2);
    chi3 = wnD;
    % chi
    chi = w;
    chie = [chi1; chi2; chi3];
    % Relative error
    kk = info.kk; DofI = info.DofI;
    kk = kk(DofI,DofI); chie = chie(DofI); chi = chi(DofI);
    ErrI(k) = sqrt(abs((chi-chie)'*kk*(chi-chie)));
end

%% Plot convergence rates and display error table
figure(2);
showrateErr(h,ErrL2,ErrH1,ErrH2,ErrI);

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1','|u-u_h|_2','||u_I-u_h||_A'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e',ErrH2,'%0.5e',ErrI,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the H2-norm (1st order), H1-norm (2nd order) and L2-norm
% (2nd order) is observed for the lowest order k = 2. 

