clc; close all; 
clear variables;

%% Parameters
nameV = [32, 64, 128, 256, 512];
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);
ErrI = zeros(maxIt,1);

%% PDE data
c = 1;
pde = Poissondata1(c);
bdNeumann = 'x==0 | x==1'; % string for Neumann

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    % get boundary information
    bdStruct = setboundary(node,elem,bdNeumann);
    % solve
    [uh,info] = PoissonVEMk2(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1); 
    showresult(node,elem,pde.uexact,uh);
    drawnow; %pause(0.1);
    % compute errors in discrete L2, H1 and energy norms
    kOrder = 2; 
    ErrL2(k) = getL2error(node,elem,uh,info,pde,kOrder);
    ErrH1(k) = getH1error(node,elem,uh,info,pde,kOrder);
    ErrI(k) = getErrI(node,elem,uh,info,pde,kOrder);
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

% figure,spy(info.kk)