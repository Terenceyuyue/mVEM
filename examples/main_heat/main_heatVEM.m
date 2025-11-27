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
t0 = 0;  tf = 1;
pde = heatData();
pde.t0 = t0;  pde.tf = tf;
bdNeumann = []; % string for Neumann

% for error
pde1.uexact = @(p) pde.uexact(p,tf);
pde1.Du = @(p) pde.Du(p,tf);

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    % get boundary information
    bdStruct = setboundary(node,elem,bdNeumann);
    % solve
    [uh,info] = heatVEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1); 
    showresult(node,elem,pde1.uexact,uh);
    pause(0.1);
    % compute errors in discrete L2, H1 and energy norms
    kOrder = 1; 

    ErrL2(k) = getL2error(node,elem,uh,info,pde1,kOrder);
    ErrH1(k) = getH1error(node,elem,uh,info,pde1,kOrder);
end

%% Plot convergence rates and display error table
figure(2);
showrateErr(h,ErrL2,ErrH1);

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');