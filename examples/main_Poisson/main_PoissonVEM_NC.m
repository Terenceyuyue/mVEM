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
c = 0; 
pde = Poissondata1(c); % c must be zero
bdNeumann = 'x==0 | x==1'; % string for Neumann

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %[node,elem] = squaremesh([0 1 0 1],1/(5*k),1/(5*k),'tri');
    % get boundary information
    bdStruct = setboundary(node,elem,bdNeumann);
    % solve the problem    
    [uh,info] = PoissonVEM_NC(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
	[uhI,nodeI,elemI] = EllipticProjection(node,elem,uh,info);
    figure(1); 
    showresult(nodeI,elemI,pde.uexact,uhI);
    drawnow; 
    % compute errors in discrete L2 and H1 norms
    kOrder = 1; 
    ErrL2(k) = getL2error(node,elem,uh,info,pde,kOrder);
    ErrH1(k) = getH1error(node,elem,uh,info,pde,kOrder);
end

%% Plot convergence rates and display error table
figure(2);
showrateErr(h,ErrL2,ErrH1);

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed for k = 1. 