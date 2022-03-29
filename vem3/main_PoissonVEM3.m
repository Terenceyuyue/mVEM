clc; close all; 
clear variables;

%% Parameters
maxIt = 5;
h = zeros(maxIt,1);      N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);  ErrH1 = zeros(maxIt,1);

%% PDE data
c = 1;
pde = Poisson3data(c);
bdNeumann = 'x==0'; % string for Neumann

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['SimpleMesh3data', num2str(k), '.mat'] ); % polyhedral mesh
    %load( ['mesh3data', num2str(k), '.mat'] ); % polyhedral mesh
    %[node3,~,elem3] = cubemesh([0 1 0 1 0 1], 1/(2*k)); % tetrahedral mesh
    % get boundary information
    bdStruct = setboundary3(node3,elem3,bdNeumann);
    % solve
    [uh,info] = PoissonVEM3(node3,elem3,pde,bdStruct);
    % record
    N(k) = length(uh);  h(k) = (1/size(elem3,1))^(1/3);
    % compute errors in discrete L2, H1 and energy norms
    kOrder = 1; 
    [ErrH1(k),ErrL2(k)] = getError3(node3,elem3,uh,info,pde,kOrder);
end

%% Plot convergence rates and display error table
figure, showrateErr(h,ErrL2,ErrH1);

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed for k = 1. 