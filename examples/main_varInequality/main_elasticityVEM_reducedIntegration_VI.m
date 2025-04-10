clc; close all; 
clear variables;

%% Parameters
nameV = [100, 200, 300, 400, 500];
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% PDE data
para.lambda = 1e10; para.mu = 1;
pde = elasticitydataLockingVI(para); %elasticitydataVI();
GammaC = 'y==0';  % frictional boundary
GammaD = 'x==1';  % Dirichlet boundary
bdString = {GammaC, GammaD}; % the remaining part is Neumann boundary
constraintType = 1;
refineType = 1;

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %[node,elem] = squaremesh([0 1 0 1],1/(5*k),1/(5*k),'tri');
    [uh,info] = elasticityVEM_reducedIntegration_VI(node,elem,pde,bdString,constraintType,refineType);
    node = info.node; elem = info.elem; % the refined mesh
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    if N(k)<4e3
        figure(1);
        showresult(node,elem,pde.uexact,uh);
        drawnow;
        %pause(0.1);
    end
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
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e');
% 
% %% Conclusion
% %
% % The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% % (2nd order) is observed for k = 1. 
% 
% figure,spy(info.kk)