clc; close all;
clear variables;

%% Parameters
nameV = [32, 64, 128, 256, 512];
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% PDE data
para.lambda = 1e8; para.mu = 1;
pde = elasticitydataLocking(para); 
bdNeumann = 'all'; % 'all' for pure traction problem
refineType = 1;

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %[node,elem] = squaremesh([0 1 0 1],1/(5*k),1/(5*k),'tri');
    %load( ['distortionSquare', num2str(5*k), '.mat'] );
    %load( ['distortionPoly', num2str(nameV(k)), '.mat'] );
    
    [uh,info] = elasticityVEM_NCUniformReducedIntegration(node,elem,pde,bdNeumann,refineType);
    node = info.node; elem = info.elem; % the refined mesh
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    [uhI,nodeI,elemI] = EllipticProjection(node,elem,uh,info);     
    if N(k)<4e3
        figure(1);
        showresult(nodeI,elemI,pde.uexact,uhI);
        drawnow;
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
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.4e',ErrH1,'%0.4e');
format shorte
ErrL2'

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (1st order) and L2-norm
% (2nd order) is observed for k = 1.