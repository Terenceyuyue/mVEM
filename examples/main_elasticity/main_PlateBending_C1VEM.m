clc; close all; 
clear variables;

%% Parameters
nameV = [1000, 2000, 3000, 4000, 5000];
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);
ErrH2 = zeros(maxIt,1);

%% PDE data
pde = PlateBendingData();

%[node,elem] = squaremesh([0 1 0 1], 0.5,0.5);
%% Virtual element method
for k = 1:maxIt
    % load mesh
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %[node,elem] = uniformrefine(node,elem);
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    option.solver = 'direct';
    [w,info] = PlateBending_C1VEM(node,elem,pde,bdStruct,option);
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
end

%% Plot convergence rates and display error table
figure(2);
showrateErr(h,ErrL2,ErrH1,ErrH2);

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1','|u-u_h|_2'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e',ErrH2,'%0.5e');

%% Conclusion
% 
% The optimal rate of convergence of the H2-norm (1st order), H1-norm (2nd order) and L2-norm (2nd order) is observed 
% for the lowest order k = 2 when meshes are fine enough or the domain is refined by triangulation.

