clc; clear; close all; 

%% Parameters
nameV = [32, 64, 128, 256, 512]; % 1000:1000:5000;
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
[ErrH1,ErrH2] = deal(zeros(maxIt,1));

%% PDE data
id = 1; % id = 1,2
pde = biharmonicdata_IPVEM(id);
pde.a = 1; 
bdNeumann = [];

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    % get boundary information
    bdStruct = setboundary(node,elem,bdNeumann);
    % solve
    [uh,info] = IPVEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1); 
    showresult(node,elem,pde.uexact,uh);
    drawnow; %pause(0.1);
    % compute errors in discrete norms
    kOrder = 2; 
    infoH1 = info; infoH1.Ph = info.Ph(:,1);
    infoH2 = info; infoH2.Ph = info.Ph(:,2);
    ErrH1(k) = getH1error(node,elem,uh,infoH1,pde,kOrder);
    ErrH2(k) = getH2error(node,elem,uh,infoH2,pde,kOrder);
end

%% Plot convergence rates and display error table
figure(2);
showrateErr(h,ErrH1,ErrH2);

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|u-u_h|_1','|u-u_h|_2'};
disptable(colname,N,[],h,'%0.3e',ErrH1,'%0.5e',ErrH2,'%0.5e');