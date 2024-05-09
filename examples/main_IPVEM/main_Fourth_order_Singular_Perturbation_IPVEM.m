clc; clear; close all; 

%% Parameters
nameV = [32, 64, 128, 256, 512];
%nameV = [100, 200, 300, 400, 500]; % Lshape
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
[ErrH1,ErrH2] = deal(zeros(maxIt,1));

%% PDE data
% epsilon should be small for id = 3,4,5: boundary layers
id = 1;  
epsilon = 10^(-0);
pde = Fourth_order_Singular_Pertubation_Data_IPVEM(id,epsilon); 
pde.a = 10;
bdNeumann = [];

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %load( ['distortionPoly', num2str(nameV(k)), '.mat'] );
    %load( ['meshdata_Lshape', num2str(nameV(k)), '.mat'] ); % Lshape
    % get boundary information
    bdStruct = setboundary(node,elem,bdNeumann);
    % solve
    [uh,info] = Fourth_order_Singular_Perturbation_IPVEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1); 
    showresult(node,elem,pde.uexact,uh);
    drawnow; %pause(0.1);
    % compute errors in discrete L2, H1 and energy norms
    kOrder = 2; 
    infoH1 = info; infoH1.Ph = info.Ph(:,1);
    infoH2 = info; infoH2.Ph = info.Ph(:,2);
    ErrH1(k) = getH1error(node,elem,uh,infoH1,pde,kOrder);
    ErrH2(k) = getH2error(node,elem,uh,infoH2,pde,kOrder);
end
% error in the mesh-dependent norm
f2 = @(x,y) sum(pde.f([x,y]).^2,2);
fsquare = sqrt(squareint(f2,[0,1,0,1]));
ErrEps = sqrt(epsilon^2*ErrH2.^2 + ErrH1.^2)/fsquare;

%% Plot convergence rates and display error table
figure(2);
showrateErr(h,ErrEps);

fprintf('\n')
fprintf('Errors for convergence rate:\n');
fprintf('\n')
fprintf('   &');
fprintf(' %0.4e   &', ErrEps);
fprintf('\n')

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|u-u_h|_1','|u-u_h|_2'};
disptable(colname,N,[],h,'%0.3e',ErrH1,'%0.5e',ErrH2,'%0.5e');