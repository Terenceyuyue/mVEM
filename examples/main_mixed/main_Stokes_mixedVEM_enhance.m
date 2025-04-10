clc; close all;
clear variables;

%% Parameters
nameV = [32, 64, 128, 256, 512];
maxIt = length(nameV);
h = zeros(maxIt,1);   N = zeros(maxIt,1);
ErruL2 = zeros(maxIt,1);
ErruH1 = zeros(maxIt,1);
ErrpL2 = zeros(maxIt,1);

%% PDE data
pde = Stokesdata;

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve the problem
    [uh,ph,info] = Stokes_mixedVEM_enhance(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1);
    showresult(node,elem,pde.uexact,uh);
    %[uhI,nodeI,elemI] = EllipticProjection(node,elem,uh,info,2);
    %showresult(nodeI,elemI,pde.uexact,uhI);
    %[phI,nodeI,elemI] = PolyProjection(node,elem,ph,1);
    %showresult(nodeI,elemI,pde.pexact,phI);
    drawnow; %pause(0.1);
    % compute errors in discrete L2 norm
    kOrder = 2; pOrder = 1;
    ErruL2(k) = getL2error(node,elem,uh,info,pde,kOrder);
    ErruH1(k) = getH1error(node,elem,uh,info,pde,kOrder);
    ErrpL2(k) = getL2error_Poly(node,elem,ph,pde,pOrder);
end

%% Plot convergence rates and display error table
figure,
subplot(1,2,1);
showrateh(h,ErruL2,'r-*','||u-u_h||', ...
            ErruH1,'b-o','||Du-Du_h||');
subplot(1,2,2);
showrateh(h,ErruH1,'r-o','||Du-Du_h||', ...
            ErrpL2, 'b-s','||p-p_h||');
        
fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||','||p-p_h||'};
disptable(colname,N,[],h,'%0.3e',ErruL2,'%0.5e',ErruH1,'%0.5e',ErrpL2,'%0.5e');