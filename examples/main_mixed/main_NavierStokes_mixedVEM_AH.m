clc; close all;
clear variables;

%% Parameters
nameV = [32, 64, 128, 256, 512];
maxIt = length(nameV);
maxj = maxIt;
h = zeros(maxj,1);   N = zeros(maxj,1);
ErruL2 = zeros(maxj,1);
ErruH1 = zeros(maxj,1);
ErrpL2 = zeros(maxj,1);

%% PDE data
nu = 1/1; 
id = 1;
pde = NavierStokesdata(nu,id);
rho = 1/(2*nu);  alpha = rho^2;
para.rho = rho; para.alpha = alpha;

%% Virtual element method
for k = 1:maxj
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    % [node,elem] = squaremesh([0 1 0 1],1/(5*k),1/(5*k),'tri');
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve the problem
    [uh,ph,info] = NavierStokes_mixedVEM_AH(node,elem,pde,bdStruct,para);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1);
    showresult(node,elem,pde.uexact,uh);
    % figure(2);
    % clf;
    % plot(info.Err, 'k', 'linewidth', 1)
    % [uhI,nodeI,elemI] = EllipticProjection(node,elem,uh,info,2);
    % showresult(nodeI,elemI,pde.uexact,uhI);
    % [phI,nodeI,elemI] = PolyProjection(node,elem,ph,1);
    % showresult(nodeI,elemI,pde.pexact,phI);
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
showrateh(h,ErruH1,'r-o','ErruH1', ...
            ErrpL2, 'b-s','ErrpL2');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||Du-Du_h||','||p-p_h||'};
format shorte
disptable(colname,N,[],h,'%0.3e',ErruL2,'%0.5e',ErruH1,'%0.5e',ErrpL2,'%0.5e');

