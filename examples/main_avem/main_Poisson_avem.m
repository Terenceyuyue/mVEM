clc;clear;close all;

%% Parameters
maxN = 1e4;     theta = 0.4;    maxIt = 30;
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Generate an initial mesh
load meshex1

%% Get the PDE data
pde = Poissondata_avem();

%% Adaptive Virtual Element Method
etaN = zeros(maxIt,1);
for k = 1:maxIt
    % Step 1: SOLVE
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    [uh,info] = PoissonVEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1);
    showresult(node,elem,pde.uexact,uh);
    
    %     figure(1); clf;
    %     showmesh(node,elem,options);
    %     figure(2); clf;
    %     showsolution(node,elem,uh);
    %     figure(3); clf;
    %     showsolution(node,elem,pde.uexact(node));
    
    %     % compute errors in discrete L2 and H1 norms
    %     kOrder = 1;
    %     ErrL2(k) = getL2error(node,elem,uh,info,pde,kOrder);
    %     ErrH1(k) = getH1error(node,elem,uh,info,pde,kOrder);
    
    % Step 2: ESTIMATE
    eta = PoissonVEM_indicator(node,elem,uh,info,pde);
    etaN(k) = sqrt(sum(eta.^2));
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    
    % Step 4: REFINE
    [node,elem] = PolyMeshRefine(node,elem,elemMarked);
    
    if (size(node,1)>maxN) || (k==maxIt)
        bdStruct = setboundary(node,elem);
        uh = PoissonVEM(node,elem,pde,bdStruct);
        step = k
        break;
    end
end

figure,
plot((1:step),etaN(1:step),'k.-');
xlabel('k'); ylabel('\eta (u_h)');