clc;clear;close all;

%% Parameters
maxN = 1e4;     theta = 0.9;    maxIt = 5;
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErruL2 = zeros(maxIt,1);
ErruH1 = zeros(maxIt,1);
ErrpL2 = zeros(maxIt,1);

%% Generate an initial mesh
load meshdata1

%% Get the PDE data
pde = Stokesdata();

%% Adaptive Virtual Element Method
etaN = zeros(maxIt,1);
for k = 1:maxIt
    % Step 1: SOLVE
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    [uh,ph,info] = Stokes_mixedVEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    figure(1);
    showresult(node,elem,pde.uexact,uh);
    %[uhI,nodeI,elemI] = EllipticProjection(node,elem,uh,info,2);
    %showresult(nodeI,elemI,pde.uexact,uhI);
    %[phI,nodeI,elemI] = PolyProjection(node,elem,ph,1);
    %showresult(nodeI,elemI,pde.pexact,phI);
    drawnow; %pause(0.1);
    
    % compute errors in discrete L2 and H1 norms
    kOrder = 2; pOrder = 1;
    ErruL2(k) = getL2error(node,elem,uh,info,pde,kOrder);
    ErruH1(k) = getH1error(node,elem,uh,info,pde,kOrder);
    ErrpL2(k) = getL2error_Poly(node,elem,ph,pde,pOrder);
    
    % Step 2: ESTIMATE
    eta = Stokes_mixedVEM_indicator(node,elem,uh,ph,info,pde);
    etaN(k) = norm(eta);
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    
    % Step 4: REFINE
    [node,elem] = PolyMeshRefine1(node,elem,elemMarked);
    
    if (size(node,1)>maxN) || (k==maxIt)
        bdStruct = setboundary(node,elem);
        uh = Stokes_mixedVEM(node,elem,pde,bdStruct);
        step = k
        break;
    end
end

figure,
plot((1:step),etaN(1:step),'k.-','linewidth',1);
xlabel('k'); ylabel('\eta (u_h)');

figure,
id = 1;
h = 1./sqrt(N(id:step));
showrateh(h,etaN(id:step),'r-*','\eta (u_h)', ...
    ErruH1(id:step), 'b-s','|u-u_h|_1', ...
    ErruL2(id:step), 'g-o','||u-u_h|', ...
    ErrpL2(id:step), 'k-d','||p-p_h||')

