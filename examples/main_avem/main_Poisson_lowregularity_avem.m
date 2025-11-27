clc;clear;close all;


%% Parameters
maxN = 2e4;     theta = 0.6;    maxIt = 100;
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);

%% Generate an initial mesh
load meshLshape1

%% Get the PDE data
pde = Poissdata_lowregularity();

%% Adaptive Virtual Element Method
etaN = zeros(maxIt,1);
for k = 1:maxIt
    % Step 1: SOLVE
    fprintf('Step %d: \n', k);
    % get boundary information
    [bdStruct,node,elem] = setboundary(node,elem);
    % solve
    [uh,info] = PoissonVEM(node,elem,pde,bdStruct);
    h = min(info.diameter)
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    if N(k)<2e3
        figure(1);
        showresult(node,elem,pde.uexact,uh);
        drawnow; %pause(0.1);
    end
    
    % compute errors in discrete L2 and H1 norms
    kOrder = 1;
    % ErrL2(k) = getL2error(node,elem,uh,info,pde,kOrder);
    ErrH1(k) = getH1error(node,elem,uh,info,pde,kOrder);
    
    % Step 2: ESTIMATE
    eta = PoissonVEM_indicator(node,elem,uh,info,pde);
    etaN(k) = norm(eta);
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    
    % Step 4: REFINE
    [node,elem] = PolyMeshRefine1(node,elem,elemMarked);
    
    if (size(node,1)>maxN) || (k==maxIt)
        bdStruct = setboundary(node,elem);
        uh = PoissonVEM(node,elem,pde,bdStruct);
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
showrateh(h,etaN(id:step),'r-*','\eta (u_h)', ErrH1(id:step), 'b-s','|u-u_h|_1');