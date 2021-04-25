clc; close all;
clear variables;

%% Parameters
maxN = 1e4;     theta = 0.4;    maxIt = 30;
h = zeros(maxIt,1);  N = zeros(maxIt,1);

%% Generate an initial mesh
load distortionPolygon %Lshape100

%% Get the PDE data
iota = 1e-0; lambda = 1; mu = 1;
para = struct('lambda',lambda, 'mu',mu, 'iota',iota);
pde = Strain_gradient_elasticity_data5(para);

%% Adaptive Virtual Element Method
etaN = zeros(maxIt,1);
for k = 1:maxIt
    % Step 1: SOLVE
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    [w,info] = Strain_gradient_elasticity_C0VEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(w);  h(k) = 1/sqrt(size(elem,1));
    if size(node,1)<2e3
        figure(1);
        showresult(node,elem,pde.uexact,w);
        pause(0.1);
    end
    
    % Step 2: ESTIMATE
    eta = Strain_gradient_elasticity_C0VEM_indicator(node,elem,w,info,pde);
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    etaN(k) = norm(eta); 
    
    % Step 4: REFINE
    [node,elem] = PolyMeshRefine(node,elem,elemMarked);
    
    if (size(node,1)>maxN) || (k==maxIt) || max(eta)<1e-4
        bdStruct = setboundary(node,elem);
        uh = Strain_gradient_elasticity_C0VEM(node,elem,pde,bdStruct);
        step = k
        break;
    end
end

figure,
plot((1:step),etaN(1:step),'k.-','linewidth',1);
% xlabel('k'); ylabel('\eta (u_h)');
