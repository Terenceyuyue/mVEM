clc;clear;close all;

tic;

%% Parameters
maxN = 1e4;  maxIt = 20;
N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
[ErrH1,ErrH2] = deal(zeros(maxIt,1));

%% Initial mesh
load meshdata32

%% PDE data
id = 2; % id = 1,2
pde = biharmonicdata_avem(id);
pde.a = .1;
theta = 0.4*(id==1) + 0.6*(id==2);    

%% Adaptive Virtual Element Method
etaN = zeros(maxIt,1);
for k = 1:maxIt
    % Step 1: SOLVE
    fprintf('Step %d: \n', k);
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    [uh,info] = IPVEM_Zhao(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  
    figure(1); 
    showresult(node,elem,pde.uexact,uh);
    drawnow; %pause(0.1);
    
    % compute errors in discrete norms
    kOrder = 2; 
    infoH1 = info; infoH1.Ph = info.Ph(:,1);
    infoH2 = info; infoH2.Ph = info.Ph(:,2);
    ErrH1(k) = getH1error(node,elem,uh,infoH1,pde,kOrder);
    ErrH2(k) = getH2error(node,elem,uh,infoH2,pde,kOrder);
    
    % Step 2: ESTIMATE
    eta = IPVEM_Zhao_indicator(node,elem,uh,info,pde,bdStruct);
    etaN(k) = norm(eta);
    
    % Step 3: MARK
    elemMarked = mark(elem,eta,theta);
    
    % Step 4: REFINE
    [node,elem] = PolyMeshRefine(node,elem,elemMarked);
    
    if (size(node,1)>maxN) || (k==maxIt)
        bdStruct = setboundary(node,elem);
        uh = IPVEM(node,elem,pde,bdStruct);
        step = k
        break;
    end
end

id = 15;
figure,
plot((1:step),etaN(1:step),'k.-','linewidth',1);
xlabel('k'); ylabel('\eta (u_h)');

figure,
h = 1./sqrt(N(id:step));
showrateh(h,etaN(id:step),'r-*','\eta (u_h)', ErrH2(id:step), 'b-s','|u-u_h|_2')

figure,
h = 1./sqrt(N(id:step));
showrateh(h,h.*etaN(id:step),'r-*','\eta (u_h)', ErrH1(id:step), 'b-s','|u-u_h|_1')


toc
