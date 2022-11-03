clc; close all;
clear variables;

%% Parameters
nameV = [200, 400, 600, 800, 1000]./2;
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);
ErrH2 = zeros(maxIt,1);
Erreps = zeros(maxIt,1); 
ErrI = zeros(maxIt,1);

%% PDE data
epsilon = 10^(-5);
pde = Fourth_order_Singular_Pertubation_Data(epsilon); % no boundary layer
% epsilon = 10^(-0);
% pde = Fourth_order_Singular_Pertubation_Data1(epsilon); % with boundary layer

%% Exact norms
% exact
Du = @(x,y) sum(pde.Du([x,y]).^2, 2); 
uH1square = squareint(Du,[0,1,0,1]); % square
DDu = @(x,y) sum(pde.DDu([x,y]).^2, 2); 
uH2square = squareint(DDu,[0,1,0,1]);

%[node,elem] = squaremesh([0 1 0 1], 0.5,0.5);
%% Virtual element method
for k = 1:maxIt
    % load mesh
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %[node,elem] = uniformrefine(node,elem);
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    [w,info] = Fourth_order_Singular_Perturbation_C0VEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(w);  h(k) = 1/sqrt(size(elem,1));
    if size(node,1)<2e3
        figure(1);
        showresult(node,elem,pde.uexact,w);
        drawnow; %pause(0.1);
    end
    % compute errors in discrete L2,H1,H2 norms
    kOrder = 2;
    ErrL2(k) = getL2error(node,elem,w,info,pde,kOrder);
    ErrH1(k) = getH1error(node,elem,w,info,pde,kOrder);
    ErrH2(k) = getH2error(node,elem,w,info,pde,kOrder);
    Erreps(k) = sqrt((ErrH1(k)^2 + epsilon^2*ErrH2(k)^2)/(uH1square+epsilon^2*uH2square));
    % compute errors in discrete energy norm
    chi = w;
    chie = getDof_C0VEM(node,elem,pde);
    kk = info.kk; DofI = info.DofI;
    kk = kk(DofI,DofI); chie = chie(DofI); chi = chi(DofI);
    ErrI(k) = sqrt(abs((chi-chie)'*kk*(chi-chie)));
end

%% Plot convergence rates and display error table
% figure;
% subplot(1,2,1); showrateErr(h,ErrL2,ErrH1);
% subplot(1,2,2); showrateErr(h,ErrH2,ErrI);
figure;
showrate(h,Erreps,'b-s','k--','|||u-uh|||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1','|u-u_h|_2','|||u-uh|||'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e',ErrH2,'%0.5e',Erreps,'%0.5e');



