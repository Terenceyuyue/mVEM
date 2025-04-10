clc; close all;
clear variables;

%% Parameters
nameV = [32, 64, 128, 256, 512];
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);
ErrH2 = zeros(maxIt,1);
ErrEps = zeros(maxIt,1);
ErrI = zeros(maxIt,1);

%% PDE data
id = 2; % id = 1,2,3 (2,3 with boundary layers)
epsilon = 10^(-8);
pde = Fourth_order_Singular_Pertubation_Data(id,epsilon); % no boundary layer

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
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %[node,elem] = uniformrefine(node,elem);
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve
    [w,info] = Fourth_order_Singular_Perturbation_C0VEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(w);  h(k) = 1/sqrt(size(elem,1));
    figure(1);
    showresult(node,elem,pde.uexact,w);
    drawnow; %pause(0.1);
    % compute errors in discrete L2,H1,H2 norms
    kOrder = 2;
    ErrL2(k) = getL2error(node,elem,w,info,pde,kOrder);
    ErrH1(k) = getH1error(node,elem,w,info,pde,kOrder);
    ErrH2(k) = getH2error(node,elem,w,info,pde,kOrder);
    ErrEps(k) = sqrt((ErrH1(k)^2 + epsilon^2*ErrH2(k)^2)/(uH1square+epsilon^2*uH2square));
    %%compute errors in discrete energy norm
    % chi = w;
    % chie = getDof_C0VEM(node,elem,pde);
    % kk = info.kk; DofI = info.DofI;
    % kk = kk(DofI,DofI); chie = chie(DofI); chi = chi(DofI);
    % ErrI(k) = sqrt(abs((chi-chie)'*kk*(chi-chie)));
end

%% Plot convergence rates and display error table
figure;
showrate(h,ErrEps,'b-s','k--','|||u-uh|||');

fprintf('\n')
fprintf('Errors for convergence rate:\n');
fprintf('\n')
fprintf('   &');
fprintf(' %0.4e   &', ErrEps);
fprintf('\n')

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1','|u-u_h|_2','|||u-uh|||'};
disptable(colname,N,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e',ErrH2,'%0.5e',ErrEps,'%0.5e');