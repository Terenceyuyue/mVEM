clc; close all;
clear variables;

%% Parameters
nameV = [32, 64, 128, 256, 512]';
maxIt = length(nameV);
h = zeros(maxIt,1);  N = zeros(maxIt,1);
ErrL2 = zeros(maxIt,1);
ErrH1 = zeros(maxIt,1);
ErrH2 = zeros(maxIt,1);
Erreps = zeros(maxIt,1); 
ErrI = zeros(maxIt,1);

%% PDE data
iota = 1e-0; lambda = 1e0; mu = 1; 
para = struct('lambda',lambda, 'mu',mu, 'iota',iota);
pde = Strain_gradient_elasticity_data2(para);

%% Exact norms
% exact
tri = eye(4);
u1x = @(x,y) pde.Du([x,y]).^2*tri(:,1);
u1y = @(x,y) pde.Du([x,y]).^2*tri(:,2);
u2x = @(x,y) pde.Du([x,y]).^2*tri(:,3);
u2y = @(x,y) pde.Du([x,y]).^2*tri(:,4);
Du = @(x,y) u1x(x,y)+u1y(x,y)+u2x(x,y)+u2y(x,y);
uH1square = squareint(Du,[0,1,0,1]); % square

tri = eye(8);
u1xx = @(x,y) pde.DDu([x,y]).^2*tri(:,1);
u1xy = @(x,y) pde.DDu([x,y]).^2*tri(:,2);
u1yx = @(x,y) pde.DDu([x,y]).^2*tri(:,3);
u1yy = @(x,y) pde.DDu([x,y]).^2*tri(:,4);
u2xx = @(x,y) pde.DDu([x,y]).^2*tri(:,5);
u2xy = @(x,y) pde.DDu([x,y]).^2*tri(:,6);
u2yx = @(x,y) pde.DDu([x,y]).^2*tri(:,7);
u2yy = @(x,y) pde.DDu([x,y]).^2*tri(:,8);
DDu = @(x,y) u1xx(x,y)+u1xy(x,y)+u1yx(x,y)+u1yy(x,y)+u2xx(x,y)+u2xy(x,y)+u2yx(x,y)+u2yy(x,y);
uH2square = squareint(DDu,[0,1,0,1]);

%[node,elem] = squaremesh([0 1 0 1], 0.5, 0.5);
%% Virtual element method
for k = 1:maxIt
    % load mesh
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    %[node,elem] = uniformrefine(node,elem);
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve the problem
    [w,info] = Strain_gradient_elasticity_C0VEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(w);  h(k) = 1/sqrt(size(elem,1));
    if size(node,1)<2e3
        figure(1);
        showresult(node,elem,pde.uexact,w);
        pause(0.1);
    end
    % compute errors in discrete H1 and H2 norms
    kOrder = 2;
    Ph = info.Ph;
    Ph1 = Ph(:,2);  info.Ph = Ph1;
    ErrL2(k) = getL2error_vector(node,elem,w,info,pde,kOrder);
    Ph1 = Ph(:,1);  info.Ph = Ph1;
    ErrH1(k) = getH1error_vector(node,elem,w,info,pde,kOrder);
    Ph2 = Ph(:,2);  info.Ph = Ph2;
    ErrH2(k) = getH2error_vector(node,elem,w,info,pde,kOrder);
    Erreps(k) = sqrt((ErrH1(k)^2 + iota^2*ErrH2(k)^2)/(uH1square+iota^2*uH2square));
    
    % compute errors in discrete energy norm
    % auxstructure
    auxT = auxstructure(node,elem);
    edge = auxT.edge;
    % chi1: evaluation at all vertices
    uexact = pde.uexact;
    chi1 = uexact(node);
    % chi2: evaluation at all mid-edge vertices
    z1 = node(edge(:,1),:); z2 = node(edge(:,2),:); zc = (z1+z2)/2;
    chi2 = uexact(zc);
    % chi3: moments of \partial_n v on edges with given orientation
    e = z1-z2;  % e = z2-z1
    Ne = [-e(:,2),e(:,1)];  % scaled ne
    Du = pde.Du;
    Dw1 = Du(z1); Dwc = Du(zc); Dw2 = Du(z2); % Dw = [u1x,u1y,u2x,u2y]
    wnD1 = sum(1/6*(Dw1(:,1:2)+4*Dwc(:,1:2)+Dw2(:,1:2)).*Ne,2);
    wnD2 = sum(1/6*(Dw1(:,3:4)+4*Dwc(:,3:4)+Dw2(:,3:4)).*Ne,2);
    wnD = [wnD1,wnD2];
    chi3 = wnD;
    % chi
    chi = w;
    chie = [chi1; chi2; chi3]; chie = chie(:);
    % Relative error in terms of discrete energy norm ------ direction ???
    kk = info.kk; freeDof = info.freeDof;
    kk = kk(freeDof,freeDof); chi = chi(freeDof); chie = chie(freeDof);
    ErrI(k) = sqrt(abs((chi-chie)'*kk*(chi-chie)));
end

%% Plot convergence rates and display error table
figure;
showrate(h,Erreps,'b-s','k--','|||u-uh|||');
% figure;
% showrate(h,ErrI,'b-s','k--','||u-uh||_A');
% figure;
% showrate(h,ErrL2,'b-s','k--','||u-uh||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','|u-u_h|_1','|u-u_h|_2','|||u-uh|||','||uI-uh||'};
disptable(colname,nameV,[],h,'%0.3e',ErrL2,'%0.5e',ErrH1,'%0.5e',ErrH2,'%0.5e',Erreps,'%0.5e',ErrI,'%0.5e');

format shorte
Erreps'
% 
% figure,spy(info.kk)
