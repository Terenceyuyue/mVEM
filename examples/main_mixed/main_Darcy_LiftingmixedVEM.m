clc; close all;
clear variables;

%% Parameters
nameV = [32, 64, 128, 256, 512];
maxIt = length(nameV);
h = zeros(maxIt,1);   N = zeros(maxIt,1);
ErruL2 = zeros(maxIt,1);
ErrpL2 = zeros(maxIt,1);

%% PDE data
pde = Darcydata;

%% Virtual element method
for k = 1:maxIt
    % load mesh
    fprintf('Mesh %d: \n', k);
    load( ['meshdata', num2str(nameV(k)), '.mat'] );
    % get boundary information
    bdStruct = setboundary(node,elem);
    % solve the problem
    [uh,ph,info] = Darcy_LiftingmixedVEM(node,elem,pde,bdStruct);
    % record and plot
    N(k) = length(uh);  h(k) = 1/sqrt(size(elem,1));
    %[uhI,phI,nodeI,elemI] = ProjectionDarcy(node,elem,uh,ph,info,pde);
    figure(1); 
    %showresult(nodeI,elemI,pde.uexact,uhI);
    [phI,nodeI,elemI] = PolyProjection(node,elem,ph,1);
    showresult(nodeI,elemI,pde.pexact,phI);
    drawnow; % pause(0.1)
    % compute errors in discrete L2 norm
    [ErruL2(k),ErrpL2(k)] = getL2error_Darcy(node,elem,uh,ph,info,pde);
end

%% Plot convergence rates and display error table
figure,
ErrpL22 = 5*ErrpL2; % widen the gap
showrateh(h,ErruL2,'r-*','||u-u_h||',  ErrpL22, 'b-s','||p-p_h||')

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||','||p-p_h||'};
disptable(colname,N,[],h,'%0.3e',ErruL2,'%0.5e',ErrpL2,'%0.5e');

%% Conclusion
%
% The optimal rate of convergence of the L2-norm for u (2nd order) and p
% (2st order) is observed for k = 1.

% %% Display exact d.o.f.s
% % auxiliary data
% aux = auxgeometry(node,elem);
% centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
% auxT = auxstructure(node,elem);
% edge = auxT.edge; NE = size(edge,1);
% u = pde.uexact;  p = pde.pexact;
% % ue
% z1 = node(edge(:,1),:); z2 = node(edge(:,2),:);  ze = (z1+z2)./2;
% e = z1-z2;  % e = z2-z1
% Ne = [-e(:,2),e(:,1)];
% ue1 = 1/6*sum((u(z1)+4*u(ze)+u(z2)).*Ne,2);
% ue2 = 1/12*sum((u(z2)-u(z1)).*Ne,2);
% ue = [ue1;ue2];
% bdEdgeIdx = bdStruct.bdEdgeIdx;
% id = [bdEdgeIdx; bdEdgeIdx+NE];
% ue(id) = uh(id); % counterclockwise
% 
% % pe
% NT = size(elem,1);
% pe = zeros(NT,3);
% for iel = 1:NT
%     index = elem{iel};   Nv = length(index);
%     xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
%     nodeT = [node(index,:);centroid(iel,:)];
%     elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
%     m1 = @(x,y) 1+0*x;
%     m2 = @(x,y) (x-xK)./hK;
%     m3 = @(x,y) (y-yK)./hK;
%     m = {m1,m2,m3};
%     mvec = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
%     H = zeros(3,3);  d = zeros(3,1);
%     for i = 1:3
%         pm = @(x,y) p([x,y]).*m{i}(x,y);
%         d(i) = integralTri(pm,4,nodeT,elemT); 
%         for j = i:3
%             mm = @(x,y) m{i}(x,y).*m{j}(x,y);
%             H(i,j) = integralTri(mm,4,nodeT,elemT);
%             H(j,i) = H(i,j);
%         end
%     end
%     pe(iel,:) = H\d;
% end
% pe = pe(:);
% 
% % plot
% pid = 1:fix(2*NE/100):2*NE;
% figure, plot(pid, uh(pid),'-or', pid, ue(pid),'-*k','linewidth',1);
% pid = 1:fix(3*NT/100):3*NT;
% figure, plot(pid,ph(pid),'-or', pid,pe(pid),'-*k','linewidth',1);