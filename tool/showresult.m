function showresult(node,elem,u,uh)
%showresult displays the mesh and the solution
%
% Copyright (C) Terence Yu.

clf;  % clear figures

set(gcf,'Units','normal');
set(gcf,'Position',[0.25,0.25,0.7,0.25]);

%% Plot mesh
subplot(1,3,1);
showmesh(node,elem); 
hold on;
plot(node(:,1),node(:,2),'k.', 'MarkerSize', 4);

%% Plot exact solution
subplot(1,3,2);
ue = u(node); ue = ue(:,1);
showsolution(node,elem,ue(1:size(node,1)));

%% Plot numerical solution
subplot(1,3,3);
showsolution(node,elem,uh(1:size(node,1),1));