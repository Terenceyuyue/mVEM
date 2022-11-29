clc;clear;close all

% Example-1
g = [2  2  2  2  2  2
    0  1  1 -1 -1  0
    1  1 -1 -1  0  0
    0  0  1  1 -1 -1
    0  1  1 -1 -1  0
    1  1  1  1  1  1
    0  0  0  0  0  0];
[pp,~,tt] = initmesh(g,'hmax',1); 
pp = pp'; tt = tt(1:3,:)';
figure, 
subplot(1,2,1)
showmesh(pp,tt);
subplot(1,2,2)
[node,elem] = dualMesh(pp,tt);
showmesh(node,elem);
findnode(node); findelem(node,elem)

% Delete hanging nodes and non-convex boundary elements by triangulation
isDelete = 1;
[node,elem] = dualMesh(pp,tt,isDelete);
figure, showmesh(node,elem);

% % Example-2
% load meshAirfoil % load meshLake
% figure,
% subplot(1,2,1)
% showmesh(node,elem);
% subplot(1,2,2)
% [node,elem] = dualMesh(node,elem);
% showmesh(node,elem);
% 
% % Example-3
% load meshLake
% figure,
% subplot(1,2,1)
% showmesh(node,elem);
% subplot(1,2,2)
% [node,elem] = dualMesh(node,elem);
% showmesh(node,elem);