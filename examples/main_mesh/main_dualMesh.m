clc;clear;close all

% Example-1
g = [2  2  2  2  2  2
    0  1  1 -1 -1  0
    1  1 -1 -1  0  0
    0  0  1  1 -1 -1
    0  1  1 -1 -1  0
    1  1  1  1  1  1
    0  0  0  0  0  0];
[node,~,elem]=initmesh(g,'hmax',1); 
node = node'; elem = elem(1:3,:)';
figure, 
subplot(1,2,1)
showmesh(node,elem);
subplot(1,2,2)
[node,elem] = dualMesh(node,elem);
showmesh(node,elem);

% Example-2
load meshAirfoil % load meshLake
figure,
subplot(1,2,1)
showmesh(node,elem);
subplot(1,2,2)
[node,elem] = dualMesh(node,elem);
showmesh(node,elem);

% Example-3
load meshLake
figure,
subplot(1,2,1)
showmesh(node,elem);
subplot(1,2,2)
[node,elem] = dualMesh(node,elem);
showmesh(node,elem);