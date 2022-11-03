clc;clear;close all

%% Example 1
load meshex1
z = linspace(0,1,4);
[node3,elem3] = PolyMesh3Simple(node,elem,z);
figure,
option.FaceAlpha = 0.4;
showmesh(node3,elem3,option);
findnode(node3,'noindex');

% %% Example 2
% g = [2  2  2  2  2  2
%     0  1  1 -1 -1  0
%     1  1 -1 -1  0  0
%     0  0  1  1 -1 -1
%     0  1  1 -1 -1  0
%     1  1  1  1  1  1
%     0  0  0  0  0  0];
% [node,~,elem]=initmesh(g,'hmax',1);
% node = node'; elem = elem(1:3,:)';
% figure,
% subplot(2,2,1)
% showmesh(node,elem);
% subplot(2,2,2)
% [node,elem] = dualMesh(node,elem);
% showmesh(node,elem);
% subplot(2,2,4)
% option.FaceAlpha = 0.5;
% z = linspace(0,1,4);
% [node3,elem3] = PolyMesh3Simple(node,elem,z);
% showmesh(node3,elem3,option); 