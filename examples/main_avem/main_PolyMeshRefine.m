clc;clear;close all

load meshex1
figure,
% 初始剖分
showmesh(node,elem);
findnode(node); 
% findelem(node,elem);

% 第 1 次加密
figure,
elemMarked = [2,5];
[node,elem] = PolyMeshRefine(node,elem,elemMarked);
showmesh(node,elem);
findnode(node,'noindex');
% findelem(node,elem);
%
% 第 2 次加密
figure,
elemMarked = [10];
[node,elem] = PolyMeshRefine(node,elem,elemMarked);
showmesh(node,elem);
findnode(node,'noindex');
% findelem(node,elem);
%
% 第 3 次加密
figure,
elemMarked = [10,15,3,4,8,5,1];
[node,elem] = PolyMeshRefine(node,elem,elemMarked);
showmesh(node,elem);
findnode(node,'noindex');
% findelem(node,elem);

