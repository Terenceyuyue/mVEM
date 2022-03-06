clc;clear;close all

load meshex1
figure,
% 初始剖分
showmesh(node,elem);
findnode(node,'noindex');
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

for i = 1:10
    % 第 3 次加密
    figure(4),
    elemMarked = randperm(15,5);
    [node,elem] = PolyMeshRefine(node,elem,elemMarked);
    showmesh(node,elem);
    findnode(node,'noindex');
    % findelem(node,elem);
    drawnow;
end

