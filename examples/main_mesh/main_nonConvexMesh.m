clc;clear;close all

[node,elem] = nonConvexMesh([0 1 0 1], 10, 10);
showmesh(node,elem)
findnode(node,'noindex')