
%Meshfun: generate basic data structure (node,elem) representing meshes
% Copyright (C) Terence Yue Yu.

clc;clear;close all

% ---------------------- Choice of the domain ------------------
% Options: Rectangle_Domain, Circle_Domain, Upper_Circle_Domain,
% Upper_Circle_Circle_Domain, Rectangle_Circle_Domain, Circle_Circle_Domain
% Michell_Domain, Horn_Domain, Suspension_Domain, Wrench_Domain
Domain = @Rectangle_Domain;  MaxIter = 500;
NT = 200;  
% Nx = 5; Ny = 5;
[node,elem] = PolyMesher(Domain,MaxIter,NT);

% ----------------- Visulization of the mesh ------------------
showmesh(node,elem);
% findnode(node); findelem(node,elem); findedge(node,elem);
% axis off

% % % ------------------- Save the mesh data --------------------
% save meshdata512 node elem
