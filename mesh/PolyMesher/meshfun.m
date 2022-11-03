
%Meshfun: generate basic data structure (node,elem) representing meshes
% Copyright (C) Terence Yu.

clc;clear;close all

% ---------------------- Choice of the domain ------------------
% Options: Lshape_Domain, Rectangle_Domain, Circle_Domain, Upper_Circle_Domain,
% Upper_Circle_Circle_Domain, Rectangle_Circle_Domain, Circle_Circle_Domain
% Michell_Domain, Horn_Domain, Suspension_Domain, Wrench_Domain
Domain = @Rectangle_Domain;  MaxIter = 400;
NT = 10000;  
[node,elem] = PolyMesher(Domain,MaxIter,NT);
% Nx = 10; Ny = 10;
% [node,elem] = PolyMesher(Domain,MaxIter,Nx,Ny);

% ----------------- Visulization of the mesh ------------------
showmesh(node,elem);
% findnode(node); findelem(node,elem); findedge(node,elem);
% axis off

% % % ------------------- Save the mesh data --------------------
meshname = sprintf('meshdata%d',NT);
save(meshname,'node','elem');
