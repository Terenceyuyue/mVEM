
%Meshfun: generate basic data structure (node,elem) representing meshes
% Copyright (C) Terence Yu.

clc;clear;close all

% ---------------------- Choice of the domain ------------------
% Options: Lshape_Domain, Rectangle_Domain, Circle_Domain, Upper_Circle_Domain,
% Upper_Circle_Circle_Domain, Rectangle_Circle_Domain, Circle_Circle_Domain
% Michell_Domain, Horn_Domain, Suspension_Domain, Wrench_Domain
Domain = @Rectangle_Domain;  MaxIter = 2000;
Nx = 64; Ny = Nx;
NT = Nx*Ny;
% NT = 16;  
[node,elem] = PolyMesher(Domain,MaxIter,NT);
% Nx = 16; Ny = Nx;
% [node,elem] = PolyMesher(Domain,MaxIter,Nx,Ny);

% ----------------- Visulization of the mesh ------------------
showmesh(node,elem);
% findnode(node); findelem(node,elem); findedge(node,elem);
% axis off

% % % ------------------- Save the mesh data --------------------
meshname = sprintf('meshdata_%d_%d',Nx, Ny);
save(meshname,'node','elem');
