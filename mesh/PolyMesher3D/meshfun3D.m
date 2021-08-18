
%Meshfun3D: generate basic data structure (node,elem) representing meshes
% Copyright (C) Terence Yu.

% --------------------------- Usage ---------------------------
% Preserved as meshdata.mat
% load('meshdata.mat') to load basic data structure: node, elem
% --------------------------------------------------------------

clc;clear;close all

% ---------------------- Choice of the domain ------------------
% Options: Cube_Domain (NT>=2), Sphere_Domain (NT>=40)
Domain = @Sphere_Domain; 
MaxIter = 100;
NT = 100;  
[node,elem]= PolyMesher3D(Domain,MaxIter,NT);
% Nx = 10; Ny = 5; Nz = 10;
% [node,elem]= PolyMesher3DNew(Domain,MaxIter,Nx,Ny,Nz);

% ----------------- Visulization of the mesh ------------------
option.FaceAlpha = 1;
option.facecolor = [0.5 0.9 0.45];
figure, 
showmesh(node,elem,option); axis equal
axis off

% hold on
iel = 1;
elemf = elem{iel};
% option.facecolor = 'yellow';
% showmesh(node,elemf,option), axis equal
% range = unique(horzcat(elemf{:}));
% findnode(node,range)

figure,
showmesh(node,elemf), axis equal
range = unique(horzcat(elemf{:}));
findnode(node,range)

% % ------------------- Save the mesh data --------------------
%save meshdata node elem
