
%Meshfun3D: generate basic data structure (node,elem) representing meshes
% Copyright (C) Terence Yue Yu.

% --------------------------- Usage ---------------------------
% Preserved as meshdata.mat
% load('meshdata.mat') to load basic data structure: node, elem
% --------------------------------------------------------------

clc;clear;close all

% ---------------------- Choice of the domain ------------------
% Options: Cube_Domain (NT>=2), Sphere_Domain (NT>=40)
Domain = @Sphere_Domain; 
MaxIter = 200;
NT = 50;  rm = false; % rm = true 仅去除极短的边
[node,elem]= PolyMesher3D(Domain,MaxIter,NT,rm);

% ----------------- Visulization of the mesh ------------------
option.FaceAlpha = 1;
figure, showmesh(node,elem,option); axis equal

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
