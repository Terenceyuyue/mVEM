
%Meshfun: genearte basic data structure (node, elem) representing meshes
% Copyright (C) Terence Yue Yu.

% ----------- Useage -------------
% Preserved as meshdata.mat
% load('meshdata.mat') to load basic data structure: node, elem
% ----------------------------------

clc;clear;close all
addpath(genpath(pwd)); % Include all folders under the current path

% -------- Choice of the domain --------
% Options: Rectangle_Domain, Circle_Domain, Upper_Circle_Domain,
% Upper_Circle_Circle_Domain, Rectangle_Circle_Domain, Circle_Circle_Domain
Domain = @Circle_Circle_Domain; 
NT = 10; MaxIter = 50;
[node,elem]= PolyMesher(Domain,NT,MaxIter);

% -------- Plot the mesh --------
showmesh(node,elem);
findnode(node);

% -------- Save the mesh data --------
save meshdata node elem
