%mpt_meshfun2D: generate basic data structure (node,elem) representing meshes
% Copyright (C) Terence Yu.

% --------------------------- Usage ---------------------------
% Preserved as meshdata.mat
% load('meshdata.mat') to load basic data structure: node, elem
% --------------------------------------------------------------

clc;clear;close all

% ---------------------- Choice of the domain ------------------
NT = 10;  MaxIter = 100;
BdBox = [0 1 0 1];
Pb = mpt_rectangle_domain(BdBox); 
[node,elem] = MPT_Mesher2D(Pb,NT,MaxIter,BdBox);

% ----------------- Visulization of the mesh ------------------
showmesh(node,elem);
%findnode(node); findelem(node,elem)

% % ------------------- Save the mesh data --------------------
%save meshdata node elem