%mpt_meshfun3D: generate basic data structure (node,elem) representing meshes
% Copyright (C) Terence Yue Yu.

% --------------------------- Usage ---------------------------
% Preserved as meshdata.mat
% load('meshdata.mat') to load basic data structure: node, elem
% --------------------------------------------------------------

clc;clear;close all
tic;
% ---------------------- Choice of the domain ------------------
NT = 10;  MaxIter = 10;
BdBox = [0 1 0 1 0 1];
Pb = mpt_rectangle_domain(BdBox); 
[node,elem,Pn,P] = MPT_Mesher3D(Pb,NT,MaxIter,BdBox);

% ----------------- Visulization of the mesh ------------------
option.FaceAlpha = 1;
figure, showmesh(node,elem,option); axis equal
% hold on
% center = P; range = (1:size(center,1))';
% plot3(center(:,1),center(:,2),center(:,3),'o','LineWidth',1,'MarkerEdgeColor','k',...
%     'MarkerFaceColor','y','MarkerSize',18);
% text(center(:,1)-0.04,center(:,2),center(:,3),int2str(range),'FontSize',12,...
%      'FontWeight','bold','Color','k');
% % findnode(node)

iel = 3; elemf = elem{iel}; 
option.FaceAlpha = 0.25;
figure,showmesh(node,elemf,option); axis equal
range = unique(horzcat(elemf{:}))'; 
findnode(node,range)

toc

% % ------------------- Save the mesh data --------------------
% save meshdata node elem