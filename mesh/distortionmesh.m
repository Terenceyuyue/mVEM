clc;clear;close all

[node,elem] = squaremesh([0 1 0 1],1/20,1/20,'rec');

tc = 0.1;
xmap = sin(2*pi*node(:,1)).*sin(2*pi*node(:,2));
ymap = xmap;
node = node + tc*[xmap, ymap];

showmesh(node,elem);