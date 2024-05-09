clc;clear;close all

%%% This is not an appropriate way. See gunelveMesh.m instead.

n = 4;


% right
P0 = [48, 44];   P1 = [48, 60];   v = P1 - P0;
L = vecnorm(v,2,2);
s = linspace(0,1,fix(L/n))';  s = [s, s];
p1 = repmat(P0, size(s,1), 1) + s.*repmat(v,size(s,1),1);
% upper
P0 = [48, 60];   P1 = [0, 44];   v = P1 - P0;
L = vecnorm(v,2,2);
s = linspace(0,1,fix(L/n))';  s = [s, s];
p2 = repmat(P0, size(s,1), 1) + s.*repmat(v,size(s,1),1);
% left
P0 = [0, 44];   P1 = [0, 0];   v = P1 - P0;
L = vecnorm(v,2,2);
s = linspace(0,1,fix(L/n))';  s = [s, s];
p3 = repmat(P0, size(s,1), 1) + s.*repmat(v,size(s,1),1);
% lower
P0 = [0, 0];   P1 = [48, 44];   v = P1 - P0;
L = vecnorm(v,2,2);
s = linspace(0,1,fix(L/n))';  s = [s, s];
p4 = repmat(P0, size(s,1), 1) + s.*repmat(v,size(s,1),1);

P1 = [p1(1:end-1,:); p2(1:end-1,:); p3(1:end-1,:);  p4(1:end-1,:)]; 
P2 = [p1(2:end,:);   p2(2:end,:);   p3(2:end,:);    p4(2:end,:)];

options.plot = true;
options.refineTimes = 0;
options.hmax = n;
[node,elem] = EdgePolyMesher(P1,P2,options);
findnode(node,'noindex')

%save meshdata_Cook1 node elem  % n = 4, 2, 1
