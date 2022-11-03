function [node,elem,pp,tt] = EdgePolyMesher(P1,P2,options)
%EdgePolyMesher Generate polygonal meshes according to the line segments for
% the boundary of a polygon.
% The line segments must be in the conventional directions.
%
%    [node,elem] = PolyEdgeMesh(P1,P2)
%    P1 = [x,y] is the starting points of edges
%    P2 = [x,y] is the ending points of edges
%    
%    options: plot is true or false
%             refineTimes is the number of refinements
%

% Copyright (C)  Terence Yu


if nargin==2
    options.plot = false; 
end

%% Decomposed geometry matrix
N = size(P1,1);
g = zeros(7,N);
g(1,:) = 2;  % 2 for line 
g(2,:) = P1(:,1); % x: starting points
g(3,:) = P2(:,1); % x: ending points
g(4,:) = P1(:,2); % y: starting points
g(5,:) = P2(:,2); % y: ending points
g(6,:) = 1; % subdomain on the left
g(7,:) = 0; % subdomain on the right

%% Triangulation
if ~isfield(options,'hmax')
    options.hmax = 1;
end
[p,e,t] = initmesh(g,'hmax',options.hmax,'Hgrad',1.1); % initial mesh
if ~isfield(options,'refineTimes')
    options.refineTimes = 1; 
end
for i = 1:options.refineTimes
    [p,e,t] = refinemesh(g,p,e,t);
end
pp = p';
tt = t(1:3,:)';

%% Dual mesh
[node,elem] = dualMesh(pp,tt);
[node,elem] = PolyMesher_Reorder(node,elem);

%% Plot mesh
if options.plot
    figure,
    subplot(1,2,1),
    showmesh(pp,tt); 
    subplot(1,2,2),
    showmesh(node,elem); 
end