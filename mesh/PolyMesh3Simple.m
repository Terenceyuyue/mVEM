function [node3,elem3] = PolyMesh3Simple(node,elem,z)
% Generate polyheral meshes with repeated top and bottom faces 
% node,elem: the polygonal mesh of the top or the bottom plane
% z: the 1-D grid points along z-axis
% 
%   Example 1:
%      load meshdata1
%      z = linspace(0,1,4);
%      [node3,elem3] = PolyMesh3Simple(node,elem,z);
%      figure,
%      option.FaceAlpha = 0.5;
%      showmesh(node3,elem3,option); 
%
% Copyright (C)  Terence Yu.

if nargin==2
    z = linspace(0,1,4);
end
N = size(node,1); NT = size(elem,1);

%% node3
J = length(z)-1;
node3 = [repmat(node,J+1,1),  kron(z(:),ones(N,1))];

%% elem3
elem3 = cell(NT*J,1);
for level = 1:J
    for iel = 1:NT
        index = elem{iel} + (level-1)*N;  Nv = length(index);
        elemf = cell(Nv+2,1);
        elemf{1} = index;  % bottom
        elemf{2} = index(end:-1:1) + N; % top
        v1 = 1:Nv;  v2 = [2:Nv,1];
        for i = 1:Nv  % side
            elemf{i+2} = index([v2(i),v1(i),v1(i),v2(i)]) + [0,0,N,N];
        end
        elem3{iel + (level-1)*NT} = elemf;
    end
end