function [node,elem] = dualMesh(pp,tt,isDelete)
% Generate polygonal mesh by establishing the dual mesh of a given
% triangulation.
%
% Example-1:
% g = [2  2  2  2  2  2
%     0  1  1 -1 -1  0
%     1  1 -1 -1  0  0
%     0  0  1  1 -1 -1
%     0  1  1 -1 -1  0
%     1  1  1  1  1  1
%     0  0  0  0  0  0];
% [pp,~,tt] = initmesh(g,'hmax',0.5); 
% pp = pp'; tt = tt(1:3,:)';
% [node,elem] = dualMesh(pp,tt);

% Example-2:
% load meshAirfoil % load meshLake
% [node,elem] = dualMesh(node,elem);
%
% Copyright (C) Terence Yu

% numbers
nt = size(tt,1);  np = size(pp,1);  

%% Auxiliary data structure
% totalEdge
totalEdge = sort([tt(:,[2,3]); tt(:,[3,1]); tt(:,[1,2])],2);
[~, i1, totalJ] = unique(totalEdge,'rows'); % first occurence
% edge, bdEdge
[i,j,rep] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
bdEdge = edge(rep==1,:);
% edge2elem
totalJelem = repmat((1:nt)',3,1);
i2(totalJ) = 1:length(totalJ); i2 = i2(:); % second occurence by overlapping
edge2elem = totalJelem([i1,i2]);
% index for boundary nodes
ne = size(edge,1);
isBoundary = false(np,1);
isBoundary(bdEdge(:)) = true;

%% Dual mesh
elem = cell(np,1);
for ip = 1:np
    % index of edges sharing vertex as an endpoint
    v2e = unique([find(edge(:,1)==ip); find(edge(:,2)==ip)]);
    Nv = length(v2e);
    % neighboring triangles form the edge of voronoi cell
    ee = edge2elem(v2e,:);  isCon = false(Nv,1);
    if isBoundary(ip)  % boundary nodes                      
        isCon(rep(v2e)==1) = true;         
        ee(isCon,2) = v2e(isCon)+nt; % connection number
        ee((1:2)+Nv,:) = [ip+nt+ne+[0;0], ee(isCon,2)]; % add edges on the boundary 
        ee = ee([(1:2)+Nv, 1:Nv], :);  % vertex as the first one
        Nv = size(ee,1);
    end
    % connectivity of voronoi cell
    index = zeros(1,Nv);
    index(1:2) = ee(1,:);
    ee(1,:) = []; % delete the old edges
    for i = 2:Nv-1
        current = index(i); % the current vertex
        for s = 1:size(ee,1)
            if ismember(current,ee(s,:))
                next = setdiff(ee(s,:),current);
                ee(s,:) = [];
                break;
            end
        end
        index(i+1) = next;
    end
    elem{ip} = index;
end

%% Nodes in the dual mesh
% centroid of triangles
v1 = tt(:,1); v2 = tt(:,2); v3 = tt(:,3);
nodeCd = 1/3*(pp(v1,:)+pp(v2,:)+pp(v3,:));
% mid-points of boundary edges
z1 = pp(bdEdge(:,1),:); z2 = pp(bdEdge(:,2),:);
nodeBdm = (z1 + z2)/2;
% boundary vertices
nodeBd = pp(unique(bdEdge),:);
% dual nodes
% node = [nodeCd; nodeBdm; nodeBd]; % the original order
node = [nodeBd; nodeBdm; nodeCd]; 

%% Update elem 
[~,~,totalid] = unique(horzcat(elem{:})');
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';

%% Reorder
nBd = size(nodeBd,1); nBdm = size(nodeBdm,1); nCd = size(nodeCd,1);
connection = [(1:nCd)+nBd+nBdm, (1:nBdm)+nBd, (1:nBd)];
elem = cellfun(@(index) connection(index), elem, 'UniformOutput', false);

%% Counterclockwise order
[node,elem] = PolyMesher_Reorder(node,elem);
% PolyMesher_Reorder.m uses the first three vertices to decide the sign
% For voronoi cell corresponding to the vertices on the boundary, the first
% three vertices are: the vertex, the mid-point and the centroid.  Hence,
% the formed triangle is inside the domain, even if the vertex is near the 
% non-convex corner.

%% Delete hanging nodes and non-convex cells
if exist('isDelete','var') && (isDelete==true)
    nBd = sum(isBoundary);
    elemBd = cell(nBd,1);
    s = 1;
    for ip = 1:np
        if isBoundary(ip)
            index = elem{ip};  Nv = length(index);
            elem{ip} = index(1:3);
            elemBd{s} = [index(1)*ones(Nv-3,1), index(3:Nv-1)', index(4:Nv)'];
            s = s + 1;
        end
    end
    elemBd = cell2mat(elemBd);
    elemBd = mat2cell(elemBd, ones(size(elemBd,1),1), 3);
    elem = [elem; elemBd];
end