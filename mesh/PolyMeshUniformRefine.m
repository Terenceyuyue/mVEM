function [node,elem] = PolyMeshUniformRefine(node,elem,refineType)
%PolyMeshRefine uniformly refines a 2-D polygonal mesh
%
% We divide elements by connecting the midpoint of each edge to its
% barycenter.
%
% Copyright (C) Terence Yu.

tol = 1e-10; % accuracy for finding midpoint

if nargin==2, refineType = 1; end

%% Get auxiliary data
if ~iscell(elem), elem = num2cell(elem,2); end
% diameter
diameter = cellfun(@(index) max(pdist(node(index,:))), elem);
if max(diameter)<tol
    disp('The mesh is too dense'); return;
end
% totalEdge
shiftfun = @(verts) [verts(2:end),verts(1)];
T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
v0 = horzcat(elem{:})'; v1 = horzcat(T1{:})';
totalEdge = sort([v0,v1],2);
% edge, elem2edge
[edge, ~, totalJ] = unique(totalEdge,'rows');
elemLen = cellfun('length',elem);
elem2edge = mat2cell(totalJ',1,elemLen)';
% number
N = size(node,1); NT = size(elem,1); NE = size(edge,1);
elemLen = cellfun('length',elem); 

%% Update node
nodeEdge = (node(edge(:,1),:) + node(edge(:,2),:))/2;
nodeCenter = zeros(NT,2);
for iel = 1:NT
    index = elem{iel}; verts = node(index,:);
    nodeCenter(iel,:) = polycentroid(verts);
end
switch refineType
    case 1
        node = [node; nodeEdge; nodeCenter];
    case 2
        node = [node; nodeEdge];
    case 3
        node = [node; nodeEdge];
end

%% Update elem
elemNew = cell(NT,1);
if refineType == 1  % center -- midpoints
    s = 1;
    for iel = 1:NT
        index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
        ide = indexEdge+N;  % connection number
        z1 = ide([Nv,1:Nv-1]);   z0 = index;
        z2 = ide;                zc = iel*ones(Nv,1)+N+NE;  % connection number
        elemNew{s} = [z1(:), z0(:), z2(:), zc(:)];
        s = s+1;
    end
    elem = num2cell(vertcat(elemNew{:}), 2);
end

if refineType == 2  % connecting midpoints
    Ntri = sum(elemLen); Npoly = NT;
    elemNew = cell(Ntri+Npoly,1);
    s = 0;
    for iel = 1:NT
        index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
        ide = indexEdge+N;  % connection number
        z1 = ide([Nv,1:Nv-1]);   z0 = index;      z2 = ide;                
        elemTri = [z1(:), z0(:), z2(:)];  % triangle
        elemNew(s+1:s+Nv) = mat2cell(elemTri,ones(1,Nv), 3);
        s = s+Nv;
        elemNew{Ntri+iel} = ide;        
    end
    elem = elemNew;
end

if refineType == 3  % two edges for hanging node 
    for iel = 1:NT
        index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
        z = zeros(1,2*Nv);
        z(1:2:end) = index; z(2:2:end) = indexEdge+N; % connection number
        elem{iel} = z;
    end
end