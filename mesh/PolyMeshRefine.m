function [node,elem] = PolyMeshRefine(node,elem,elemMarked)
%PolyMeshRefine refines a 2-D polygonal mesh satisfying one-hanging node rule
%
% We divide elements by connecting the midpoint of each edge to its
% barycenter.
% We remove small edges by further partitioning some adjacent elements.
%
% Copyright (C) Terence Yu.

if nargin==2  % refine all the elements
    [node,elem] = PolyMeshUniformRefine(node,elem);
    return;  
end  

idElemMarked = unique(elemMarked(:)); % in ascending order
tol = 1e-10; % accuracy for finding midpoint

%% Get auxiliary data
NT = size(elem,1);
if ~iscell(elem), elem = num2cell(elem,2); end
% diameter
diameter = cellfun(@(index) max(pdist(node(index,:))), elem(idElemMarked));
if max(diameter)<tol
    disp('The mesh is too dense'); return; 
end 
% totalEdge
shiftfun = @(verts) [verts(2:end),verts(1)];
T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
v0 = horzcat(elem{:})'; v1 = horzcat(T1{:})'; 
totalEdge = sort([v0,v1],2);
% edge, elem2edge
[edge, i1, totalJ] = unique(totalEdge,'rows');
elemLen = cellfun('length',elem); 
elem2edge = mat2cell(totalJ',1,elemLen)';
% edge2elem
Num = num2cell((1:NT)');    Len = num2cell(elemLen);
totalJelem = cellfun(@(n1,n2) n1*ones(n2,1), Num, Len, 'UniformOutput', false);
totalJelem = vertcat(totalJelem{:});
i2(totalJ) = 1:length(totalJ); i2 = i2(:); 
edge2elem = totalJelem([i1,i2]);
% neighbor
neighbor = cell(NT,1); 
% number
N = size(node,1); NE = size(edge,1);

ismElem = cell(NT,1); % store the elementwise mid-point locations
%% Determine the trivial and nontrivial marked elements
% nontrivial marked elements: marked elements with hanging nodes
nMarked = length(idElemMarked);
isT = false(nMarked,1);
for s = 1:nMarked
    % current element
    iel = idElemMarked(s);
    % local logical index of elements with hanging nodes
    p0 = node(elem{iel},:); p1 = circshift(p0,1,1); p2 = circshift(p0,-1,1);
    err = vecnorm(p0-0.5*(p1+p2),2,2);
    ism = (err<tol);  ismElem{iel} = ism;  % is midpoint    
    if sum(ism)<1, isT(s) = true; end  % trivial    
end
idElemMarkedT = idElemMarked(isT);
idElemMarkedNT = idElemMarked(~isT);

%% Find the additional elements to be refined
% initialized as marked elements
idElemMarkedNew = idElemMarked; % marked and all new elements
idElemNew = idElemMarked; % new elements generated in current step
isEdgeMarked = false(NE,1); 
while ~isempty(idElemNew)
    % adjacent polygons of new elements
    for iel = idElemNew(:)'
        if ~isempty(neighbor{iel}); continue; end
        index = elem2edge{iel};
        ia = edge2elem(index,1); ib = edge2elem(index,2);
        ia(ia==iel) = ib(ia==iel);
        neighbor{iel} = ia';
    end
    idElemNewNeighbor = unique(horzcat(neighbor{idElemNew}));
    idElemNewNeighbor = setdiff(idElemNewNeighbor,idElemMarkedNew); % delete the ones in the new marked set
    % edge set of new marked elements
    isEdgeMarked(horzcat(elem2edge{idElemNew})) = true;
    % find the adjacent elements to be refined
    nElemNewAdj = length(idElemNewNeighbor);
    isRefine = false(nElemNewAdj,1);
    for s = 1:nElemNewAdj
        % current element
        iel = idElemNewNeighbor(s);
        index = elem{iel};  indexEdge = elem2edge{iel}; Nv = length(index);
        % local logical index of elements with hanging nodes
        v1 = [Nv,1:Nv-1]; v0 = 1:Nv;  % left,current
        p0 = node(elem{iel},:); p1 = circshift(p0,1,1); p2 = circshift(p0,-1,1);
        err = vecnorm(p0-0.5*(p1+p2),2,2);
        ism = (err<tol);  ismElem{iel} = ism;  % is midpoint    
        if sum(ism)<1, continue; end  % start the next loop if no hanging nodes exist
        % index numbers of edges connecting hanging nodes in the adjacent elements to be refined
        idEdgeDg = indexEdge([v1(ism),v0(ism)]);
        % whether or not the above edges are in the edge set of new marked elements
        if sum(isEdgeMarked(idEdgeDg)), isRefine(s) = true; end
    end
    idElemNew = idElemNewNeighbor(isRefine);
    idElemMarkedNew = unique([idElemMarkedNew(:); idElemNew(:)]);
end
idElemRefineAddL = setdiff(idElemMarkedNew,idElemMarked);

%% Partition the nontrivial elements to be refined
% these elements are composed of: 
%  - nontrivial marked elements 
%  - additional elements to be refined
idElemRefineNT = [idElemMarkedNT; idElemRefineAddL];
nRefineNT = length(idElemRefineNT);
elemRefineNT = cell(nRefineNT,1);
elem2edgeRefineNT = cell(nRefineNT,1);
for s = 1:nRefineNT
    % current element
    iel = idElemRefineNT(s);
    index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
    % find midpoint
    v1 = [Nv,1:Nv-1]; v0 = 1:Nv; 
    ism = ismElem{iel};
    % modify the edge number
    ide = indexEdge+N;  % the connection number
    ide(v1(ism)) = index(ism); ide(ism) = index(ism);
    % elem (with or without hanging nodes)
    nsub = Nv-sum(ism);
    z1 = ide(v1(~ism));  z0 = index(~ism);
    z2 = ide(~ism);      zc = iel*ones(nsub,1)+N+NE;
    elemRefineNT{s} = [z1(:), z0(:), z2(:), zc(:)];
    % elem2edge
    ise = false(Nv,1); ise([v1(ism),v0(ism)]) = true;
    idg = zeros(Nv,1); idg(ise) = indexEdge(ise);
    e1 = idg(v1(~ism));     e0 = idg(~ism);
    elem2edgeRefineNT{s} = [e1(:), e0(:), zeros(nsub,2)]; % e2 = ec = 0
end
addElemRefineNT = num2cell(vertcat(elemRefineNT{:}), 2);
addElemRefineNT2edge = num2cell(vertcat(elem2edgeRefineNT{:}), 2); 

%% Determine the elements to be expanded
% these elements are composed of 
% - adjacent polygonals of elements to be refined
% - subcells of nontrivial elements
% adjacent polygons of elements to be refined
idElemRefine = unique([idElemRefineAddL(:); idElemMarked(:)]); % ascending order for update
for iel = idElemRefine(:)'
    if ~isempty(neighbor{iel}); continue; end
    index = elem2edge{iel};
    ia = edge2elem(index,1); ib = edge2elem(index,2);
    ia(ia==iel) = ib(ia==iel);
    neighbor{iel} = ia';
end
idElemRefineNeighbor = unique(horzcat(neighbor{idElemRefine}));
idElemRefineNeighbor = setdiff(idElemRefineNeighbor,idElemRefine);
% basic data structure of elements to be extended
elemExtend = [elem(idElemRefineNeighbor); addElemRefineNT]; 
elem2edgeExtend = [elem2edge(idElemRefineNeighbor); addElemRefineNT2edge];

%% Extend elements by adding hanging nodes 
% natural numbers of trivial edges w.r.t some element to be refined
isEdgeCut = false(NE,1);
for s = 1:length(idElemRefine)
    iel = idElemRefine(s);
    index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
    v1 = [Nv,1:Nv-1]; 
    ism = ismElem{iel};
    idx = true(Nv,1); idx(v1(ism)) = false; idx(ism) = false;
    isEdgeCut(indexEdge(idx)) = true;
end
% extend the elements
for s = 1:length(elemExtend)
    index = elemExtend{s}; indexEdge = elem2edgeExtend{s};
    id = find(indexEdge>0); 
    id = id(isEdgeCut(indexEdge(id)));
    idvec = zeros(1,2*length(index));
    idvec(1:2:end) = index;  idvec(2*id) = indexEdge(id)+N;
    elemExtend{s} = idvec(idvec>0);
end

%% Partition the trivial marked elements
nMarkedT = length(idElemMarkedT);
addElemMarkedT = cell(nMarkedT,1);
for s = 1:nMarkedT
    iel = idElemMarkedT(s); 
    index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
    ide = indexEdge+N;  % connection number    
    z1 = ide([Nv,1:Nv-1]);   z0 = index;
    z2 = ide;                zc = iel*ones(Nv,1)+N+NE;
    addElemMarkedT{s} = [z1(:), z0(:), z2(:), zc(:)];
end
% addElem
addElemMarkedT = num2cell(vertcat(addElemMarkedT{:}), 2);
addElem = [elemExtend; addElemMarkedT];

%% Update node and elem
% node
nodeEdgeCut = (node(edge(isEdgeCut,1),:) + node(edge(isEdgeCut,2),:))/2;
nodeCenter = zeros(length(idElemRefine),2);  
for s = 1:length(idElemRefine)
    iel = idElemRefine(s);  index = elem{iel};
    verts = node(index(~ismElem{iel}),:); verts1 = verts([2:end,1],:); 
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    ar = 0.5*abs(sum(area_components));
    nodeCenter(s,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar); 
end
node = [node; nodeEdgeCut; nodeCenter];
% elem
elem([idElemRefine(:); idElemRefineNeighbor(:)]) = []; % delete old
elem = [elem; addElem];

%% Reorder the vertices
[~,~,totalid] = unique(horzcat(elem{:})');
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';