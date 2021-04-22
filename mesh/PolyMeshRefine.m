function [node,elem] = PolyMeshRefine(node,elem,elemMarked)
%PolyMeshRefine refines a 2-D polygonal mesh satisfying one-hanging node rule
%
% We divide elements by connecting the midpoint of each edge to its
% barycenter.
% We remove small edges by further partitioning some adjacent elements.
%
% Copyright (C) Terence Yu.

idElemMarked = unique(elemMarked(:)); % in ascending order
eps = 1e-10; % accuracy for finding midpoint

%% Get auxiliary data
NT = size(elem,1);
if ~iscell(elem), elem = mat2cell(elem,ones(NT,1),length(elem(1,:))); end
% centroid
centroid = zeros(NT,2); diameter = zeros(NT,1); 
s = 1;
for iel = 1:NT
    index = elem{iel};
    verts = node(index,:); verts1 = verts([2:end,1],:); 
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    ar = 0.5*abs(sum(area_components));
    centroid(s,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
    diameter(s) = max(pdist(verts));    
    s = s+1;
end
if max(diameter)<4*eps, error('The mesh is too dense'); end
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
[~, i2] = unique(totalJ(end:-1:1),'rows');
i2 = length(totalEdge)+1-i2;
edge2elem = totalJelem([i1,i2]);
% neighbor
neighbor = cell(NT,1);
for iel = 1:NT
    index = elem2edge{iel};  
    ia = edge2elem(index,1); ib = edge2elem(index,2);
    ia(ia==iel) = ib(ia==iel);
    neighbor{iel} = ia';
end
% number
N = size(node,1); NE = size(edge,1);

%% Determine the trivial and nontrivial marked elements
% nontrivial marked elements: marked elements with hanging nodes
nMarked = length(idElemMarked);
isT = false(nMarked,1);
for s = 1:nMarked
    % current element
    iel = idElemMarked(s);
    index = elem{iel};  Nv = length(index);
    % local logical index of elements with hanging nodes
    v1 = [Nv,1:Nv-1]; v0 = 1:Nv; v2 = [2:Nv,1]; % left,current,right
    p1 = node(index(v1),:); p0 = node(index(v0),:); p2 = node(index(v2),:);
    err = sqrt(sum((p0-0.5*(p1+p2)).^2,2));
    ism = (err<eps);  % is midpoint
    if sum(ism)<1, isT(s) = true; end  % trivial    
end
idElemMarkedT = idElemMarked(isT);
idElemMarkedNT = idElemMarked(~isT);

%% Find the additional elements to be refined
% initialized as marked elements
idElemMarkedNew = idElemMarked; % marked and all new elements
idElemNew = idElemMarked; % new elements generated in current step
while ~isempty(idElemNew)
    % adjacent polygons of new elements
    idElemNewAdj = unique(horzcat(neighbor{idElemNew}));
    idElemNewAdj = setdiff(idElemNewAdj,idElemMarkedNew); % delete the ones in the new marked set
    % edge set of new marked elements
    idEdgeMarkedNew = unique(horzcat(elem2edge{idElemMarkedNew}));
    % find the adjacent elements to be refined
    nElemNewAdj = length(idElemNewAdj);
    isRefine = false(nElemNewAdj,1);
    for s = 1:nElemNewAdj
        % current element
        iel = idElemNewAdj(s);
        index = elem{iel};  indexEdge = elem2edge{iel}; Nv = length(index);
        % local logical index of elements with hanging nodes
        v1 = [Nv,1:Nv-1]; v0 = 1:Nv; v2 = [2:Nv,1]; % left,current,right
        p1 = node(index(v1),:); p0 = node(index(v0),:); p2 = node(index(v2),:);
        err = sqrt(sum((p0-0.5*(p1+p2)).^2,2));
        ism = (err<eps);  % is midpoint
        if sum(ism)<1, continue; end  % start the next loop if no hanging nodes exist
        % index numbers of edges connecting hanging nodes in the adjacent elements to be refined
        idEdgeDg = unique(indexEdge([v1(ism),v0(ism)]));
        % whether or not the above edges are in the edge set of new marked elements
        if intersect(idEdgeDg, idEdgeMarkedNew), isRefine(s) = true; end
    end
    idElemNew = idElemNewAdj(isRefine);
    idElemMarkedNew = unique([idElemMarkedNew(:); idElemNew(:)]);
end
idElemRefineAddL = setdiff(idElemMarkedNew,idElemMarked);

%% Partition the additional elements
% additional elements are composed of: 
%  - nontrivial marked elements 
%  - additional elements to be refined
idElemAddL = [idElemMarkedNT; idElemRefineAddL];
nAddL = length(idElemAddL);
elemAddL = cell(nAddL,1);
elem2edgeAddL = cell(nAddL,1);
for s = 1:nAddL
    % current element
    iel = idElemAddL(s);
    index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
    % find midpoint
    v1 = [Nv,1:Nv-1]; v0 = 1:Nv; v2 = [2:Nv,1];
    p1 = node(index(v1),:); p0 = node(index(v0),:); p2 = node(index(v2),:);
    err = sqrt(sum((p0-0.5*(p1+p2)).^2,2));
    ism = (err<eps);
    % modify the edge number
    ide = indexEdge+N;  % the connection number
    ide(v1(ism)) = index(ism); ide(ism) = index(ism);
    % elem (with or without hanging nodes)
    nsub = Nv-sum(ism);
    z1 = ide(v1(~ism));  z0 = index(~ism);
    z2 = ide(~ism);      zc = iel*ones(nsub,1)+N+NE;
    elemAddL{s} = [z1(:), z0(:), z2(:), zc(:)];
    % elem2edge
    ise = false(Nv,1); ise([v1(ism),v0(ism)]) = true;
    idg = zeros(Nv,1); idg(ise) = indexEdge(ise);
    e1 = idg(v1(~ism));     e0 = idg(~ism);
    e2 = zeros(nsub,1);     ec = zeros(nsub,1);  % numbered as 0
    elem2edgeAddL{s} = [e1(:), e0(:), e2(:), ec(:)];
end
addElemAddL = vertcat(elemAddL{:});
addElemAddL2edge = vertcat(elem2edgeAddL{:}); % transform to cell arrays
if ~isempty(addElemAddL) % may be empty
    addElemAddL = mat2cell(addElemAddL,ones(size(addElemAddL,1),1),4);
    addElemAddL2edge = mat2cell(addElemAddL2edge,ones(size(addElemAddL2edge,1),1),4);
end

%% Determine the elements to be expanded
% these elements are composed of 
% - adjacent polygonals of elements to be refined
% - additional elements (see the def in Line 107)
% adjacent polygons of elements to be refined
idElemRefine = [idElemRefineAddL(:); idElemMarked(:)];
idElemRefineAdj = unique(horzcat(neighbor{idElemRefine}));
idElemRefineAdj = setdiff(idElemRefineAdj,idElemRefine);
% basic data structure of elements to be extended
elemExtend = [elem(idElemRefineAdj); addElemAddL]; 
elem2edgeExtend = [elem2edge(idElemRefineAdj); addElemAddL2edge];

%% Extend elements by adding hanging nodes
% natural numbers of trivial edges w.r.t some element to be refined
isEdgeCut = false(NE,1);
for s = 1:length(idElemRefine)
    iel = idElemRefine(s);
    index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
    v1 = [Nv,1:Nv-1]; v0 = 1:Nv; v2 = [2:Nv,1];
    p1 = node(index(v1),:); p0 = node(index(v0),:); p2 = node(index(v2),:);
    err = sqrt(sum((p0-0.5*(p1+p2)).^2,2));
    ism = (err<eps);
    idx = true(Nv,1);
    idx(v1(ism)) = false; idx(ism) = false;
    isEdgeCut(indexEdge(idx)) = true;
end
idEdgeCut = find(isEdgeCut);
% extend the elements
for s = 1:length(elemExtend)
    index = elemExtend{s}; indexEdge = elem2edgeExtend{s};
    Nv = length(index);
    [idm,id] = intersect(indexEdge,idEdgeCut);
    idvec = zeros(1,2*Nv);
    idvec(1:2:end) = index;  idvec(2*id) = idm+N;
    elemExtend{s} = idvec(idvec>0);
end
% replace the old elements
idElemOld = [idElemRefineAdj(:); idElemRefineAddL(:); idElemMarkedNT(:)];
nOld = length(idElemOld);  
elem(idElemOld) = elemExtend(1:nOld);
addElemAddL = elemExtend(nOld+1:end);

%% Partition the trivial marked elements
nMarkedT = length(idElemMarkedT);
addElemMarkedT = cell(nMarkedT,1);
for s = 1:nMarkedT
    iel = idElemMarkedT(s);
    index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
    v1 = [Nv,1:Nv-1];
    ide = indexEdge+N;  % connection number    
    z1 = ide(v1);   z0 = index;
    z2 = ide;       zc = iel*ones(Nv,1)+N+NE;
    addElemMarkedT{s} = [z1(:), z0(:), z2(:), zc(:)];
end
% replace the old elements 
addElemMarkedT = vertcat(addElemMarkedT{:});
addElemMarkedT = mat2cell(addElemMarkedT, ones(size(addElemMarkedT,1),1), 4);
elem(idElemMarkedT) = addElemMarkedT(1:nMarkedT);
addElemMarkedT = addElemMarkedT(nMarkedT+1:end);

%% Update node and elem
idElemRefine = unique(idElemRefine); % in ascending order
z1 = node(edge(idEdgeCut,1),:); z2 = node(edge(idEdgeCut,2),:);
nodeEdgeCut = (z1+z2)/2;
nodeCenter = centroid(idElemRefine,:);
node = [node; nodeEdgeCut; nodeCenter];
elem = [elem; addElemAddL; addElemMarkedT];

%% Reorder the vertices
[~,~,totalid] = unique(horzcat(elem{:})');
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';