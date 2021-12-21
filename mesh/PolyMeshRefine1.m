function [node,elem] = PolyMeshRefine1(node,elem,elemMarked)
% a new implementation of PolyMeshRefine
%
% Copyright (C) Dohyun, modified by Terence Yu.

tol = 1e-10; % tolerance for detecting colinearity
if ~iscell(elem), elem = num2cell(elem,2); end

N = size(node,1); NT = numel(elem);
elemLen = cellfun('length', elem);

idElemMarked = unique(elemMarked(:))';
isMarkedElem = false(1,NT);
isMarkedElem(idElemMarked) = true;


%% Get auxiliary data
% diameter
diameter = cellfun(@(index) max(pdist(node(index,:))), elem(idElemMarked));
if max(diameter)<tol
    disp('The mesh is too dense');  return; 
end 
% edge
edge = cellfun(@(index) [index;index([2:end,1])], elem, 'un', false);
edge = [edge{:}]';
[~, i1, totalJ] = unique(sort(edge,2), 'rows');
edge = edge(i1,:);
% elem2edge
elem2edge = mat2cell(totalJ', 1, elemLen)';
NE = size(edge,1);
% edge2elem
% find whether a local edge is negative edge or not
% by observing whether it is second edge or not
isSecondEdge = true(length(totalJ), 1);
isSecondEdge(i1) = false; % if it appears in the unique list, it is the first edge
% get element index of local edges
EDGE2elem = arrayfun(@(i,n) repmat(i,1,n), 1:NT, elemLen.', 'un', 0);
EDGE2elem = [EDGE2elem{:}];
% element index of global edges
edge2elem = zeros(NE, 2);
edge2elem(isSecondEdge*NE + totalJ) = EDGE2elem;
% if there is no negative element (edge2elem(i,2)==0)
% it is boundary
isBoundary = edge2elem(:,2) == 0;
% set negative element as itself for boundary to avoid 0 indexing
edge2elem(isBoundary,2) = edge2elem(isBoundary,1);

%% Mark edges from current marked elements
isRefinedEdge = false(1,NE); % flag for refined edge
isMid4elem = cell(1,NT); % flag for vertex being midpoint
for iel = idElemMarked
    % local information
    index = elem{iel}; indexEdge = elem2edge{iel}; Nv = length(index);
    p0 = node(index, :); pL = circshift(p0,1,1); pR = circshift(p0,-1,1);
    % check midpoint
    isMidPoint = (vecnorm(0.5*(pL+pR)-p0,2,2)<tol); 
    isMid4elem{iel} = isMidPoint;    
    % find local edges refined due to current element refinement 
    % if edge contains midpoint, do not refine
    isLocalRefinedEdge = true(Nv,1);  
    v1 = [Nv,1:Nv-1]; 
    isLocalRefinedEdge(v1(isMidPoint)) = false; 
    isLocalRefinedEdge(isMidPoint) = false;  
    isRefinedEdge(indexEdge(isLocalRefinedEdge)) = true;
end

%% Propagate refinement
isCandidate = true(1,NT);
while any(isCandidate)
    % update candidates by not marked and neighbor of marked
    isCandidate = false(1,NT); % candidate for refinement  
    isCandidate(edge2elem(isRefinedEdge,:)) = true; % include neighbor of marked
    isCandidate = (isCandidate & ~isMarkedElem); % delete marked

    for iel = find(isCandidate) % for each candidate
        % get local info
        index = elem{iel}; indexEdge = elem2edge{iel};  Nv = length(index);

        % get midpoint info
        isMidPoint = isMid4elem{iel}; % default for checked before
        if isempty(isMidPoint) % check midpoint
            p0 = node(index,:); 
            pL = circshift(p0,1,1); 
            pR = circshift(p0,-1,1);
            isMidPoint = (vecnorm(0.5*(pL+pR)-p0,2,2) < tol);
            isMid4elem{iel} = isMidPoint;
        end
        
        % do not refine if no hanging node and move on to next element
        if ~any(isMidPoint), isCandidate(iel) = false; continue; end         
        % else: find hanging edges
        isHangingEdge = false(1,Nv);
        v1 = [Nv,1:Nv-1]; 
        isHangingEdge(v1(isMidPoint)) = true; % left edge
        isHangingEdge(isMidPoint) = true; % right edge

        % do not refine if there is no refined edges among hanging edges
        if ~any(isRefinedEdge(indexEdge(isHangingEdge)))
            isCandidate(iel) = false; continue; % move on to next element
        end
        % else: it should be refined
        % find local edges to be refined
        isLocalRefinedEdge = true(Nv,1);
        % if edge contains midpoint, do not refine
        isLocalRefinedEdge(v1(isMidPoint)) = false; 
        isLocalRefinedEdge(isMidPoint) = false;
        
        % mark edges that should be refined
        isRefinedEdge(indexEdge(isLocalRefinedEdge)) = true;
    end
    
    % mark candidates for refinement
    isMarkedElem(isCandidate) = true;
end

%% Refine polygons
idElemMarked = find(isMarkedElem);
nMarked = length(idElemMarked);
elemIdx = cumsum([0, cellfun(@(v) nnz(~v), isMid4elem(isMarkedElem))]);
newElem = cell(elemIdx(end), 1);
centroid = zeros(nMarked, 2);

cyc = @(i, n) mod(i-1, n) + 1;
for id = 1:nMarked
    iel = idElemMarked(id);

    % local information
    index = elem{iel}; indexEdge = elem2edge{iel};
    isMidPoint = isMid4elem{iel};
    midpoints = find(isMidPoint);
    % next/previos of midpoint
    prevPoint = midpoints-1; % midpoint cannot be the first point due to arrangement
    nextPoint = cyc(midpoints+1, elemLen(iel)); % cyclic


    % make subdivision
    curelem = [
        index; % 1-current
        repmat(N+NT+NE+1, 1, elemLen(iel)); % 2-dummy for hanging node
        N + indexEdge; % 3-next
        repelem(N+NE + iel, 1, elemLen(iel)); % 4-centroid
        N + indexEdge([end,1:end-1]) % 5-previous
        repmat(N+NT+NE+1, 1, elemLen(iel)); % 6-dummy for hanging nodes
        ];
    curelem = curelem';

    % update elements additional vertex by midpoint
    % note that mid cannot be the first/last vertex
    curelem(prevPoint,3) = index(isMidPoint); % update previous element
    curelem(nextPoint,5) = index(isMidPoint); % update next element

    % check left and right edge of midpoint for further refinement
    isNeighborMarked = isRefinedEdge(elem2edge{iel});
    isMarkedL = isNeighborMarked(midpoints-1);
    isMarkedR = isNeighborMarked(midpoints);
    % update hanging nodes
    curelem(prevPoint(isMarkedL), 2) = N + indexEdge(midpoints(isMarkedL)-1);
    curelem(nextPoint(isMarkedR), 6) = N + indexEdge(midpoints(isMarkedR));

    % make it cell
    curelem = num2cell(curelem, 2);
    % remove 0's and gather
    newElem(elemIdx(id)+1 : elemIdx(id+1)) = curelem(~isMidPoint);

    verts = node(index(~isMidPoint),:); verts1 = verts([2:end,1],:); 
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    ar = 0.5*abs(sum(area_components));
    centroid(id,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
end

%% Update neighboring elements 
% if not marked & edge is refined -> neighbor element
isNeighbor = false(1, NT);
isNeighbor(edge2elem(isRefinedEdge, :)) = true;
neighbors = find(isNeighbor & ~isMarkedElem);
% update neighbors
for iel = neighbors
    elem{iel} = reshape([elem{iel}; N + elem2edge{iel}], 1, []);
end

%% Update node
node = [node; ...
       (node(edge(isRefinedEdge, 1),:) + node(edge(isRefinedEdge,2),:))/2; ...
       centroid];

%% Remove dummy nodes 
isAlive = true(N+NT+NE+1, 1);
isAlive(N+1:N+NE) = isRefinedEdge;
isAlive(N+NE+1:N+NE+NT) = isMarkedElem;
isAlive(end) = false;

newNumbering = zeros(1,length(isAlive));
newNumbering(isAlive) = 1:nnz(isAlive);

elem(isMarkedElem) = [];
elem(end + 1 : end + length(newElem)) = newElem;
elem = cellfun(@(v) newNumbering(v(isAlive(v))), elem, 'un', 0);

end