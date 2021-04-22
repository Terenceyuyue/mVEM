function [node,elem] = PolyMeshRefine2(node,elem,isMarkedElem)

tol = 1e-10; % tolerance for detecting colinearity
if ~iscell(elem)
    elem = num2cell(elem,2);
end

% number of elements
NT = numel(elem);
% number of vertices
N = size(node, 1);
% polygon shape
elemLen = cellfun('length', elem);

if ~islogical(isMarkedElem)
    markedElem = isMarkedElem;
    isMarkedElem = false(1,NT);
    isMarkedElem(markedElem) = true;
    markedElem = unique(markedElem);
else
    markedElem = find(isMarkedElem);
end
isMarkedElem = isMarkedElem(:).';
markedElem = markedElem(:).';

%% MAKE EDGE LIST
% make local edges(EDGE) from polygon (edges are repeated for interior edges)
edge = cellfun(@(c) [c; c([2:end,1])], elem, 'un', false);
edge = [edge{:}].';
% remove duplication
% store mapping between unique list and original list
[~, ind, bak] = unique(sort(edge,2), 'rows');
edge = edge(ind,:);
elem2edge = mat2cell(bak.', 1, elemLen).';

% number of edges
NE = size(edge,1);

%% MAPPING FROM EDGE TO ELEMENT
% find whether a local edge is negative edge or not
% by observing whether it is second edge or not
isSecondEdge = true(length(bak), 1);
isSecondEdge(ind) = false; % if it appears in the unique list, it is the first edge

% get element index of local edges
EDGE2elem = arrayfun(@(i,n) repmat(i, 1, n), 1:NT, elemLen.', 'un', 0);
EDGE2elem = [EDGE2elem{:}];

% locid2elem = arrayfun(@(n) 1:n, elemLen, 'un', 0);
% locid2elem = [locid2elem{:}];

% element index of global edges
edge2elem = zeros(NE, 2);
edge2elem(isSecondEdge*NE + bak) = EDGE2elem;

% if there is no negative element (edge2elem(i,2)==0)
% it is boundary
isBoundary = edge2elem(:,2) == 0;
% set negative element as itself for boundary to avoid 0 indexing
edge2elem(isBoundary,2) = edge2elem(isBoundary,1);

%% MARK EDGES FROM CURRENT MARKED ELEMENTS

isRefinedEdge = false(1, NE); % flag for refined edge
isMid4elem = cell(1,NT); % flag for vertex being midpoint
for i = markedElem
    if ~isMarkedElem(i) % if not marked,
        continue % move to next
    end
    
    % local information
    verts = elem{i};
    edgeidx = elem2edge{i};
    p0 = node(verts, :);
    
    % Get neighboring point information
    % left and right
    pL = circshift(p0,  1, 1);
    pR = circshift(p0, -1, 1);
    
    % check midpoint by |mid(pL,pR) - p0| / |edge(pLpR)|
    % divide error by left edge length for scale invariant error check
    isMidPoint = vecnorm((pL + pR)/2 - p0, 2, 2)./vecnorm((pL + pR)/2, 2, 2) < tol;
    isMid4elem{i} = isMidPoint;
    midpoints = find(isMidPoint);
    % find local edges refined
    % due to current element refinement
    isLocalRefinedEdge = true(length(verts),1);
    if ~isempty(midpoints)
        % if edge contains midpoint, do not refine
        isLocalRefinedEdge(midpoints) = false;
        isLocalRefinedEdge(midpoints-1) = false;
    end
    isRefinedEdge(edgeidx(isLocalRefinedEdge)) = true;
end

%% PROPAGATE REFINEMENT
isCandidate = false(1,NT); % candidate for refinement
while true
    % update candidates by not marked and neighbor of marked
    isCandidate(:) = false;
    isCandidate(edge2elem(isRefinedEdge, :)) = true;
    isCandidate = isCandidate & ~isMarkedElem;
    % isCandidate will be updated in the loop.
    % if isCandidate(i) is true after this loop,
    % it should be refined
    
    candidates = find(isCandidate);
    for i = candidates % for each candidate
        
        % get local info
        verts = elem{i};
        edgeidx = elem2edge{i};
        
        % get midpoint info
        if ~isempty(isMid4elem{i}) % if checked before
            isMidPoint = isMid4elem{i}; % get it
        else % check midpoint
            % local information
            p0 = node(verts, :);

            % Get neighboring point information
            % left and right
            pL = circshift(p0,  1, 1);
            pR = circshift(p0, -1, 1);

            % check midpoint by |mid(pL,pR) - p0| / |edge(pLpR)|
            % divide error by left edge length for scale invariant error check
            isMidPoint = vecnorm((pL + pR)/2 - p0, 2, 2)./vecnorm((pL + pR)/2, 2, 2) < tol;
            isMid4elem{i} = isMidPoint;
        end
        midpoints = find(isMidPoint);
        
        % if there is no midpoints in the current element
        if isempty(midpoints) % then it cannot have two-level hanging node
            isCandidate(i) = false; % should not be refined
            continue; % move on to next element
        end
        
        
        % find hanging edges
        isHangingEdge = false(1,length(verts));
        isHangingEdge(midpoints-1) = true; % mark left edge
        isHangingEdge(midpoints) = true; % mark right edge
        
        % if there is no refined edges among hanging edges
        if ~any(isRefinedEdge(edgeidx(isHangingEdge)))
            isCandidate(i) = false; % it should not be refined
            continue; % move on to next element
        end
        
        
        % If you are here, then it should be refined.
        % find local edges to be refined
        % due to current element refinement
        isLocalRefinedEdge = true(length(verts),1);
        % if edge contains midpoint, do not refine
        isLocalRefinedEdge(midpoints-1) = false; % left is false
        isLocalRefinedEdge(midpoints) = false; % right is false
        
        % mark edges should be refined
        isRefinedEdge(edgeidx(isLocalRefinedEdge)) = true;
    end
    
    if ~any(isCandidate) % if there is no new elements to be refined
        break;
    end
    % mark candidates for refinement
    isMarkedElem(isCandidate) = true;
end

%% REFINE POLYGON
markedElem = find(isMarkedElem);
elemIdx = cumsum([0, cellfun(@(v) nnz(~v), isMid4elem(isMarkedElem))]);
newElem = cell(elemIdx(end), 1);
centroid = zeros(length(markedElem), 2);

idxEdge0 = N;
idxCent0 = N + NE;
cyc = @(i, n) mod(i-1, n) + 1;
for id = 1 : length(markedElem)
    i = markedElem(id);

    % local information
    verts = elem{i};
    edgeidx = elem2edge{i};
    isMidPoint = isMid4elem{i};
    midpoints = find(isMidPoint);
    % next/previos of midpoint
    prevPoint = midpoints-1; % midpoint cannot be the first point. So, there is no need for cyclic indexing
    nextPoint = cyc(midpoints + 1, elemLen(i)); % cyclic indexing to avoid out of range
    
    
    % make subdivision
    curelem = [
        verts; % 1 current
        repmat(N + NT + NE + 1, 1, elemLen(i)); % 2 dummy for hanging node
        idxEdge0 + edgeidx; % 3 next
        repelem(idxCent0 + i, 1, elemLen(i)); % 4 centroid
        idxEdge0 + edgeidx([end,1:end-1]) % 5 previous
        repmat(N + NT + NE + 1, 1, elemLen(i)); % 6 dummy for hanging nodes
        ].';

    % update elements additional vertex by midpoint
    % note that mid cannot be the first/last vertex
    curelem(prevPoint,3) = verts(isMidPoint); % update previous element
    curelem(nextPoint,5) = verts(isMidPoint); % update next element

    % check left and right edge of midpoint for further refinement
    isNeighborMarked = isRefinedEdge(elem2edge{i});
    isMarkedL = isNeighborMarked(midpoints-1);
    isMarkedR = isNeighborMarked(midpoints);
    % update hanging nodes
    curelem(prevPoint(isMarkedL), 2) = idxEdge0 + edgeidx(midpoints(isMarkedL)-1);
    curelem(nextPoint(isMarkedR), 6) = idxEdge0 + edgeidx(midpoints(isMarkedR));
    
    % make it cell
    curelem = num2cell(curelem, 2);
    % remove 0's and gather
    newElem(elemIdx(id) +  1 : elemIdx(id+1)) = curelem(~isMidPoint);
    
    verts = node(verts(~isMidPoint),:); verts1 = verts([2:end,1],:); 
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    ar = 0.5*abs(sum(area_components));
    centroid(id,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
end

%% UPDATE NEIGHBORING ELEMENTS
% if not marked & edge is refined -> neighbor element
isNeighbor = false(1, NT);
isNeighbor(edge2elem(isRefinedEdge, :)) = true;
neighbors = find(isNeighbor & ~isMarkedElem);
flat = @(v) v(:);
% update neighbors
for i = neighbors
    elem{i} = flat([elem{i};
                    idxEdge0 + elem2edge{i}]).';
end

%% NODE UPDATE

node = [node; (node(edge(isRefinedEdge, 1),:) + node(edge(isRefinedEdge,2),:))/2; centroid];

%% REMOVE DUMMY NODES
isAlive = true(N + NT + NE + 1, 1);
isAlive(N+1:N+NE) = isRefinedEdge;
isAlive(N+NE+1:N+NE+NT) = isMarkedElem;
isAlive(end) = false;

newNumbering = zeros(1,length(isAlive));
newNumbering(isAlive) = 1:nnz(isAlive);

elem(isMarkedElem) = [];
elem(end + 1 : end + length(newElem)) = newElem;
elem = cellfun(@(v) newNumbering(v(isAlive(v))), elem, 'un', 0);

end