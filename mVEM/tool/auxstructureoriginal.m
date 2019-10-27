function aux = auxstructure(node,elem)

if iscell(elem)
    % totalEdge
    shiftfun = @(verts) [verts(2:end),verts(1)];  % or shiftfun = @(verts) circshift(verts,-1);
    T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
    v0 = horzcat(elem{:})'; % the starting points of edges
    v1 = horzcat(T1{:})'; % the ending points of edges
    totalEdge = sort([v0,v1],2);
    
    % elem2edge
    [~, ~, totalJ] = unique(totalEdge,'rows');
    elemLen = cellfun('length',elem); % length of each elem
    elem2edge = mat2cell(totalJ,elemLen,1);
    elem2edge = cellfun(@transpose, elem2edge, 'UniformOutput', false);
    
else % Triangulation
    totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
    [~, ~, totalJ] = unique(totalEdge,'rows');
    NT = size(elem,1);
    elem2edge = reshape(totalJ,NT,3);
end

% -------- edge, bdEdge --------
[i,j,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
bdEdge = edge(s==1,:);


% ------- elemp2 -------
NT = size(elem,1); NV = size(node,1); NE = size(edge,1);
if iscell(elem)
    centroid = mat2cell((1:NT)',ones(NT,1));
    montage = @(x,y,z) [x,NV+y,NV+NE+z];
    elemp2 = cellfun(montage, elem, elem2edge, centroid, 'UniformOutput', false);
else
    elemp2 = [elem, NV+elem2edge, NV+NE+(1:NT)'];
end

% ------ nodeEdge -------
nodeEdge = 0.5*(node(edge(:,1),:) + node(edge(:,2),:));

% ------ nodeCenter, area, diameter -------
nodeCenter = zeros(NT,2); area = zeros(NT,1); diameter = zeros(NT,1);
s = 1;
for iel = 1:NT
    if iscell(elem)
        index = elem{iel};
    else
        index = elem(iel,:);
    end
    verts = node(index, :); verts1 = verts([2:end,1],:);
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    ar = 0.5*abs(sum(area_components));
    area(iel) = ar;
    nodeCenter(s,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
    diameter(s) = max(pdist(verts));
    s = s+1;
end

if ~iscell(elem)
    elem = mat2cell(elem,ones(NT,1),3);
    elem2edge = mat2cell(elem2edge,ones(NT,1),3);
    elemp2 = mat2cell(elemp2,ones(size(elemp2,1),1),size(elemp2,2));
end

aux.node = node; aux.elem = elem;
aux.elem2edge = elem2edge;
aux.edge = edge;
aux.bdEdge = bdEdge;
aux.elemp2 = elemp2;
aux.nodeEdge = nodeEdge;
aux.nodeCenter = nodeCenter;
aux.area = area;
aux.diameter = diameter;


