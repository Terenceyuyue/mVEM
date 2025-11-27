function aux = auxstructure(node,elem)
%auxstructure gets auxiliary mesh data structure
%
% Copyright (C) Terence Yu.

NT = size(elem,1); N = size(node,1);
if ~iscell(elem) % transform to cell
    elem = mat2cell(elem,ones(NT,1),length(elem(1,:)));
end

% totalEdge
shiftfun = @(verts) [verts(2:end),verts(1)];  % or shiftfun = @(verts) circshift(verts,-1);
T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
v0 = horzcat(elem{:})'; % the starting points of edges
v1 = horzcat(T1{:})'; % the ending points of edges
allEdge = [v0,v1];
totalEdge = sort(allEdge,2);

% -------- elem2edge: elementwise edges -------
[~, i1, totalJ] = unique(totalEdge,'rows'); % first occurence
elemLen = cellfun('length',elem); % length of each elem
elem2edge = mat2cell(totalJ',1,elemLen)';

% -------- edge, bdEdge, bdNode --------
[i,j,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
%bdEdge = edge(s==1,:); % not counterclockwise
bdEdge = allEdge(i1(s==1),:); % counterclockwise
bdEdgeIdx = find(s==1);      % index of all boundary edges
bdNodeIdx = unique(bdEdge);
bdNode = node(bdNodeIdx,:);

% ------- edge2elem --------
Num = num2cell((1:NT)');    Len = num2cell(elemLen);
totalJelem = cellfun(@(n1,n2) n1*ones(n2,1), Num, Len, 'UniformOutput', false);
totalJelem = vertcat(totalJelem{:});
i2(totalJ) = 1:length(totalJ); i2 = i2(:); % second occurence by overlapping
edge2elem = totalJelem([i1,i2]);

% --------- neighbor ---------
neighbor = cell(NT,1);
for iel = 1:NT
    index = elem2edge{iel};
    ia = edge2elem(index,1); ib = edge2elem(index,2);
    ia(ia==iel) = ib(ia==iel);
    neighbor{iel} = ia';
end

% -------- edge2elemInd: add local edge indices in edge2elem -------
totalJLocalEdge = cellfun(@(nv) (1:nv)', Len, 'UniformOutput', false);
totalJLocalEdge = vertcat(totalJLocalEdge{:});
edge2elemLocalIndex = totalJLocalEdge([i1,i2]);

% --------- node2elem ---------
ii = totalJelem; jj = v0; ss = ones(length(ii),1);
t2v = sparse(ii,jj,ss,NT,N);
node2elem = cellfun(@(x) find(x), num2cell(t2v,1) ,'UniformOutput', false);
node2elem = node2elem';

% ------- elem2sgn: sign of edges on each element ----------------
% the positive sign satisfies: e(k,1) < e(k,2)
% on the boundary of domain: the sign is set as zero
elem2sgn = cell(NT,1);
NE = size(edge,1);
E = false(NE,1); E(bdEdgeIdx) = 1;
for iel = 1:NT
    index = elem{iel}; Nv = length(index);
    sgnedge = sign(index([2:Nv,1]) - index);
    id = elem2edge{iel}; sgnbd = E(id); sgnedge(sgnbd) = 0;  % on the domain boundary
    elem2sgn{iel} = sgnedge;
end


aux.node = node; aux.elem = elem;
aux.elem2edge = elem2edge;
aux.edge = edge;
aux.bdEdge = bdEdge;  aux.bdEdgeIdx = bdEdgeIdx;
aux.bdNode = bdNode;  aux.bdNodeIdx = bdNodeIdx;
aux.edge2elem = edge2elem;
aux.edge2elemLocalIndex = edge2elemLocalIndex;
aux.neighbor = neighbor; aux.node2elem = node2elem;
aux.elem2sgn = elem2sgn;
aux.elemLen = cellfun('length',elem);