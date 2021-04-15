function [node,elem] = PolyMesher_rm_smalledge(node,elem)
% remove small edges
[node,elem] = PolyMesher_Reorder(node,elem);
T1 = cellfun(@(verts) [verts(2:end),verts(1)], elem, 'UniformOutput', false);
v0 = horzcat(elem{:})'; v1 = horzcat(T1{:})'; 
totalEdge = sort([v0,v1],2);
[i,j] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:);
he = sqrt(sum((z1-z2).^2,2));
ir = find(he < 0.1*max(he));
nv1 = edge(ir,1); nv2 = edge(ir,2); % starting and ending indices
[id,~,totalid] = unique(horzcat(elem{:})');
totalid = id(totalid);
for i = 1:length(nv2)
    v1 = nv1(i); v2 = nv2(i);
    totalid(totalid==v2) = v1;
    nv1(nv1==v2) = v1;  nv2(nv2==v2) = v1;
end
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';
for iel = 1:size(elem,1)
    index = elem{iel};
    [~,i1] = unique(index);
    elem{iel} = index(sort(i1));
end
[node,elem] = PolyMesher_Reorder(node,elem);