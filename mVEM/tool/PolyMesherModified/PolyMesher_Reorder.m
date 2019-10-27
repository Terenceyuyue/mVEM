function [node,elem] = PolyMesher_Reorder(NT,node,elem)
% Extract node list
elem = elem(1:NT)';
[id,~,totalid] = unique(horzcat(elem{:})');
node = node(id,:);
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';