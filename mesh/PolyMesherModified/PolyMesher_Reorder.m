function [node,elem] = PolyMesher_Reorder(NT,node,elem)

elem = elem(1:NT);
[id,~,totalid] = unique(horzcat(elem{:})');
node = node(id,:);
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';
for iel = 1:NT
    index = elem{iel};
    z1 = node(index(:,1),:); 
    z2 = node(index(:,2),:); 
    z3 = node(index(:,3),:);
    e2 = z3-z1; e3 = z1-z2;
    area = 0.5*(-e3(:,1).*e2(:,2)+e3(:,2).*e2(:,1));
    if area<0,  elem{iel} = index(end:-1:1);  end    
end