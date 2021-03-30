function d = dCube(P,BdBox)
x1 = BdBox(1); x2 = BdBox(2); 
y1 = BdBox(3); y2 = BdBox(4);
z1 = BdBox(5); z2 = BdBox(6);
d = [x1-P(:,1), P(:,1)-x2, y1-P(:,2), P(:,2)-y2, z1-P(:,3), P(:,3)-z2];
d = [d,max(d,[],2)];
