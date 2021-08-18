function d = dLshape(P,BdBox)

% BdBox = [0 2 0 1 0 2];

BdBox2 = BdBox;
BdBox2(1) = (BdBox(1)+BdBox(2))/2;
BdBox2(5) = (BdBox(5)+BdBox(6))/2;

d1 = dCube(P,BdBox);
d2 = dCube(P,BdBox2);
d = dDiff(d1,d2);