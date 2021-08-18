function [x] = Lshape_Domain(Demand,P)
BdBox = [0 2 0 1 0 2];
switch(Demand)
    case('Dist');  x = DistFnc(P,BdBox);
    case('BdBox'); x = BdBox;
end
%--------- the signed distance function ----------
function Dist = DistFnc(P,BdBox)
Dist = dLshape(P,BdBox);