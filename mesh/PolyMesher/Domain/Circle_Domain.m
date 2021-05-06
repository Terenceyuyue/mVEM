%------------- PolyMesher -------------------%
% Ref: 圆形区域
% 修改 BdBox,dCircle 即可
%--------------------------------------------%
function [x] = Circle_Domain(Demand,P)
BdBox = [-1 1 -1 1];
switch(Demand)
    case('Dist');  x = DistFnc(P,BdBox);
    case('BdBox'); x = BdBox;
    case('PFix'); x = [];
end
%--------- the signed distance function ----------
function Dist = DistFnc(P,BdBox)
Dist = dCircle(P,0,0,1);