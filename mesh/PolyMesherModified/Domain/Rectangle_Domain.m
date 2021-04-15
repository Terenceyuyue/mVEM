%------------- PolyMesher -------------------%
% Ref: 矩形区域
% 修改 BdBox 即可
%--------------------------------------------%
function x = Rectangle_Domain(Demand,P)
  BdBox = [0 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(P,BdBox);
    case('BdBox'); x = BdBox;
    case('PFix'); x = [];
  end
%------------- the signed distance function --------------
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
