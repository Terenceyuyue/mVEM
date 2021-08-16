%------------- PolyMesher -------------------%
% Ref: 矩形区域
% 修改 BdBox 即可
%--------------------------------------------%
function x = Cube_Domain(Demand,P)
  BdBox = [0 1 0 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(P,BdBox);
    case('BdBox'); x = BdBox;
  end
%------------- the signed distance function --------------
function Dist = DistFnc(P,BdBox)
  Dist = dCube(P,BdBox);
