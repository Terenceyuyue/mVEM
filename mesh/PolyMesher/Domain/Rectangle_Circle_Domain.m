%------------- PolyMesher -------------------%
% Ref: 矩形与内部圆组成的“环形”区域
% 修改 BdBox,dCircle 即可
%--------------------------------------------%
function [x] = Rectangle_Circle_Domain(Demand,Arg)
  BdBox = [-2 2 -2 2];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  d2 = dCircle(P,0,0,1);
  Dist = dDiff(d1,d2);
%------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
x = cell(2,1);
%------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%------------------------------------%