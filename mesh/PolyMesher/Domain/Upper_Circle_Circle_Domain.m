%------------------------------ PolyMesher -------------------------------%
% Ref: 上半圆环区域
% dLine: 左右点连线 (线的左侧区域是区域内部，左侧按逆时针理解)
%-------------------------------------------------------------------------%
function [x] = Upper_Circle_Circle_Domain(Demand,Arg)
  BdBox = [-1 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dCircle(P,0,0,1);
  d2 = dCircle(P,0,0,0.5);
  d3 = dLine(P,0,0,1,0);
  Dist = dIntersect(d3,dDiff(d1,d2));

%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  x = cell(2,1); %No boundary conditions specified for this problem
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%