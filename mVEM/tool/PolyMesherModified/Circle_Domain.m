%------------- PolyMesher -------------------%
% Ref: 圆形区域
% 修改 BdBox,dCircle 即可
%--------------------------------------------%
function [x] = Circle_Domain(Demand,Arg)
BdBox = [-1 1 -1 1];
switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
end
%------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
d = dCircle(P,0,0,1);
Dist = d;
%------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
x = cell(2,1);
%------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
PFix = [];
%------------------------------------%