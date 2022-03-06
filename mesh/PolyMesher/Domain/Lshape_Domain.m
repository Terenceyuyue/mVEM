%------------- PolyMesher -------------------%
% Ref: L ÐÍÇøÓò
%--------------------------------------------%
function x = Lshape_Domain(Demand,Arg)
  BdBox = [0 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  x1 = BdBox(1); x2 = BdBox(2); 
  y1 = BdBox(3); y2 = BdBox(4);
  xm = (x1+x2)/2; ym = (y1+y2)/2;
  d1 = dRectangle(P,x1,x2,y1,y2);
  d2 = dRectangle(P,xm,x2,y1,ym);
  Dist = dDiff(d1,d2);
% %------------- SPECIFY BOUNDARY CONDITIONS
% function [x] = BndryCnds(Node,Element,BdBox)
% x = cell(2,1);
%------------- SPECIFY FIXED POINTS
function PFix = FixedPoints(BdBox)
  x1 = BdBox(1); x2 = BdBox(2); 
  y1 = BdBox(3); y2 = BdBox(4);
  xm = (x1+x2)/2; ym = (y1+y2)/2;
  NT = 200; h = 1/NT;
  PFix = [ xm-h, ym-h;
           xm+h, ym+h;
           xm-h, ym+h];
%------------------------------------%