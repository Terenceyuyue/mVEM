function Pc = circumcenter3(varargin)
% Compute the circumcenter of a 3-D triangle
%
% Reference: https://www.zhihu.com/question/40422123/answer/95637252
%
% Copyright (C)  Terence Yu.

if nargin==1
    P = varargin{1};
    P1 = P(1,:);  P2 = P(2,:);  P3 = P(3,:);
end
if nargin==3
    P1 = varargin{1};  P2 = varargin{2};  P3 = varargin{3};
end

x1 = P1(1);  y1 = P1(2);  z1 = P1(3);
x2 = P2(1);  y2 = P2(2);  z2 = P2(3);
x3 = P3(1);  y3 = P3(2);  z3 = P3(3);

a1 = (y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2);
b1 = -(x1*z2 - x2*z1 - x1*z3 + x3*z1 + x2*z3 - x3*z2);
c1 = (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
d1 = -(x1*y2*z3 - x1*y3*z2 - x2*y1*z3 + x2*y3*z1 + x3*y1*z2 - x3*y2*z1);

a2 = 2 * (x2 - x1);
b2 = 2 * (y2 - y1);
c2 = 2 * (z2 - z1);
d2 = x1 * x1 + y1 * y1 + z1 * z1 - x2 * x2 - y2 * y2 - z2 * z2;

a3 = 2 * (x3 - x1);
b3 = 2 * (y3 - y1);
c3 = 2 * (z3 - z1);
d3 = x1 * x1 + y1 * y1 + z1 * z1 - x3 * x3 - y3 * y3 - z3 * z3;

x = -(b1*c2*d3 - b1*c3*d2 - b2*c1*d3 + b2*c3*d1 + b3*c1*d2 - b3*c2*d1) / ...
    (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1);
y = (a1*c2*d3 - a1*c3*d2 - a2*c1*d3 + a2*c3*d1 + a3*c1*d2 - a3*c2*d1) / ...
     (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1);
z = -(a1*b2*d3 - a1*b3*d2 - a2*b1*d3 + a2*b3*d1 + a3*b1*d2 - a3*b2*d1) / ...
    (a1*b2*c3 - a1*b3*c2 - a2*b1*c3 + a2*b3*c1 + a3*b1*c2 - a3*b2*c1);

Pc = [x,y,z];