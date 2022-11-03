function area = polyarea3(P)
% Compute the area of a polygon in 3-D
% See simplexvolume.m
%
% Copyright (C)  Terence Yu.

nv = size(P,1);
tri = [ones(nv-2,1), (2:nv-1)', (3:nv)'];
d12 = P(tri(:,2),:)-P(tri(:,1),:);
d13 = P(tri(:,3),:)-P(tri(:,1),:);
normal = mycross(d12,d13,2);
areaTri = sqrt(sum(normal.^2,2))/2;
area = sum(abs(areaTri));