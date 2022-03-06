function centroid = polycentroid(V)
% polycentroid returns the x,y coordinates of centroid of polygon in 2-D
% The vertices must be in cyclic order
%
%    V = [x,y]
%
% Copyright (C) Terence Yu.

V1 = circshift(V,-1);
area_components = V(:,1).*V1(:,2) - V1(:,1).*V(:,2);
ar = 0.5*(sum(area_components)); 
centroid = sum((V+V1).*repmat(area_components,1,2))/(6*ar);


