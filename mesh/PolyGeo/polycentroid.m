function centroid = polycentroid(V)
% polycentroid returns the coordinates of the centroid of a polygon in 2-D
% or a polygonal face embedded in 3-D.
% The vertices must be in cyclic order
%
% Copyright (C) Terence Yu.

if size(V,2)==2
    V1 = circshift(V,-1);
    area_components = V(:,1).*V1(:,2) - V1(:,1).*V(:,2);
    ar = 0.5*(sum(area_components));
    centroid = sum((V+V1).*repmat(area_components,1,2))/(6*ar);
    return;
end

% if size(V,2)==3  % 3-D face
Nv = size(V,1);
% ---- derive local coordinates -----
% ne,te
e1 = V(2,:)-V(1,:);
en = V(1,:)-V(end,:);
he = sqrt(sum(e1.^2,2));
nf = cross(e1,en); nf = nf./norm(nf);
Ne = -cross(e1, nf);
ne = Ne./repmat(he,1,3);
te = e1./repmat(he,1,3);
% transformation matrix
Tmat = [ne; te];
% node, elem for the 2-D polygon
node = (V - repmat(V(1,:),Nv,1))/Tmat;
%elem = {1:Nv};
% ----- get auxiliary data -----
verts = node; verts1 = verts([2:end,1],:);
area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
area = 0.5*abs(sum(area_components));
centroidf = sum((verts+verts1).*repmat(area_components,1,2))/(6*area);
centroid = centroidf(1)*ne + centroidf(2)*te + V(1,:); % global coordinates (x,y,z)