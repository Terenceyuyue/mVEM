function poly = localPolygon3(P)

Nv = size(P,1);
%% Derive local coordinates
% ne,te
e1 = P(2,:)-P(1,:);
en = P(1,:)-P(end,:);
he = sqrt(sum(e1.^2,2));
nf = cross(e1,en); nf = nf./norm(nf);
Ne = -cross(e1, nf);
ne = Ne./repmat(he,1,3);
te = e1./repmat(he,1,3);
% transformation matrix
Tmat = [ne; te];
% node, elem for the 2-D polygon
node = (P - repmat(P(1,:),Nv,1))/Tmat;
%elem = {1:Nv};

%% Get auxiliary data
verts = node; verts1 = verts([2:end,1],:);
area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
area = 0.5*abs(sum(area_components));
centroid = sum((verts+verts1).*repmat(area_components,1,2))/(6*area);
diameter = max(pdist(verts));

%% Store info
poly.nodef = node;  
poly.Coord = @(s,t) s*ne + t*te + P(1,:); % (s,t) --> (x,y,z)
poly.areaf = area;  
poly.centroidf = centroid; % local coordinates (s,t)
poly.centroid = centroid(1)*ne + centroid(2)*te + P(1,:); % global coordinates (x,y,z)
poly.diameterf = diameter;
poly.nf = nf;