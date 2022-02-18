function poly = localPolygon3(P)

Nv = size(P,1);
%% Derive local coordinates
% ne,te
e1 = P(2,:)-P(1,:);
en = P(1,:)-P(end,:);
he = sqrt(sum(e1.^2,2));
nf = mycross(e1,en); nf = nf./norm(nf);
Ne = -mycross(e1, nf);
ne = Ne./repmat(he,1,3);
te = e1./repmat(he,1,3);
% transformation matrix
Tmat = [ne; te];
% node, elem for the 2-D polygon
node = (P - repmat(P(1,:),Nv,1))/Tmat;
elem = {1:Nv};

%% Get auxiliary data
aux = auxgeometry(node,elem);
centroid = aux.centroid;  
diameter = aux.diameter; 
area = aux.area;

%% Store info
poly.nodef = node;  
poly.Coord = @(s,t) s*ne + t*te + P(1,:); % (s,t) --> (x,y,z)
poly.areaf = area;  
poly.centroidf = centroid; 
poly.diameterf = diameter;
poly.nf = nf;