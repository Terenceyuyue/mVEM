function Pis = faceEllipticProjection(P)
% Compute the matrices of elliptic projection on a face embedded in 3-D
%
% Copyright (C)  Terence Yu.

Nv = size(P,1);
%% Derive local coordinates
% ne,te
e1 = P(2,:)-P(1,:);  en = P(1,:)-P(end,:);
he = sqrt(sum(e1.^2,2));
nf = cross(e1,en); nf = nf./norm(nf);
Ne = -cross(e1, nf);
ne = Ne./repmat(he,1,3);  te = e1./repmat(he,1,3);
% transformation matrix
Tmat = [ne; te];
% node for the 2-D polygon
node = (P - repmat(P(1,:),Nv,1))/Tmat;

%% Get auxiliary data
verts = node; verts1 = verts([2:end,1],:);
area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
ar = 0.5*abs(sum(area_components));
centroid = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
diameter = max(pdist(verts));

%% Compute the elliptic projection matrices
% --- element information (only one element) ---
xK = centroid(1,1); yK = centroid(1,2);
hK = diameter(1);
x = node(:,1); y = node(:,2);  
% --- transition matrix ---
D = [1+0*x, (x-xK)/hK, (y-yK)./hK];   
% --- elliptic projection ---
% first term  = 0
% second term
rotid1 = [Nv,1:Nv-1]; rotid2 = [2:Nv,1]; % ending and starting indices
Gradm = [0 0; 1./hK*[1, 0]; 1./hK*[0, 1]];
normVec = 0.5*[y(rotid2) - y(rotid1), x(rotid1)-x(rotid2)]'; % a rotation of edge vector
B = Gradm*normVec;
% constraint
Bs = B;  Bs(1,:) = 1/Nv;
% consistency relation
Gs = Bs*D;
% elliptic projection matrices
Pis = Gs\Bs;  %Pi = D*Pis;