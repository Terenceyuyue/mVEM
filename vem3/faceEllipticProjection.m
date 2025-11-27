function [Intf,Pifs] = faceEllipticProjection(P,gradmMat)
% Compute the matrices of elliptic projection on a face embedded in 3-D
%
% Copyright (C)  Terence Yu.

nv = size(P,1);
%% Derive local coordinates
% ne,te
e1 = P(2,:)-P(1,:);  en = P(1,:)-P(end,:);
he = sqrt(sum(e1.^2,2));
nf = cross(e1,en); nf = nf./norm(nf);
te = e1./repmat(he,1,3);
ne = cross(nf, te);  % it is the outer normal whenever the vertices are arranges in cyclic order
% transformation matrix
Tmat = [ne; te];
% node for the 2-D polygon
node = (P - repmat(P(1,:),nv,1))/Tmat;

%% Get auxiliary data
verts = node; verts1 = verts([2:end,1],:);
area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
areaf = 0.5*abs(sum(area_components));
centroid = sum((verts+verts1).*repmat(area_components,1,2))/(6*areaf);
diameter = max(pdist(verts));

%% Compute the elliptic projection matrices
% --- element information (only one element) ---
sf = centroid(1,1); tf = centroid(1,2); hf = diameter;
s = node(:,1); t = node(:,2);  

% --- transition matrix ---
Df = [1+0*s, (s-sf)/hf, (t-tf)/hf]; 

% --- elliptic projection ---
% first term  = 0
% second term
rotid1 = [nv,1:nv-1]; rotid2 = [2:nv,1]; % ending and starting indices
gradmf = [0 0; 1./hf*[1, 0]; 1./hf*[0, 1]];
normVec = 0.5*[t(rotid2) - t(rotid1), s(rotid1)-s(rotid2)]'; % a rotation of edge vector
Bf = gradmf*normVec;
% constraint
Bfs = Bf;  Bfs(1,:) = 1/nv;
% consistency relation
Gfs = Bfs*Df;
% elliptic projection matrices
Pifs = Gfs\Bfs;  %Pif = Df*Pifs;

if isempty(gradmMat)
    Intf = 0;
    return;
end

%% Compute the face integral for B
% I1 = 0
% I2: gradm*nf * mf^T
intProj = [areaf,0,0]*Pifs;  % local
Intf = dot(gradmMat, repmat(nf,4,1), 2)*intProj;