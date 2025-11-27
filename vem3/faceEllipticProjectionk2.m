function Intf = faceEllipticProjectionk2(P,gradmc)
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
% gemoetric quantities
verts = node; verts1 = verts([2:end,1],:);
area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
area = 0.5*abs(sum(area_components));
centroid = sum((verts+verts1).*repmat(area_components,1,2))/(6*area);
diameter = max(pdist(verts));
% number
Nmf = 6;

%% Compute the elliptic projection matrices on face
% --- element information (only one element) ---
sf = centroid(1,1); tf = centroid(1,2);  hf = diameter;
s = node(:,1); t = node(:,2); 
v1 = 1:nv; v2 = [2:nv,1]; % loop index for vertices or edges
s0 = (s(v1)+s(v2))/2;  y0 = (t(v1)+t(v2))/2; % mid-edge points
Ne = [t(v2)-t(v1), s(v1)-s(v2)]; % he*ne
nodeT = [node;centroid];
elemT = [(nv+1)*ones(nv,1),(1:nv)',[2:nv,1]'];

% ------------ scaled monomials --------
mf1 = @(s,t) 1+0*s;                   gradmf1 = @(s,t) [0+0*s, 0+0*s];
mf2 = @(s,t) (s-sf)/hf;               gradmf2 = @(s,t) [1+0*s, 0+0*s]/hf;
mf3 = @(s,t) (t-tf)/hf;               gradmf3 = @(s,t) [0+0*s, 1+0*s]/hf;
mf4 = @(s,t) (s-sf).^2/hf^2;          gradmf4 = @(s,t) [2*(s-sf), 0+0*s]/hf^2;
mf5 = @(s,t) (s-sf).*(t-tf)/hf^2;     gradmf5 = @(s,t) [(t-tf), (s-sf)]/hf^2;
mf6 = @(s,t) (t-tf).^2/hf^2;          gradmf6 = @(s,t) [0+0*s, 2*(t-tf)]/hf^2;
mf = @(s,t) [mf1(s,t), mf2(s,t), mf3(s,t), mf4(s,t), mf5(s,t), mf6(s,t)];
gradmfc = {gradmf1, gradmf2, gradmf3, gradmf4, gradmf5, gradmf6};
    
% -------- transition matrix ----------
ndof = 2*nv+1;
Df = zeros(ndof,Nmf);
Df(1:2*nv,:) = [mf(s,t); mf(s0,y0)];  % vertices and mid-edge points
Df(end,:) = 1/area*integralTri(mf,5,nodeT,elemT);

% --------- elliptic projection -----------
% first term
Lapmf = zeros(6,1); Lapmf([4,6]) = 2/hf^2;
Dof = [zeros(1,2*nv), area];
I1 = Lapmf*Dof;
% second term
I2 = zeros(Nmf, ndof);
elem1 = [v1(:), v2(:), v1(:)+nv];
for im = 1:Nmf
    gradmf = gradmfc{im};
    F1 = 1/6*sum(gradmf(s(v1), t(v1)).*Ne, 2);
    F2 = 1/6*sum(gradmf(s(v2), t(v2)).*Ne, 2);
    F3 = 4/6*sum(gradmf(s0, y0).*Ne, 2);
    F = [F1, F2, F3];
    I2(im, :) = accumarray(elem1(:), F(:), [ndof 1]);
end
Bf = -I1 + I2;
% constraint
Bfs = Bf;  Bfs(1,:) = Dof;
% consistency relation
Gfs = Bfs*Df;
% elliptic projection matrices
Pifs = Gfs\Bfs;  %Pif = Df*Pifs;

%% Compute the face integral for B
% triangulation
tri = [ones(nv-2,1), (2:nv-1)', (3:nv)'];
% Gauss quadrature points and weights
[lambda,weight] = quadpts(5);
% area of triangles
d12 = P(tri(:,2),:)-P(tri(:,1),:);
d13 = P(tri(:,3),:)-P(tri(:,1),:);
normal = mycross(d12,d13,2);
areaTri = sqrt(sum(normal.^2,2))/2;
% integral
Nm = length(gradmc);  
Dmm = zeros(Nm,Nmf); % gradm*nf * mf^T
for im = 1:Nm
    Int = zeros(nv-2,Nmf);
    for p = 1:length(weight)
        % quadrature points in (x,y,z)
        pxyz = lambda(p,1)*P(tri(:,1),:) ...
            + lambda(p,2)*P(tri(:,2),:) ...
            + lambda(p,3)*P(tri(:,3),:);
        % quadrature points in (s,t)
        pst = lambda(p,1)*node(tri(:,1),:) ...
            + lambda(p,2)*node(tri(:,2),:) ...
            + lambda(p,3)*node(tri(:,3),:);
        f1 = gradmc{im}(pxyz(:,1), pxyz(:,2), pxyz(:,3))*nf(:);
        f2 = mf(pst(:,1),pst(:,2));
        Int = Int + weight(p)*repmat(f1,1,Nmf).*f2;
    end
    Dmm(im,:) = sum(repmat(areaTri,1,Nmf).*Int);    
end
Intf = Dmm*Pifs;