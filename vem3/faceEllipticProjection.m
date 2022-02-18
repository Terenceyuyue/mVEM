function Pis = faceEllipticProjection(P)
% Compute the matrices of elliptic projection on a face embedded in 3-D
%
% Copyright (C)  Terence Yu.

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
Nm = 3;

%% Compute the elliptic projection matrices
% element information (only one element)
index = elem{1};  
xK = centroid(1,1); yK = centroid(1,2);
hK = diameter(1);
x = node(index,1); y = node(index,2);
v1 = 1:Nv;  v2 = [2:Nv,1];
Ne = [y(v2)-y(v1), x(v1)-x(v2)]; % he*ne

% scaled monomials
m1 = @(x,y) 1+0*x;                gradm1 = @(x,y) [0+0*x, 0+0*x];
m2 = @(x,y) (x-xK)./hK;           gradm2 = @(x,y) [1+0*x, 0+0*x]./hK;
m3 = @(x,y) (y-yK)./hK;           gradm3 = @(x,y) [0+0*x, 1+0*x]./hK;
m = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];
Gradmc = {gradm1, gradm2, gradm3};

% transition matrix
D = m(x,y);   

% ----- elliptic projection ------
% first term  = 0
B = zeros(Nm, Nv);
% second term
elem1 = [v1(:), v2(:)];
for im = 1:Nm
    gradmc = Gradmc{im};
    F1 = 0.5*sum(gradmc(x(v1), y(v1)).*Ne, 2);
    F2 = 0.5*sum(gradmc(x(v2), y(v2)).*Ne, 2);
    F = [F1, F2];
    B(im, :) = accumarray(elem1(:), F(:), [Nv 1]);
end
% constraint
Bs = B;  Bs(1,:) = 1/Nv;
% consistency relation
Gs = Bs*D;
% elliptic projection matrices
Pis = Gs\Bs;  %Pi = D*Pis;