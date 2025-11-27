function aux = auxgeometry(node,elem,isTri)
%auxgeometry gets geometry data
%
% Copyright (C) Terence Yu.

if nargin==2
    isTri = false;
end

%% Geometric quantities
NT = size(elem,1);
if ~iscell(elem) % transform to cell
    elem = mat2cell(elem,ones(NT,1),length(elem(1,:)));
end
% ------ centroid, area, diameter -------
centroid = zeros(NT,2); area = zeros(NT,1); diameter = zeros(NT,1);
s = 1;
for iel = 1:NT
    index = elem{iel};
    verts = node(index, :); verts1 = circshift(verts,-1); % verts([2:end,1],:);
    area_components = verts(:,1).*verts1(:,2)-verts1(:,1).*verts(:,2);
    ar = 0.5*abs(sum(area_components));
    area(iel) = ar;
    centroid(s,:) = sum((verts+verts1).*repmat(area_components,1,2))/(6*ar);
    diameter(s) = max(pdist(verts));    
    s = s+1;
end

%% Triangulation of elements
if isTri
    nodeTri = cell(1,NT); elemTri = cell(1,NT);
    for iel = 1:NT
        % element info
        index = elem{iel};   Nv = length(index);
        P1 = node(index,:);  P2 = node(index([2:Nv,1]),:);
        % Decomposed geometry matrix
        N = size(P1,1);
        g = zeros(7,N);
        g(1,:) = 2;  % 2 for line
        g(2,:) = P1(:,1); % x: starting points
        g(3,:) = P2(:,1); % x: ending points
        g(4,:) = P1(:,2); % y: starting points
        g(5,:) = P2(:,2); % y: ending points
        g(6,:) = 1; % subdomain on the left
        g(7,:) = 0; % subdomain on the right
        % triangulation
        [p,~,t] = initmesh(g);
        nodeTri{iel} = p';
        elemTri{iel} = t(1:3,:)';
    end
end

%% Store the info
aux.node = node; aux.elem = elem;
aux.centroid = centroid;
aux.area = area;
aux.diameter = diameter;

if isTri
    aux.nodeTri = nodeTri;
    aux.elemTri = elemTri;
end
