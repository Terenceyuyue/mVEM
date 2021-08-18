function [node,elem] = PolyMesherBd(Pb,maxIter,varargin)
% This function generates polygonal/polyhedral meshes using a very different
% method from PolyMesher.m
% 
% Copyright (C) Terence Yu

d = size(Pb,2);
nvar = length(varargin);
if nvar==1, NT = varargin{1}; end
if nvar==2, Nx = varargin{1}; Ny = varargin{2}; end
if nvar==3, Nx = varargin{1}; Ny = varargin{2}; Nz = varargin{3}; end

%% Generate random seeds
tol = 1e-07;
if d == 2    
    BdBox = [min(Pb(:,1)),max(Pb(:,1)),min(Pb(:,2)),max(Pb(:,2))];
    if nvar==1
        P = zeros(NT,2);  s = 0;
        while s < NT
            p(:,1) = (BdBox(2)-BdBox(1))*rand(NT,1)+BdBox(1);
            p(:,2) = (BdBox(4)-BdBox(3))*rand(NT,1)+BdBox(3);
            in = inhull(p,Pb,[],tol);
            NumAdded = min(NT-s,length(in));  % number of seeds that can be added
            P(s+1:s+NumAdded,:) = p(in(1:NumAdded),:);
            s = s+NumAdded;
        end
    end
    if nvar==2
        x = linspace(BdBox(1),BdBox(2),Nx+1)';
        y = linspace(BdBox(3),BdBox(4),Ny+1)';
        xc = (x(1:end-1)+x(2:end))/2; yc = (y(1:end-1)+y(2:end))/2;
        [X,Y] = ndgrid(xc,yc); P = [X(:),Y(:)];
        in = inhull(P,Pb,[],tol);
        P = P(in,:);
        NT = size(P,1);
    end
end

if d == 3
  BdBox = [min(Pb(:,1)),max(Pb(:,1)),...
           min(Pb(:,2)),max(Pb(:,2)),...
           min(Pb(:,3)),max(Pb(:,3))];
    if nvar==1
        P = zeros(NT,3);  s = 0;
        while s < NT
            p(:,1) = (BdBox(2)-BdBox(1))*rand(NT,1)+BdBox(1);
            p(:,2) = (BdBox(4)-BdBox(3))*rand(NT,1)+BdBox(3);
            p(:,3) = (BdBox(6)-BdBox(5))*rand(NT,1)+BdBox(5);
            in = inhull(p,Pb,[],tol);
            NumAdded = min(NT-s,sum(in));  % number of seeds that can be added
            Pin = p(in,:);
            P(s+1:s+NumAdded,:) = Pin(1:NumAdded,:);
            s = s+NumAdded;
        end
    end
    if nvar==3
        x = linspace(BdBox(1),BdBox(2),Nx+1)';
        y = linspace(BdBox(3),BdBox(4),Ny+1)';
        z = linspace(BdBox(5),BdBox(6),Nz+1)';
        xc = (x(1:end-1)+x(2:end))/2; 
        yc = (y(1:end-1)+y(2:end))/2;
        zc = (z(1:end-1)+z(2:end))/2;
        [X,Y,Z] = ndgrid(xc,yc,zc); P = [X(:),Y(:),Z(:)];
        in = inhull(P,Pb,[],tol);
        P = P(in,:);
        NT = size(P,1);
    end   
end


%% Lloyd iteration for 2D and return [node,elem]
Iter = 0; 
if d == 2
    Err = 1; Tol = 1e-4;
    while (Iter<=maxIter && Err>Tol)
        % Construct the Voronoi diagram
        [node,elem] = voronoinBd(P,Pb);
        % Compute the centroid of Voronoi cell (Lloyd's update)
        [P,~,Err] = PolyMesher_VoroCentroid(P,node,elem);
        Iter = Iter+1;
        fprintf('Iter: %3d   Error: %1.3e\n',Iter,Err);
    end
    % remove small edges
    [node,elem] = PolyMesher_rm_smalledge(node,elem);
    if NT<=1000
        clf; showmesh(node,elem); 
    end
    return;
end

% ------ The remaining part is only for 3D case -------- %

%% Lloyd iteration for 3D
elemTri = cell(NT,1);
while(Iter<=maxIter)
    [node,C] = voronoinBd(P,Pb);
    P = zeros(NT,d);
    for iel = 1:NT
        index = C{iel};  X = node(index,:);
        Tri = convhulln(X);
        P(iel,:) = polycentroid3(X,Tri);
        if Iter==maxIter, elemTri{iel} = index(Tri); end
    end
    Iter = Iter+1;
end

%% Determine elemf consisting of coplanar triangles
cp = 1e-5; % coplanar parameter
elem = cell(NT,1);
for iel = 1:NT
    Tri = elemTri{iel};  nTri = size(Tri,1);  s = 0;
    elemf = cell(nTri,1);
    while size(Tri,1)>=1
        rows = false(size(Tri,1),1); rows(1) = true;
        T = Tri(1,:);  % current triangle
        % indexf
        index = unique(Tri); nindex = length(index); % µ¥Ôª¶¥µã±àºÅ
        z1 = node(T(1),:); z2 = node(T(2),:); z3 = node(T(3),:);
        nvec = cross(z2-z1, z3-z1);
        nvec = nvec/norm(nvec); % normalization
        pvec = node(index,:)-repmat(z1,nindex,1);
        pvec = pvec/norm(pvec); % normalization
        on = abs(dot(repmat(nvec,nindex,1)',pvec'))<=cp;
        indexf = index(on);
        % elemf
        for i = 2:size(Tri,1)
            T1 = Tri(i,:);
            con = intersect(T1,indexf); ncon = length(con);
            if ncon==3, rows(i) = true; end
        end
        s = s+1;
        elemf{s} = Tri(rows,:);
        Tri = Tri(~rows,:);  % triangles remaining
    end
    elem{iel} = elemf(1:s);
end

%% Derive data structure elem
for iel = 1:NT
    elemf = elem{iel};
    for s = 1:size(elemf,1)
        % face_s
        face = elemf{s};
        % boundary edges
        totalEdge = sort([face(:,[2,3]); face(:,[3,1]); face(:,[1,2])],2);
        [i,j,sb] = find(sparse(totalEdge(:,2),totalEdge(:,1) ,1));
        edge = [j,i]; edge = edge(sb==1,:);
        % adjVf
        indexf = unique(edge); nindexf = length(indexf);
        repEdge = [edge; edge(:,[2,1])];
        repEdge = sortrows(repEdge,1);
        adjVf = reshape(repEdge(:,2),2,[])';
        % cyclic order
        order = zeros(1,nindexf);
        pre_id = 1;  previous = indexf(pre_id);
        adj = adjVf(pre_id,:);
        current = adj(1); cur_id = find(indexf==current);
        order(1) = previous; order(2) = current;
        for p = 2:nindexf-1
            adj = adjVf(cur_id,:);  next = adj(adj~=previous);
            order(p+1) = next;
            previous = current; current = next;
            cur_id = find(indexf==current);
        end
        % counterclockwise order
        volumn = det([ones(4,1), [node(order(1:3),:); P(iel,:)]])/6;
        if volumn<0, order = order(end:-1:1); end
        elemf{s} = order;
    end
    elem{iel} = elemf;
end