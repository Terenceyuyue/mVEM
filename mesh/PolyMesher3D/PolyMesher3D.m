function [node,elem] = PolyMesher3D(Domain,MaxIter,NT,rm)

% ------------------- Generate intial pointset ----------------------
P = PolyMesher_init_Pointset3D(Domain,NT);

% ------------------------ Lloyd iteration --------------------------
BdBox = Domain('BdBox');
volume = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))*(BdBox(6)-BdBox(5));
elemTri = cell(NT,1); iter = 1;
while(iter<=MaxIter)
    % Compute the reflection pointset
    R_P = PolyMesher_Reflect3D(P,Domain,volume);
    % Construct the Voronoi diagram
    [V,C] = voronoin([P;R_P],{'Qbb','Qz'});
    % Compute the centroid of Voronoi cell (Lloyd's update)
    P = zeros(NT,3); volume = 0;
    for iel = 1:NT
        index = C{iel};
        X = V(index,:); [Tri,vol] = convhulln(X);
        P(iel,:) = polycentroid3(X,Tri); volume = volume+vol;
        if iter==MaxIter, elemTri{iel} = index(Tri); end
    end
    iter = iter+1;
end

% ------------ Remove small edges in triangulation -------------
cr = 1e-6;
% find short edges
totalTri = vertcat(elemTri{:}); % NTri * 3
totalTri = sort(totalTri,2);
Tri = unique(totalTri,'rows');
totalEdge = sort([Tri(:,[2,3]); Tri(:,[3,1]); Tri(:,[1,2])],2);
[i,j] = find(sparse(totalEdge(:,2),totalEdge(:,1) ,1));
edge = [j,i];
% replace vertices of short edges
z1 = V(edge(:,1),:); z2 = V(edge(:,2),:);
he = sqrt(sum((z1-z2).^2,2));
if rm
    ir = find(he < 1e-6);
else
    ir = find(he < cr*max(he));
end
nv1 = edge(ir,1); nv2 = edge(ir,2); % starting and ending indices
for i = 1:length(nv2)
    v1 = nv1(i); v2 = nv2(i); iv = ir(i);
    ax = min(abs([z1(iv,1)-BdBox(1), BdBox(2)-z1(iv,1)]));
    ay = min(abs([z1(iv,2)-BdBox(3), BdBox(4)-z1(iv,2)]));
    az = min(abs([z1(iv,3)-BdBox(5), BdBox(6)-z1(iv,3)]));
    bx = min(abs([z2(iv,1)-BdBox(1), BdBox(2)-z2(iv,1)]));
    by = min(abs([z2(iv,2)-BdBox(3), BdBox(4)-z2(iv,2)]));
    bz = min(abs([z2(iv,3)-BdBox(5), BdBox(6)-z2(iv,3)]));
    if (ax+ay+az)<(bx+by+bz)
        totalTri(totalTri==v2) = v1; % replaced by starting indices
        nv1(nv1==v2) = v1;  nv2(nv2==v2) = v1;
    else
        totalTri(totalTri==v1) = v2;
        nv1(nv1==v1) = v2;  nv2(nv2==v1) = v2;
    end
end
% elemTri
elemTriLen = cellfun('length',elemTri);
elemTri = mat2cell(totalTri, elemTriLen, 3);
for iel = 1:NT
    index = elemTri{iel};  nindex = size(index,1);
    rows = true(nindex,1);
    for i = 1:size(index,1)
        id = index(i,:); nid = length(unique(id));
        if nid<3, rows(i) = false; end
    end
    elemTri{iel} = index(rows,:);
end
% node, elemTri (remove redundant nodes)
totalTri = vertcat(elemTri{:});  % NTri * 3
totalTri = sort(totalTri,2);
totalTri = totalTri'; totalTri = totalTri(:); % 按行拉直
[id,~,totalid] = unique(totalTri);
node = V(id,:);
totalTri = reshape(totalid,3,[])';
elemTriLen = cellfun('length',elemTri);
elemTri = mat2cell(totalTri, elemTriLen, 3);
% restore the coplanarity
for iel = 1:NT
    index = unique(elemTri{iel});
    X = node(index,:); Tri = convhulln(X);
    elemTri{iel} = index(Tri);
end

% --------- Determine elemf consisting of coplanar triangles -------------
cp = 1e-5; % coplanar parameter
elem = cell(NT,1);
for iel = 1:NT
    Tri = elemTri{iel}; elemf = []; s = 1;
    while size(Tri,1)>=1
        rows = false(size(Tri,1),1); rows(1) = true;
        T = Tri(1,:);  % current triangle
        % indexf
        index = unique(Tri); nindex = length(index); % 单元顶点编号
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
        elemf{s} = Tri(rows,:);
        s = s+1;
        Tri = Tri(~rows,:);  % triangles remaining
    end
    elem{iel} = elemf';
end

% ------------------------ elem -----------------------------
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