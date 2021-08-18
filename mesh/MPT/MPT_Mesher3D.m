function [node,elem,Pn,P] = MPT_Mesher3D(Pb,NT,MaxIter,BdBox)

% ---------- Convex hull approximatation of boundary ---------
Options.pbound = Pb;

% ----------------- Generation of intial pointset -----------------
P = zeros(NT,3); Nbox = length(BdBox); s = 1;
while s<=NT
    x0 = (BdBox(2:2:Nbox)-BdBox(1:2:Nbox-1)).*rand(1,3) + BdBox(1:2:Nbox);
    if isinside(Pb,x0(:))
        P(s,:) = x0;  s = s+1;
    end
end

% -------------------- Lloyd's iteration ---------------
Pn = mpt_voronoi(P,Options);
for iter = 1:MaxIter
    for iel = 1:NT
        V = extreme(Pn(iel)); Tri = convhulln(V);
        P(iel,:) = polycentroid3(V,Tri);
    end
    Pn = mpt_voronoi(P,Options);
end

% ------------------ elementwise nodeElem, elem ------------------
nodeElem = cell(NT,1); elem = cell(NT,1);
for iel = 1:NT
    [V,~,~,adjV] = extreme(Pn(iel));
    [H,K] = double(Pn(iel));  nodeElem{iel} = V;    nf = length(K);
    for i = 1:nf % loop of faces
        % find all points on the current hyperplane
        indexf = find(abs(V*H(i,:)'-K(i))<=1e-6);
        % adjV of the current face
        adjVf = cellfun(@(x)intersect(x,indexf), adjV(indexf),'UniformOutput',false);
        % cyclic order
        nindexf = length(indexf);
        order = zeros(1,nindexf);
        pre_id = 1;  previous = indexf(pre_id);
        adj = adjVf{pre_id}; 
        current = adj(1); cur_id = find(indexf==current);
        order(1) = previous; order(2) = current;
        for p = 2:nindexf-1
            adj = adjVf{cur_id};  next = adj(adj~=previous);
            order(p+1) = next;
            previous = current; current = next;
            cur_id = find(indexf==current);
        end
        % counterclockwise order
        volumn = det([ones(4,1), [V(order(1:3),:); P(iel,:)]])/6;
        if volumn<0, order = order(end:-1:1); end
        elem{iel}{i} = order; % natrual numbers
    end
end

% ---------------- node, elem -----------------
% node
node = vertcat(nodeElem{:}); totalN = size(node,1);
totalid = (1:totalN)';
for i = 1:totalN
    for j = i+1:totalN
        p1 = node(i,:); p2 = node(j,:);
        if norm(p1-p2)<=1e-4 && totalid(j)>i
            totalid(j) = i;
        end
    end
end
[id, ~, totalJ] = unique(totalid,'rows');
node = node(id,:);
% elem
elemLen = cellfun('length',nodeElem);
indexElem = mat2cell(totalJ',1,elemLen)';
for iel = 1:NT
    index = indexElem{iel};
    elemf = cellfun(@(ic) index(ic),elem{iel},'UniformOutput',false);
    elem{iel} = elemf';
end

% ----------------------- Remove small edges ------------------------------
% find short edges
totalfaces = vertcat(elem{:});
T1 = cellfun(@(verts) [verts(2:end),verts(1)], totalfaces, 'UniformOutput', false);
v0 = horzcat(totalfaces{:})'; v1 = horzcat(T1{:})';
totalEdge = sort([v0,v1],2);
[i,j] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:);
he = sqrt(sum((z1-z2).^2,2));
ir = find(he < 0.1*max(he));
nv1 = edge(ir,1); nv2 = edge(ir,2); % starting and ending indices
% replace vertices of short edges
nodeElem = cell(NT,1);
for iel = 1:NT
    elemf = elem{iel};
    [id,~,totalid] = unique(horzcat(elemf{:})');
    totalid = id(totalid);
    for i = 1:length(nv2)
        v1 = nv1(i); v2 = nv2(i); iv = ir(i);
        ax = min(abs([z1(iv,1)-BdBox(1), BdBox(2)-z1(iv,1)]));
        ay = min(abs([z1(iv,2)-BdBox(3), BdBox(4)-z1(iv,2)]));
        az = min(abs([z1(iv,3)-BdBox(5), BdBox(6)-z1(iv,3)]));
        bx = min(abs([z2(iv,1)-BdBox(1), BdBox(2)-z2(iv,1)]));
        by = min(abs([z2(iv,2)-BdBox(3), BdBox(4)-z2(iv,2)]));
        bz = min(abs([z2(iv,3)-BdBox(5), BdBox(6)-z2(iv,3)]));
        if (ax+ay+az)<(bx+by+bz)
            totalid(totalid==v2) = v1; % replaced by starting indices
            nv1(nv1==v2) = v1; nv2(nv2==v2) = v1;
        else
            totalid(totalid==v1) = v2;
            nv1(nv1==v1) = v2; nv2(nv2==v1) = v2;
        end
    end
    elemfLen = cellfun('length',elemf);
    elemf = mat2cell(totalid', 1, elemfLen)';
    for i = 1:length(elemf)
        index = elemf{i};
        [~,i1] = unique(index);
        elemf{i} = index(sort(i1));
    end
    [id,~,totalid] = unique(horzcat(elemf{:})');
    nodeElem{iel} = node(id,:);
    elemfLen = cellfun('length',elemf);
    elemf = mat2cell(totalid', 1, elemfLen)';  % 单元自然序号
    elem{iel} = elemf';
end
% node
node = vertcat(nodeElem{:}); totalN = size(node,1);
totalid = (1:totalN)';
for i = 1:totalN
    for j = i+1:totalN
        p1 = node(i,:); p2 = node(j,:);
        if norm(p1-p2)<=1e-6 && totalid(j)>i
            totalid(j) = i;
        end
    end
end
[~, ~, totalJ] = unique(totalid,'rows');
%node = node(id,:);
% elem
elemLen = cellfun('length',nodeElem);
indexElem = mat2cell(totalJ',1,elemLen)';
for iel = 1:NT
    index = indexElem{iel};
    elemf = cellfun(@(ic) index(ic),elem{iel},'UniformOutput',false);
    elemfLen = cellfun('length',elemf);
    elemf = elemf(elemfLen>=3); % 去掉不构成面的
    elem{iel} = elemf';
end
% Pn with small edges removed
for iel = 1:NT
    Pn(iel) = polytope(nodeElem{iel});
end

% ------ Restore the coplanrity by repeating the previous procedure ------
% elementwise nodeElem, elem 
nodeElem = cell(NT,1); elem = cell(NT,1);
for iel = 1:NT
    [V,~,~,adjV] = extreme(Pn(iel));
    [H,K] = double(Pn(iel));  nodeElem{iel} = V;    nf = length(K);
    for i = 1:nf % loop of faces
        % find all points on the current hyperplane
        indexf = find(abs(V*H(i,:)'-K(i))<=1e-6);
        % adjV of the current face
        adjVf = cellfun(@(x)intersect(x,indexf), adjV(indexf),'UniformOutput',false);
        % cyclic order
        nindexf = length(indexf);
        order = zeros(1,nindexf);
        pre_id = 1;  previous = indexf(pre_id);
        adj = adjVf{pre_id}; 
        current = adj(1); cur_id = find(indexf==current);
        order(1) = previous; order(2) = current;
        for p = 2:nindexf-1
            adj = adjVf{cur_id};  next = adj(adj~=previous);
            order(p+1) = next;
            previous = current; current = next;
            cur_id = find(indexf==current);
        end
        % counterclockwise order
        volumn = det([ones(4,1), [V(order(1:3),:); P(iel,:)]])/6;
        if volumn<0, order = order(end:-1:1); end
        elem{iel}{i} = order; % natrual numbers
    end
end

% node
node = vertcat(nodeElem{:}); totalN = size(node,1);
totalid = (1:totalN)';
for i = 1:totalN
    for j = i+1:totalN
        p1 = node(i,:); p2 = node(j,:);
        if norm(p1-p2)<=1e-4 && totalid(j)>i
            totalid(j) = i;
        end
    end
end
[id, ~, totalJ] = unique(totalid,'rows'); node = node(id,:);
% elem
elemLen = cellfun('length',nodeElem);
indexElem = mat2cell(totalJ',1,elemLen)';
for iel = 1:NT
    index = indexElem{iel};
    elemf = cellfun(@(ic) index(ic),elem{iel},'UniformOutput',false);
    elem{iel} = elemf';
end