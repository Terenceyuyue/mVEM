function [node,elem] = MPT_Mesher2D(Pb,NT,MaxIter,BdBox)

% ---------- Convex hull approximatation of boundary ---------
Options.pbound = Pb;

% ----------------- Generation of intial pointset -----------------
P = zeros(NT,2); Nbox = length(BdBox); s = 1;
while s<=NT
    x0 = (BdBox(2:2:Nbox)-BdBox(1:2:Nbox-1)).*rand(1,2) + BdBox(1:2:Nbox);
    if isinside(Pb,x0(:))
        P(s,:) = x0;  s = s+1;
    end
end

% -------------------- Lloyd's iteration ---------------
Pn = mpt_voronoi(P,Options); nodeElem = cell(NT,1); iter = 1;
while (iter<=MaxIter)
    for iel = 1:NT
        % cyclic order
        [V,~,~,adjV] = extreme(Pn(iel));
        nv = size(V,1);   order = zeros(nv,1);
        previous = 1; adj = adjV{previous};
        current = adj(1);
        order(1) = previous; order(2) = current;
        for p = 2:nv-1
            adj = adjV{current};  next = adj(adj~=previous);
            order(p+1) = next;
            previous = current; current = next;
        end
        V = V(order,:); nodeElem{iel} = V;
        % centroid
        P(iel,:) = polycentroid(V);
    end    
    Pn = mpt_voronoi(P,Options);     iter = iter+1;
end

% ---------- Counterclockwise order of elements in Pn -----------
elemLen = cellfun('length',double(Pn)); totalN = sum(elemLen);
node = zeros(totalN,2); s = 1;
for iel = 1:NT  
    V = nodeElem{iel}; nV = size(V,1);
    z1 = V(1,:); z2 = V(2,:); z3 = V(3,:);
    e2 = z3-z1; e3 = z1-z2;
    area = 0.5*(-e3(:,1).*e2(:,2)+e3(:,2).*e2(:,1));
    if area<0, V = V(end:-1:1,:); end
    node(s:s+nV-1,:) = V;  s = s+nV;
end

% ---------------------- node, elem ----------------------
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
elem = mat2cell(totalJ',1,elemLen)';

% ---------------------- Remove small edges --------------------------
% tol = 5e-3;
% for iel = 1:1
%     index = elem{iel};  nv = length(index);
%     if nv<4, continue; end  % cannot collapse triangles
%     z1 = node(index, :); z2 = z1([2:end,1],:);
%     area_components = z1(:,1).*z2(:,2)-z2(:,1).*z1(:,2);
%     ar = 0.5*abs(sum(area_components));
%     zc = sum((z1+z2).*repmat(area_components,1,2))/(6*ar);
%     ea = z1-zc; eb = z2-zc; ec = z1-z2;
%     a2 = sum(ea.^2,2); b2 = sum(eb.^2,2); c2 = sum(ec.^2,2);
%     a = sqrt(a2); b = sqrt(b2); c = sqrt(c2);
%     beta = acos( (a2+b2-c2)./(2*a.*b) );
%     betaId = find(beta<tol*2*pi/nv)
%     
% end


T1 = cellfun(@(verts) [verts(2:end),verts(1)], elem, 'UniformOutput', false);
v0 = horzcat(elem{:})'; v1 = horzcat(T1{:})';
totalEdge = sort([v0,v1],2);
[i,j] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:);
he = sqrt(sum((z1-z2).^2,2));
ir = find(he < 0.1*max(he));
nv1 = edge(ir,1); nv2 = edge(ir,2); % starting and ending indices
[id,~,totalid] = unique(horzcat(elem{:})');
totalid = id(totalid); % actually not needed since id = 1,2,...
for i = 1:length(nv2)
    v1 = nv1(i); v2 = nv2(i); iv = ir(i);
    ax = min(abs([z1(iv,1)-BdBox(1), BdBox(2)-z1(iv,1)]));
    ay = min(abs([z1(iv,2)-BdBox(3), BdBox(4)-z1(iv,2)]));
    bx = min(abs([z2(iv,1)-BdBox(1), BdBox(2)-z2(iv,1)]));
    by = min(abs([z2(iv,2)-BdBox(3), BdBox(4)-z2(iv,2)]));
    if (ax+ay)<(bx+by)
        totalid(totalid==v2) = v1; % replaced by starting indices
        nv1(nv1==v2) = v1;
        nv2(nv2==v2) = v1;
    else
        totalid(totalid==v1) = v2;
        nv1(nv1==v1) = v2;
        nv2(nv2==v1) = v2;
    end
end
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';
for iel = 1:NT
    index = elem{iel};
    [~,i1] = unique(index);
    elem{iel} = index(sort(i1));
end
[id,~,totalid] = unique(horzcat(elem{:})');
node = node(id,:);
elemLen = cellfun('length',elem);
elem = mat2cell(totalid', 1, elemLen)';