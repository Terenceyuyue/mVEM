function [pv,ev] = removeEdges(pv,ev,op)
% remove overly short edges by merging adj. nodes.

% calc. edge lengths
el = sum((pv(ev(:,2),:)-pv(ev(:,1),:)).^2,2);

% calc. cell measure
nc = max(max(ev(:,3:4)));
cl = zeros(nc,1);
for ei = 1 : size(ev,1)
    for ii = 3 : 4
        % adj. cell
        ci = ev(ei,ii) ;
        % null cell
        if (ci <= +1), continue; end
        % keep max so far
        if (el(ei) > cl(ci)), cl(ci) = el(ei); end
    end
end
cl = cl * op.etol^2;

% collapse small edges/merge nodes

ki = zeros(size(pv,1),1);
for ei = 1 : size(ev,1)
    % adj. node
    ni = ev(ei,1);    nj = ev(ei,2);
    % chase node re-indexing
    while (ki(ni) > +0), ni = ki(ni); end
    while (ki(nj) > +0), nj = ki(nj); end
    % skip already collapsed
    if (ni == nj), continue; end
    % scan adj. cells
    for ii = 3 : 4
        % adj. cell
        ci = ev(ei,ii) ;
        % null cell
        if (ci <= +1), continue; end
        % collapse if sufficiently small
        ll = sum((pv(nj,:)-pv(ni,:)).^2,2);
        if (ll < cl(ci))
            % merge node position (È¡ÖÐµã)
            pv(ni,:) = (pv(ni,:)+pv(nj,:))/2;
            % merge node indexing
            ki(nj) = ni; break ;
        end
    end
end

% remove any collapsed edges
ne = +1;
for ei = 1 : size(ev,1)
    % adj. node
    ni = ev(ei,1); nj = ev(ei,2);
    % adj. cell
    ci = ev(ei,3); cj = ev(ei,4);
    % chase node re-indexing
    while (ki(ni) > +0), ni = ki(ni); end
    while (ki(nj) > +0), nj = ki(nj); end
    % keep un-collapsed edge
    if (ni ~= nj)
        ev(ne,1) = ni; ev(ne,2) = nj;
        ev(ne,3) = ci; ev(ne,4) = cj;
        ne = ne + 1;
    end
end
ev = ev(1:ne-1,:);