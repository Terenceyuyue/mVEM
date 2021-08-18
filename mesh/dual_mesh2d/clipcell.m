function [cp,ce] = clipcell(ev)
% build dual cells - clip about non-manifold edges.

% triangulate dual cells
tv = zeros(size(ev,1)*2,3); nt = +1 ;
for ei = 1 : size(ev,1)
    for ii = 3 : 4
        % cell node
        ni = ev(ei,ii);
        % skip if constraint
        if (ni <= +0 ), continue; end
        % edge node
        nj = ev(ei,1); nk = ev(ei,2);
        % skip if degenerate
        if (ni == nj || ni == nk ), continue; end
        % push tria
        tv(nt,1) = ni; tv(nt,2) = nj; tv(nt,3) = nk;
        nt = nt+1;
    end
end
tv = tv(1:nt-1,:);

%------------------------------------------ edge-tria conn.
[ee,te,ep,et] = triaconn2(tv);

%------------------------------------------ mark cell edges
[ef,em] = ismember(sort(ee,2),sort(ev(:,1:2),2),'rows') ;

%----------------------------------- assemble clipped cells
cp = zeros(size(tv,1)*1,3); nc = +1;
ce = zeros(size(tv,1)*1,1); ne = +1;
kt = zeros(size(tv,1)*1,1);
ts = zeros(size(tv,1)*1,1);
for ti = 1 : size(tv)
    if (kt(ti) ~= +0), continue; end
    %-------------------------------------------- cell node
    ni = tv(ti,1);
    %-------------------------------------------- cell ptrs
    cp(nc,1) = ni;
    cp(nc,2) = ne;
    %-------------------------------- assemble cell via bfs
    ns = +1;
    ts(ns) = ti;
    kt(ti) = +1;
    while (ns ~= +0)
        %----------------------------------- _pop bfs stack
        tj = ts(ns);
        ns = ns - 1;
        %----------------------------------- scan tria adj.
        for ei = 1 : 3
            ej = te(tj,ei) ;
            if (em(ej)~=0)
                %---------------------- push external cell edge
                ce(ne) = em(ej); ne = ne+1;
            else
                %----------------------------- must be manifold
                if (tj ~= et(ep(ej,1)))
                    tk  = et(ep(ej,1));
                else
                    tk  = et(ep(ej,2));
                end
                %----------------------------- walk to adjacent
                if (kt(tk)==0)
                    ns = ns+1;
                    ts(ns)=tk;
                    kt(tk)=+1;
                end
            end
        end
    end
    %------------------------------------------ update ptrs
    cp(nc,3) = ne - 1; nc = nc + 1;
end
%---------------------------------------------- trim alloc.
cp = cp(1:nc-1,:) ;
ce = ce(1:ne-1,:) ;