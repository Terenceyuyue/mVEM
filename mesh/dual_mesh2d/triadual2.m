function [tp,tv] = triadual2(cp,ce,ev)
%TRIADUAL2 triangulate dual mesh cells.
%   [TP,TV] = TRIADUAL2(CP,CE,EV) returns a conforming tria-
%   ngulation of the dual cells created by MAKEDUAL2. The
%   triangles associated with the dual cell CI are given by 
%   TV(TP(CI,1):TP(CI,2),:).

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 29/11/2014

%--------------------------------------- count tria indexing
    ne = cp(:,3)-cp(:,2)+1; nt = +1;
%--------------------------------------- alloc tria indexing
    tp = zeros(size(cp,1),2); 
    tv = zeros( sum(ne)*1,3);
    for ci = 1 : size(cp,1)        
        tp(ci,1) = nt+0;
    %----------------------------------------- "centre" node
        ni = cp(ci,1);
    %------------------ triangulate cell about "centre" node
        for ii = cp(ci,2):cp(ci,3)
        %----------------------------------------- adj. edge
            ei = ce(ii,1) ;
        %----------------------------------------- adj. node
            nj = ev(ei,1) ;
            nk = ev(ei,2) ;
        %-------------------------------- skip if degenerate
            if (ni == nj || ...
                ni == nk ), continue; end
        %---------------------------------------- local tria
            tv(nt,1) = ni ;
            tv(nt,2) = nj ;
            tv(nt,3) = nk ;
            nt = nt+1;
        end        
        tp(ci,2) = nt-1;
    end
%------------------------------------------------ trim alloc
    tv = tv(1:nt-1,:);

end
