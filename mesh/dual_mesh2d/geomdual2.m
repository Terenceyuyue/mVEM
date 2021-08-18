function [pc,ac] = geomdual2(cp,ce,pv,ev)
%GEOMDUAL2 calc. geometric properties for a dual mesh.
%   [PC,AC] = GEOMDUAL2(CP,CE,PV,EV) returns the cell bary-
%   centres and areas for the dual mesh constructed using
%   MAKEDUAL2.

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 29/11/2014

%----------------------------- alloc. cell areas/barycentres
    [ac] = zeros(size(cp,1),1) ;
    [pc] = zeros(size(cp,1),3) ;
%------------------------------------ triangulate dual cells
    [tp,tv] = triadual2(cp,ce,ev) ;
%------------------------------------------------ tria norms
    [nt] = cross(pv(tv(:,2),:)-pv(tv(:,1),:), ...
                 pv(tv(:,3),:)-pv(tv(:,1),:)) ;
%------------------------------------------------ tria areas
    [at] = +.5 * sqrt(sum(nt.^2,2)) ;
%------------------------------------------ tria barycentres
    [pt] =(pv(tv(:,1),:) + ...
           pv(tv(:,2),:) + ...
           pv(tv(:,3),:))/+3.;
%-------------------------------- tria-to-cell contributions    
    for ci = 1 : size(cp,1)
        for ti = tp(ci,1) : tp(ci,2)
            ac(ci,1) = ac(ci,1) + at(ti,1) ;
            pc(ci,1) = pc(ci,1) + pt(ti,1) * at(ti,1) ;
            pc(ci,2) = pc(ci,2) + pt(ti,2) * at(ti,1) ;
            pc(ci,3) = pc(ci,3) + pt(ti,3) * at(ti,1) ;
        end
        pc(ci,:) = pc(ci,:) / ac(ci,1) ;
    end

end
