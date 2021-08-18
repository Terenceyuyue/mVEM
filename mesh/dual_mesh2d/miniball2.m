function cc = miniball2(pp,tt)
%MINIBALL2 compute the min-enclosing balls associated with a 
%2-simplex triangulation embedded in R^3.
%   [CC] = MINIBALL2(PP,TT) returns the smallest enclosing 
%   balls associated with the triangles in [PP,TT], such th-
%   at CC = [XC,YC,ZC,RC*RC]. Such balls never lie outside 
%   the hull of their associated 2-simplexes.

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 29/11/2014

%----------------------------------------- calc. circumballs
    cc = triaball2(pp,tt);
%------------------------ replace with face-balls if smaller
    cc = miniface2(cc,pp(tt(:,1),:), ...
                      pp(tt(:,2),:), ...
                      pp(tt(:,3),:)) ;
    cc = miniface2(cc,pp(tt(:,2),:), ...
                      pp(tt(:,3),:), ...
                      pp(tt(:,1),:)) ;
    cc = miniface2(cc,pp(tt(:,3),:), ...
                      pp(tt(:,1),:), ...
                      pp(tt(:,2),:)) ;
    
end

function cc = miniface2(cc,pi,pj,pk)
%------------------------------------------------ edge balls
    bc = (pi + pj) * +.5 ;
%------------------------------------------------ edge radii
    br = (sum((bc(:,1:3)-pi).^2,2) + ...
          sum((bc(:,1:3)-pj).^2,2))* .5;
%------------------------------------------- enclosing radii
    ll =  sum((bc(:,1:3)-pk).^2,2);
%------------------------------------------- replace if min.
    ki = (br >= ll) & (br <= cc(:,4))  ;
%------------------------------------------- replace is min.
    cc(ki,1:3) = bc(ki,:);
    cc(ki,  4) = br(ki,:);
end
