function [cc] = triaball2(pp,tt)
%TRIABALL2 calc. the circumballs associated with a 2-simplex
%triangulation embedded in R^3.
%   [CC] = TRIABALL2(PP,TT) returns the circumscribing balls
%   associated with the triangles in [PP,TT], such that CC = 
%   [XC,YC,ZC,RC*RC].

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 15/12/2014

    cc = zeros(size(tt,1),4);
    rv = zeros(size(tt,1),3);
%------------------------------------------------- tria edge
    ab = pp(tt(:,1),:)-pp(tt(:,2),:);
    ac = pp(tt(:,1),:)-pp(tt(:,3),:);
%------------------------------------------------- tria norm
    nv = cross(ab,ac);
%--------------------------------------- setup coeff. vector
    rv(:,1) = sum(pp(tt(:,1),:).^2,2) - ...
              sum(pp(tt(:,2),:).^2,2) ;
    rv(:,2) = sum(pp(tt(:,1),:).^2,2) - ...
              sum(pp(tt(:,3),:).^2,2) ;
    rv(:,3) = sum(pp(tt(:,1),:).*nv,2);
    rv = rv';
%--------------------------------------- calc. circum-centre    
    aa = zeros(3,3,size(tt,1));
    aa(1,:,:) = +2.*ab';
    aa(2,:,:) = +2.*ac';
    aa(3,:,:) = +1.*nv';

    cc(:,1:3) = blocksolve(aa,rv)';
%--------------------------------------- calc. circum-radius
    cc(:,  4) =(sum((cc(:,1:3)-pp(tt(:,1),:)).^2,2) + ...
                sum((cc(:,1:3)-pp(tt(:,2),:)).^2,2) + ...
                sum((cc(:,1:3)-pp(tt(:,3),:)).^2,2))/+3.; 
    
end

