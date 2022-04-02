function centroid = polycentroid3(V,Tri)
% polycentroid3 returns the x,y,z coordinates of centroid
% of surface triangulated polyhedron.
%
% INPUT:
%     vertex: Point Cloud of Shape
%         V(:,1) : x coordinates
%         V(:,2) : y coordinates
%         V(:,3) : z coordinates
%     Tri: connectivity list of face triangulation

% AUTHOR:
%   Isfandiyar RASHIDZADE
%   Email : irashidzade@gmail.com
%   Web Site: isfzade.info
%   Year: 2016
%
%  Modified by Terence Yu.

if iscell(Tri) % cell for face
    elemf = Tri; 
    Tri = faceTriangulation(elemf);
end

vector1 = V(Tri(:,2), :) - V(Tri(:,1), :);
vector2 = V(Tri(:,3), :) - V(Tri(:,1), :);

triangAreasTmp = 0.5*cross(vector1,vector2);
triangAreas = (triangAreasTmp(:,1).^2+triangAreasTmp(:,2).^2 ...
    +triangAreasTmp(:,3).^2).^(1/2); %area of each triangle

totArea = sum(triangAreas); %total area

point1 = V(Tri(:, 1), :);
point2 = V(Tri(:, 2), :);
point3 = V(Tri(:, 3), :);

centroidTri = 1/3*(point1 + point2 + point3); % cent. of each triangle

% mg(:,1) = triangAreas .*  centroidTri(:,1);
% mg(:,2) = triangAreas .*  centroidTri(:,2);
% mg(:,3) = triangAreas .*  centroidTri(:,3);
% mg2 = [mg(:,1), mg(:,2), mg(:,3)];
mg2 = repmat(triangAreas, 1,3).* centroidTri;

centroid = sum(mg2)./totArea;