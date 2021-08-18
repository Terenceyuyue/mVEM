clc;clear;close all

Domain = @Cube_Domain;
NT = 10;
% Generate intial pointset
P = PolyMesher_init_Pointset3D(Domain,NT);
% Compute the reflection pointset
BdBox = Domain('BdBox');
volume = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))*(BdBox(6)-BdBox(5));
R_P = PolyMesher_Reflect3D(P,Domain,volume);
% Construct the Voronoi diagram
[V,C] = voronoin([P;R_P],{'Qbb','Qz'});
% Generate a particular cell of the Voronoi diagram
iel = 1;
X = V(C{iel},:);
Tri = convhulln(X);
% Show the cell
option.FaceAlpha = 1;
figure,showmesh(X,Tri,option); axis equal

elemTri = cell(NT,1);
for iel = 1:NT
    index = C{iel}; % global indices of element
    X = V(index,:); Tri = convhulln(X);
    elemTri{iel} = index(Tri);
end
totalTri = vertcat(elemTri{:});
option.FaceAlpha = 1;
figure, showmesh(V,totalTri,option); axis equal