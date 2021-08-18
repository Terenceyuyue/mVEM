
clc; clear;close all;

%% Convex polytopal boundary
Pb = Polybnd_Rectangle(0,1);
Pb = Polybnd_Cube(0,1,0,1,0,1);

[~,Pb] = bucky();
Pb = uniquetol(Pb,1e-5,'ByRows',true);

d = size(Pb,2);  % dimension

%% Mesh
maxIter = 10; NT = 20;
% 2-D
if d == 3
    [node,elem] = PolyMesherBd(Pb,maxIter,NT);
    %Nx = 5; Ny = 5;
    %[node,elem] = PolyMesherBd(Pb,maxIter,Nx,Ny);
end

% 3-D
if d == 3
    [node,elem] = PolyMesherBd(Pb,maxIter,NT);
    %Nx = 5; Ny = 5; Nz = 5;
    %[node,elem] = PolyMesherBd(Pb,maxIter,Nx,Ny,Nz);
end

%% Plot
% 2-D
if d==2
    figure, showmesh(node,elem); axis equal
    %findnode(node); %findelem(node,elem);
    return;
end

% 3-D
option.FaceAlpha = 1;
figure, 
showmesh(node,elem,option); axis equal

hold on
iel = 1;
elemf = elem{iel};
option.facecolor = 'yellow';
showmesh(node,elemf,option), axis equal
range = unique(horzcat(elemf{:}));
% findnode(node,range)

figure,
iel = 1;
elemf = elem{iel};
showmesh(node,elemf), axis equal
findnode(node,elemf)
