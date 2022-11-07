clc;clear;close all


%% rectangle mesh
square = [0 1 0 1];  Nx = 5; Ny = Nx;
% node
x1 = square(1); x2 = square(2); y1 = square(3); y2 = square(4);
x = linspace(x1,x2,Nx+1);  y = linspace(y1,y2,Ny+1);
[x,y] = ndgrid(x,y);
node = [x(:),y(:)];
% elem
nx = size(x,1); ny = size(y,2); 
N = size(node,1);
k = (1:N-nx)';  cut = nx*(1:ny-1); k(cut) = [];
elem = [k k+1 k+1+nx k+nx];
% edge, elem2edge
totalEdge = sort([elem(:,[1,2]); elem(:,[2,3]); elem(:,[3,4]); elem(:,[4,1])], 2);
[edge, i1, totalJ] = unique(totalEdge ,'rows');
NT = size(elem,1); NE = size(edge,1);
elem2edge = reshape(totalJ ,NT, 4);
% nodeEdge
nodeEdge = (node(edge(:,1),:)+node(edge(:,2),:))/2;

%% add interior points
nodeg = zeros(16*NT,2);  elemg = cell(NT*5,1);
s = 0;  t = 0;
for iel = 1:NT
    % info of currect element
    index = elem(iel,:);  indexEdge = elem2edge(iel,:);  Nv = 4;
    a = node(index(1),1);   b = node(index(2),1); 
    c = node(index(1),2);   d = node(index(4),2);
    % interior points on reference element [-3,3]*[-3,3]
    Pts = [0 -1; 1 -2; 1 -1; 2 -1; 1 0; 2 1; 1 1; 1 2; 0 1; -1 2; -1 1; 
        -2 1; -1 0; -2 -1; -1 -1; -1 -2];
    % nodeg
    Ptg = zeros(16,2);
    Ptg(:,1) = (b-a)/6*(Pts(:,1)+3) + a;
    Ptg(:,2) = (d-c)/6*(Pts(:,2)+3) + c;
    nodeg((1:16)+s,:) = Ptg;
    % elemg
    indexve = zeros(1,2*Nv);  
    indexve(1:2:end) = index; indexve(2:2:end) = indexEdge+N;
    indexg = [indexve, (1:16)+N+NE+s];
    elemg{1+t} = indexg(9:24);
    elemg{2+t} = indexg([1 2 9 24 23 22 21 8]);
    elemg{3+t} = indexg([3 4 13 12 11 10 9 2]);
    elemg{4+t} = indexg([5 6 17 16 15 14 13 4]);
    elemg{5+t} = indexg([7 8 21 20 19 18 17 6]);
    s = s + 16;  t = t + 5;
end

%% Gunelve mesh on the square
node = [node; nodeEdge; nodeg];
elem = elemg;

figure,
showmesh(node,elem)
title('Gunelve mesh on the square');

%% Gunelve mesh on the Cook's membrane
p = [x1 y1; x2 y1; x2 y2; x1 y2];
q = [0 0; 48 44; 48 60; 0 44];
tform = fitgeotrans(p,q,'projective');
node = transformPointsForward(tform, node);

figure,
showmesh(node,elem)
title('Gunelve mesh on the Cook''s membrane');