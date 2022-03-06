function [node,elem] = nonConvexMesh(square,Nx,Ny)
%Generate a typical nonconvex polygonal mesh
%
% square = [a1,b1,a2,b2] for rectangle [a1,b1]*[a2,b2]
%
% Copyright (C)  Terence Yu.

if nargin == 2, Ny = Nx; end

% ----------- Generate nodes ---------
a1 = square(1); b1 = square(2); a2 = square(3); b2 = square(4);
x = linspace(a1,b1,2*Nx+1);  h1 = x(2)-x(1);
y = linspace(a2,b2,2*Ny+1);  h2 = y(2)-y(1);
[x,y] = ndgrid(x,y);
node = [x(:),y(:)];
nx = size(x,1); ny = size(y,2); % number of columns and rows

% -------- Generate elements ---------
k = reshape(1:size(node,1), nx, ny);
k = k(1:2:nx-2, 1:2:ny-2); k = k(:);
elem4 = [k k+1 k+2 k+2+nx k+2+2*nx k+1+2*nx k+2*nx k+nx];

% -------- shift midpoints ---------
node(2:2:end, 1) = node(2:2:end, 1) +  h1/2;  % even index
node(2:2:end, 2) = node(2:2:end, 2) +  h2/2;

% -------- update node and elem ---------
% logical vector for deleting redundant nodes
isNode = true(nx,ny);
isNode(2:2:nx, 2:2:ny) = false;   % delete centroid
isNode([1,end], 2:2:end) = false; % delete mid-points on boundary
isNode(2:2:end, [1,end]) = false;
isNode = isNode(:);
% connection number
N = sum(isNode);
idDof = zeros(size(node,1),1);
idDof(isNode) = 1:N;
% update node and elem
node = node(isNode,:);
NT = size(elem4,1); 
elem = cell(NT,1);
for iel = 1:NT
    index = elem4(iel,:);
    isElem = isNode(index);
    elem{iel} = reshape(idDof(index(isElem)),1,[]);
end