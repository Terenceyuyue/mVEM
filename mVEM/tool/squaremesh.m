function [node,elem] = squaremesh(square,h1,h2)
%Squaremesh uniform mesh of a square
%
% square = [a1,b1,a2,b2] for rectangle [a1,b1]*[a2,b2]
%
% Copyright (C) Terence YUE Yu. 

if nargin == 2, h2 = h1; end

% ----------- Generate nodes ---------
a1 = square(1); b1 = square(2); a2 = square(3); b2 = square(4);
[x,y] = ndgrid(a1:h1:b1,a2:h2:b2);
node = [x(:),y(:)];

% -------- Generate elements ---------
nx = size(x,1); ny = size(y,2); % number of columns and rows

% 4 k+nx --- k+1+nx 3  
%    |        |
% 1  k  ---  k+1    2

% indices of k
N = size(node,1); 
k = (1:N-nx)';  cut = nx*(1:ny-1); k(cut) = [];

elem = [k+1 k+1+nx k; k+nx k k+1+nx];