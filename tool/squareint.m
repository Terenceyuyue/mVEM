function I = squareint(f,square)
% squareint computes the integral on a square domain
%
% Copyright (C) Terence Yu.

%% parameter
a1 = square(1); b1 = square(2); 
a2 = square(3); b2 = square(4); 

%% Gauss quadrature
[r,w] = Gaussquad;
xl = a1; xr = b1; h = xr-xl; x = xl+(1+r)./2*h; wx = h/2*w;
xl = a2; xr = b2; h = xr-xl; y = xl+(1+r)./2*h; wy = h/2*w;

%% double integral
[xx,yy] = ndgrid(x,y); wxy = wx*wy';
fwxy = wxy(:).*f(xx(:),yy(:));
I = sum(fwxy(:));