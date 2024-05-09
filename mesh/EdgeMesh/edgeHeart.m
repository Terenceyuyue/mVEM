function [P1,P2] = edgeHeart(N)

theta = linspace(0,2*pi,N+1)';
r = 1-cos(theta);  

x = r.*cos(theta) + 1;
y = r.*sin(theta);  


P1 = [x(1:end-1), y(1:end-1)];
P2 = [x(2:end), y(2:end)];
P2(end,:) = [x(1), y(1)];
