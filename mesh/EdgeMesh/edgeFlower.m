function [P1,P2] = edgeFlower(N)

a = 5;
r0 = 0.5; r1 = 0.2;
theta = linspace(0,2*pi,N+1)';
r = r0 + r1*sin(a*theta);
x = r.*cos(theta);  
y = r.*sin(theta);  
P1 = [x(1:end-1), y(1:end-1)];
P2 = [x(2:end), y(2:end)];
P2(end,:) = [x(1), y(1)];