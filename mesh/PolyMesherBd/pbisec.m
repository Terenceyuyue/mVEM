function [A,b] = pbisec(x, P)
% x= [x,y] or [x,y,z]
% This function finds perpendicular bisector between two points in 2D/3D
% Hyongju Park / hyongju@gmail.com
% input: two points in 2D/3D
% output: inequality Ax <= b

n = size(P,1);  d = size(P,2);
A = zeros(n,d); b = zeros(n,1);
for i = 1:n
    mPts = mean([x;P(i,:)],1);
    n_vec = (P(i,:)-x)/norm(P(i,:)-x);
    Ad = n_vec;
    bd = dot(n_vec,mPts);
    if Ad*x(:)>bd, Ad = -Ad; bd = -bd; end
    A(i,:) = Ad;  b(i) = bd;
end

