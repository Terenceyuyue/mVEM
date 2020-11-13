clc
echo on
%*********************************************************
%
% %Second order cone programming
%
%*********************************************************
% 
% The problem is to find the point in the intersection of
% two circles which is closest to the origin
pause % Strike any key to continue. 

% The two circles are centered at a and b with radius 1,
% and the point we are looking for will be denoted x
a = [0;1];
b = [1;1];
x = sdpvar(2,1);
pause % Strike any key to continue. 

% The point lies in both circles
F = set('||x-a||<1')+set('||x-b||<1');
pause

% An alternative way to define the SOCPs is 
F = set(cone(x-a,1))+set(cone(x-b,1));
pause


% The distance to the origin is limited
% by a constant t
t = sdpvar(1,1);
F = F+set('||x||<t');
pause % Strike any key to continue. 


% And we wish to minimize the distance
sol = solvesdp(F,t);
pause % Strike any key to continue. 

% Optimal point
double(x)

% Distance from origin
double(t)

% Which of-course equals ||x||
norm(double(x))

pause % Strike any key to continue. 
echo off