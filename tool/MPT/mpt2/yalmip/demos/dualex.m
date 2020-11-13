yalmip('clear'); 
clc
echo on
%*********************************************************
%
% Dual variables
%
%*********************************************************
%
% Dual variables can easily be extracted from the
% optimal solution using YALMIP
% 
pause

% Let us solve a Lyapuov stability problem
A = [-1 2;0 -2];
P = sdpvar(2,2); 
F = set(P>eye(2),'Normalize P') + set(A'*P+P*A<0,'Lyapunov stability');
pause % Strike any key to continue. 

% Find minimum trace solution
solution = solvesdp(F,trace(P));
pause % Strike any key to continue. 
clc

% The dual solutions are readily obtained using simple indexing
dual(F(1))
dual(F(2))
pause
clc

% ... or with the tags
dual(F('Normalize P'))
dual(F('Lyapunov stability'))
pause
clc

% Complementary slackness ...
trace(double(F(1))*dual(F(1)))
trace(double(F(2))*dual(F(2)))
pause
echo off