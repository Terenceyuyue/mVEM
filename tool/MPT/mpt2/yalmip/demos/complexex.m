clc
echo on
%*********************************************************
%
% Complex-valued problems
%
%*********************************************************
%
% YALMIP deals with complex-valued problems as easy as 
% standard problems
pause % Strike any key to continue. 
clc
% As an example, let us study the example on "Toeplitz
% covariance estimation" in the SeDuMi manual.
%
% We are given a complex matrix P, and the goal
% is to find a (possibly complex) Toeplitz matrix
% Z that minimizes the Frobenious norm of Z-P
pause % Strike any key to continue. 

% Data given
i=sqrt(-1);
P = [4 1+2*i 3-i;1-2*i 3.5 0.8+2.3*i;3+i 0.8-2.3*i 4];
pause % Strike any key to continue. 

% Define a complex Toeplitz matrix
Z = sdpvar(3,3,'toeplitz','complex');
% Warning : Be careful and note that this
% is not the same as Z = toeplitz(sdpvar(1,3))+i*toeplitz(sdpvar(1,3))
% or Z = sdpvar(3,3,'toeplitz')+sdpvar(3,3,'toeplitz')
% We will later use the command "see" to understand the difference.
pause % Strike any key to continue. 

% The Toeplitz matrix is not Hermitian, since it has got complex
% elements on the diagonal. These can be removed with the following code
Z = Z-sqrt(-1)*diag(imag(diag(Z)));
pause % Strike any key to continue. 

% An alternative way to accomplish the same thing is
Z = sdpvar(3,3,'toeplitz','complex');
Z = replace(Z,imag(diag(Z)),[0;0;0]);
pause % Strike any key to continue. 

% Now use the fact that the Frobenious norm of Z-P
% equals the 2-norm of Z(:)-P(:)
%
% Hence, we define e = Z(:)-P(:) and minimize ||e||
% Write this as     min t
%                   s.t ||e|| < t
%                           Z > 0
pause % Strike any key to continue. 
t = sdpvar(1,1);
e = Z(:)-P(:);
F = set(Z > 0);
F = F+set('||e||<t');
pause % Strike any key to continue. 

% Solve!
solvesdp(F,t);
pause % Strike any key to continue. 

% Check the result
double(Z)

% and compare with result in SeDuNi manual
Zmanual = toeplitz([4.2827,0.8079+1.7342*sqrt(-1) 2.5574-0.7938*i])
pause % Strike any key to continue. 

clc
% More convenient, use the overloaded norm operator!
solvesdp(set(Z > 0),norm(e));
pause

clc
% The problem we solved above was a complex SOCP. Let us 
% now solve a complex LMI. We solve the same problem, but by
% first applying a Schur complement to obtain 
%                   min t
%                   s.t [t e';e I] > 0
%                                Z > 0
pause % Strike any key to continue. 
F = set(Z);
F = F+set([t e';e eye(9)]);
pause % Strike any key to continue. 

% Solve!
solvesdp(F,t);
pause % Strike any key to continue. 

% Check the result
double(Z)

% and compare with result in SeDuNi manual
Zmanual = toeplitz([4.2827,0.8079+1.7342*sqrt(-1) 2.5574-0.7938*i])
pause % Strike any key to continue. 

clc

% We should also mention that complex linear constraints are interpreted
% just as they are stated, i.e. in the complex domain (there are some
% alternatives in SeDuMi).
%
% Let us solve a simple linear complex problem to illustrate how
% complex constraints are used
pause % Strike any key to continue. 

p = sdpvar(1,1,'full','complex');      % A complex scalar (4 arguments necessary)
s = sdpvar(1,1)+sqrt(-1)*sdpvar(1,1);  % Alternativ definition
F = set('0.9>imag(p)');
F = F+set('0.01>real(p)');
F = F+set('0.1+0.5*sqrt(-1)>p');
F = F+set('s+p==2+4*sqrt(-1)');

pause % Strike any key to continue. 
% If we display the set, we see which constraints that are complex
F
pause % Strike any key to continue. 

% And solve...
solvesdp(F,-real(p)-imag(p));

double(s+p)
double(s)
double(p)
pause % Strike any key to continue. 


clc
% Let us go back to the discussion earlier on how to construct
% a complex Toeplitz matrix.
%
% Since YALMIP supports Toeplitz as a standard matrix in sdpvar, 
% an easily made mistake is to belive that the following code
% would generate a complex Toeplitz
pause % Strike any key to continue. 
Z = sdpvar(3,3,'toeplitz')+i*sdpvar(3,3,'toeplitz');
pause % Strike any key to continue. 

% If we look at the base matrices, we see that this is not the case
% (the complex matrices are not complex Toeplitz)
see(Z,'full')
pause % Strike any key to continue. 

clc

% Another bad idea is the following piece of code
x = sdpvar(1,3);
y = sdpvar(1,3);
Z = toeplitz(x)+i*toeplitz(y);
pause % Strike any key to continue. 
% We see once again that this is not a complex Toeplitz
see(Z,'full')
pause % Strike any key to continue. 

clc

% One way to do it correctly is to write it as we did above
Z = sdpvar(3,3,'toeplitz','complex');
pause % Strike any key to continue. 
% ... complex Toeplitz matrices 
% (notice the diagonal complex terms that we took care of earlier)
see(Z,'full')
pause % Strike any key to continue. 

clc
% As a final remark, remember that the common functions 
% for complex algebra is available for sdpvar objects.
pause % Strike any key to continue. 
Z = sdpvar(3,3,'toeplitz','complex');

real(Z)
imag(Z)
conj(Z)
isreal(Z)
isreal(Z-sqrt(-1)*imag(Z))

pause % Strike any key to continue. 
echo off
