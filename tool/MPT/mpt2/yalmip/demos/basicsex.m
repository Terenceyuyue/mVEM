clc
echo on
%*********************************************************
%
% Basic manipulation of variables and constraints
%
%*********************************************************
%
% Let us first introduce the central part in YALMIP,
% the sdpvar object
%
% The sdpvar object is typically the decision variable
% in the optimization problem
pause % Strike any key to continue. 

% To define a symmetric 4x4 matrix, we write
P = sdpvar(4,4,'symmetric');
pause % Strike any key to continue. 

% Square matrices are assumed to be symmetric by default,
% so the same result is obtained with
P = sdpvar(4,4);
pause % Strike any key to continue. 

% If we want a fully parameterized matrix, 
% i.e. not neceserally symmetric, we write
P = sdpvar(4,4,'full');
pause % Strike any key to continue. 

% Since non-square matrices are full by default,
% the following two commands are equivalent
P = sdpvar(3,4,'full');
P = sdpvar(3,4);
pause % Strike any key to continue. 

% Of course, a scalar is defined as a 1x1 matrix
x = sdpvar(1,1);
pause % Strike any key to continue. 

% In addition to these two kind of matrices, YALMIP
% also support Toeplitz, Hankel and skew-symmetric matrices
% (all of these must be square)
P = sdpvar(4,4,'toeplitz');
P = sdpvar(4,4,'hankel');
P = sdpvar(4,4,'skew');
pause % Strike any key to continue. 

clc

% Now we know how to define symbolic variables, but how do we work
% with them? Simple, just apply standard Matlab code! All linear operators
% are overloaded.
v = trace(P)+P(1,1)+sum(sum(P))*5;
X = [P diag(P) v*eye(4)];
pause % Strike any key to continue. 

clc

% Finally, a command that can be good in order to understand the basics
% of the sdpvar object, the command see. All symbolic variables can be
% written as a sum of base matrices
%
% P = P0+p1*P1+p2*P2+...+pn*Pn
%
% With see, it is possible to see what these base matrices look like.
% Let us look at the base matrices needed to define a 2x2 symmetric matrix
P = sdpvar(2,2);
see(P)
pause % Strike any key to continue. 
clc
% If we perform any operation on a sdpvar object, this automatically changes
% the base matrices
see(eye(2)+4*P)
pause % Strike any key to continue. 

clc
%*********************************************************
%
% The second most important concept in YALMIP is the set 
% object. A set object is basically a collection of sdpvar 
% objects constrained to have some property, to be part 
% of a set, such as the set of positive semi-definite 
% matrices or the positive orthant.
%
% NOTE. 
% 1) The command SET when operating on sdpvar should not
% be confused with the command SET used to alter properties
% of, e.g, plots, in MATLAB.
% 2) The command SET was previously called LMI in YALMIP, 
% and is based on the LMI object.
%
%*********************************************************

pause % Strike any key to continue. 

% As a first example, constraining a Hermitian matrix to be 
% positive semi-definite is done by
P = sdpvar(3,3);
F = set(P > 0);
pause % Strike any key to continue. 

% Element-wise inequalities are assumed if the argument is not Hermitian
F = set(diag(P) > 0);
pause % Strike any key to continue. 

% Element-wise constraints on a Hermitian matrix is most easily
% done using the following code
F = set(P(:) > 0);
pause % Strike any key to continue. 

% A convenient extension in YALMIP 3 is constraint lists which
% can be used, e.g., to define bounds easily
F = set(0 < diag(P) < 5);
pause % Strike any key to continue. 

% Finally, equality constraints are obtained with ==. 
% Constraining the diagonal elements to be zero is written as
F = set(diag(P) == zeros(3,1));
pause % Strike any key to continue. 

% ...or (consistent with MATLABs notation)
F = set(diag(P) == 0);
pause % Strike any key to continue. 


% Actually, comparasion between a scalar and a matrix is done in a 
% matlab fashion so the following more simple code is valid
F = set(diag(P) == 0);
pause % Strike any key to continue. 

% Alternatively, string notation can be used in all cases
F = set('diag(P) == 0');

pause % Strike any key to continue. 

clc

% Okey, we now know how to define one constraint. What if we have many?
% Simple, just add the set objects (this gives the intersection of two constraints)
F = set('P>0') + set(0 < diag(P) < 5)

pause % Strike any key to continue. 


clc
% We can tag a set with a name
F = set(P>0,'Positive definite')+set(diag(P)==0,'Zero diagonal');
pause % Strike any key to continue. 

% If we display the set, we see the tags
F
pause % Strike any key to continue. 

% The tags can be used for example when we want index the list
% of sets. For example, let us extract the equality constraint
F_equality = F('Zero diagonal')
pause

clc

% To summarize, a constraint can be defined in two ways
P = sdpvar(3,3);

F = set('P>0');

F = set(P > 0);

% ... are most easily built using overloaded +

F = F + set(diag(P)==0)

pause % Strike any key to continue. 
echo off