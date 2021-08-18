clc
echo on
%*********************************************************
%
% Global bilinear programming
%
%*********************************************************
%
% YALMIP comes with a simple branch-and-bound code
% for solving problems with bilinear/quadratic constraints.
%
% The code is under development, and is currently only 
% applicable to small problems. 
%
% For the examples below to work, you must have an efficient
% LP solver (cplex,clp,nag,...) insrtalled and a solver for 
% general polynomial problems (currently only fmincon or PENBMI)


pause
yalmip('clear')
clc
%*********************************************************
% The first example is the problem from the moment 
% relaxation demo (concave quadratic constraint).
%*********************************************************
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);

objective = -2*x1+x2-x3;

F = set(x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24>0);
F = F + set(4-(x1+x2+x3)>0);
F = F + set(6-(3*x2+x3)>0);
F = F + set(x1>0);
F = F + set(2-x1>0);
F = F + set(x2>0);
F = F + set(x3>0);
F = F + set(3-x3>0);
pause

% Explicitly specify the global solver (this is currently recommended)
ops = sdpsettings('solver','bmibnb');                 % Global solver
pause

% Solve!
solvesdp(F,objective,ops)
pause
yalmip('clear')
clc
%***********************************************************
% The second example is slightly larger. It is a problem 
% with a concave quadratic objective function, and two 
% concave quadratic constraints.
%***********************************************************
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);

objective = -25*(x1-2)^2-(x2-2)^2-(x3-1)^2-(x4-4)^2-(x5-1)^2-(x6-4)^2;

F = set((x3-3)^2+x4>4) + set((x5-3)^2+x6>4);
F = F + set(x1-3*x2<2) + set(-x1+x2<2) + set(x1+x2>2);
F = F + set(6>x1+x2>2);
F = F + set(1<x3<5) + set(0<x4<6)+set(1<x5<5)+set(0<x6<10)+set(x1>0)+set(x2>0);

pause

% The global solver in YALMIP can currently only deal with linear objectives.
% A simple epi-graph formulation solves this.
t = sdpvar(1,1);
F = F + set(objective<t);
pause

% Solve!
solvesdp(F,t,ops);
pause
clc
%***********************************************************
% Random bound-constrained indefinte quadratic program.
%***********************************************************

n = 5;
x = sdpvar(n,1);
t = sdpvar(1,1);

A = randn(n,n);A = A*A';
B = randn(n,n);B = B*B';

objective = x'*A*x-x'*B*x;

F = set(-1 < x < 1) + set(objective<t);
pause

% Solve!
solvesdp(F,t,ops);
pause
yalmip('clear')
clc
%***********************************************************
% Random boolean quadratic programming.
%
% By adding a nonlinear constraint x(x-1)==0, we can solve
% boolean programming using global nonlinear programming.
%***********************************************************

n = 10;
x = sdpvar(n,1);
t = sdpvar(1,1);

Q = randn(n,n);Q = Q*Q'/norm(Q)^2;
c = randn(n,1);

objective = x'*Q*x+c'*x;

F = set(x.*(x-1)==0) + set(objective < t);
pause

% Note, all complicating variables (variables involved in nonlinear terms)
% must be bounded (either explicitely or implicetely via linear and relaxed 
% nonlinear constraints).
%
% The easible set for the relaxed problem is not bounded
% so YALMIP will not be able to derive variable bounds.
% We add some valid bounds to fix this
F = F + set(0<x<1);
pause

% Solve!
solvesdp(F,t,ops);
double(x'*Q*x+c'*x)
pause

% Let us compare with an integer solver (use we use the internal MIQP solver)
solvesdp(set(binary(x)),x'*Q*x+c'*x,sdpsettings(ops,'solver','bnb'));
double(x'*Q*x+c'*x)
pause

yalmip('clear')
clc
%***********************************************************
% This example is a somewhat larger indefinite quadratic
% programming problem, taken from the GAMS problem library.
%***********************************************************
pause

% The problem has 11 variables
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);
x7 = sdpvar(1,1);
x8 = sdpvar(1,1);
x9 = sdpvar(1,1);
x10 = sdpvar(1,1);
x11 = sdpvar(1,1);
t = sdpvar(1,1);
x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11];

pause

% ...and 22 linear constraints
e1 =  - x1 - 2*x2 - 3*x3 - 4*x4 - 5*x5 - 6*x6 - 7*x7 - 8*x8 - 9*x9 - 10*x10      - 11*x11 	- 0;
e2 =   - 2*x1 - 3*x2 - 4*x3 - 5*x4 - 6*x5 - 7*x6 - 8*x7 - 9*x8 - 10*x9 - 11*x10      - x11 	- 0;
e3 =   - 3*x1 - 4*x2 - 5*x3 - 6*x4 - 7*x5 - 8*x6 - 9*x7 - 10*x8 - 11*x9 - x10      - 2*x11 	- 0;
e4 =   - 4*x1 - 5*x2 - 6*x3 - 7*x4 - 8*x5 - 9*x6 - 10*x7 - 11*x8 - x9 - 2*x10      - 3*x11 	- 0;
e5 =   - 5*x1 - 6*x2 - 7*x3 - 8*x4 - 9*x5 - 10*x6 - 11*x7 - x8 - 2*x9 - 3*x10      - 4*x11 	- 0;
e6 =   - 6*x1 - 7*x2 - 8*x3 - 9*x4 - 10*x5 - 11*x6 - x7 - 2*x8 - 3*x9 - 4*x10      - 5*x11 	- 0;
e7 =   - 7*x1 - 8*x2 - 9*x3 - 10*x4 - 11*x5 - x6 - 2*x7 - 3*x8 - 4*x9 - 5*x10      - 6*x11 	- 0;
e8 =   - 8*x1 - 9*x2 - 10*x3 - 11*x4 - x5 - 2*x6 - 3*x7 - 4*x8 - 5*x9 - 6*x10      - 7*x11 	- 0;
e9 =   - 9*x1 - 10*x2 - 11*x3 - x4 - 2*x5 - 3*x6 - 4*x7 - 5*x8 - 6*x9 - 7*x10      - 8*x11 	- 0;
e10 =   - 10*x1 - 11*x2 - x3 - 2*x4 - 3*x5 - 4*x6 - 5*x7 - 6*x8 - 7*x9 - 8*x10       - 9*x11 	- 0;
e11 =   - 11*x1 - x2 - 2*x3 - 3*x4 - 4*x5 - 5*x6 - 6*x7 - 7*x8 - 8*x9 - 9*x10       - 10*x11 	- 0;
e12 =     x1 + 2*x2 + 3*x3 + 4*x4 + 5*x5 + 6*x6 + 7*x7 + 8*x8 + 9*x9 + 10*x10       + 11*x11 	- 66;
e13 =     2*x1 + 3*x2 + 4*x3 + 5*x4 + 6*x5 + 7*x6 + 8*x7 + 9*x8 + 10*x9 + 11*x10       + x11 	- 66;
e14 =     3*x1 + 4*x2 + 5*x3 + 6*x4 + 7*x5 + 8*x6 + 9*x7 + 10*x8 + 11*x9 + x10       + 2*x11 	- 66;
e15 =     4*x1 + 5*x2 + 6*x3 + 7*x4 + 8*x5 + 9*x6 + 10*x7 + 11*x8 + x9 + 2*x10       + 3*x11 	- 66;
e16 =     5*x1 + 6*x2 + 7*x3 + 8*x4 + 9*x5 + 10*x6 + 11*x7 + x8 + 2*x9 + 3*x10       + 4*x11 	- 66;
e17 =     6*x1 + 7*x2 + 8*x3 + 9*x4 + 10*x5 + 11*x6 + x7 + 2*x8 + 3*x9 + 4*x10       + 5*x11 	- 66;
e18 =     7*x1 + 8*x2 + 9*x3 + 10*x4 + 11*x5 + x6 + 2*x7 + 3*x8 + 4*x9 + 5*x10       + 6*x11 	- 66;
e19 =     8*x1 + 9*x2 + 10*x3 + 11*x4 + x5 + 2*x6 + 3*x7 + 4*x8 + 5*x9 + 6*x10       + 7*x11 	- 66;
e20 =     9*x1 + 10*x2 + 11*x3 + x4 + 2*x5 + 3*x6 + 4*x7 + 5*x8 + 6*x9 + 7*x10       + 8*x11 	- 66;
e21 =     10*x1 + 11*x2 + x3 + 2*x4 + 3*x5 + 4*x6 + 5*x7 + 6*x8 + 7*x9 + 8*x10       + 9*x11 	- 66;
e22 =     11*x1 + x2 + 2*x3 + 3*x4 + 4*x5 + 5*x6 + 6*x7 + 7*x8 + 8*x9 + 9*x10       + 10*x11 	- 66;

pause


% Define the objective function and the whole program
obj =  (0.5*x1*x2 - x1*x1 + 0.5*x2*x1 - x2*x2 + 0.5*x2*x3 + 0.5*x3*x2 - x3*x3       + 0.5*x3*x4 + 0.5*x4*x3 - x4*x4 + 0.5*x4*x5 + 0.5*x5*x4 - x5*x5 + 0.5*x5      *x6 + 0.5*x6*x5 - x6*x6 + 0.5*x6*x7 + 0.5*x7*x6 - x7*x7 + 0.5*x7*x8 + 0.5      *x8*x7 - x8*x8 + 0.5*x8*x9 + 0.5*x9*x8 - x9*x9 + 0.5*x9*x10 + 0.5*x10*x9       - x10*x10 + 0.5*x10*x11 + 0.5*x11*x10 - x11*x11);
e = -[e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20;e21;e22];
F = set(10>x>0);
F = F + set(e>0);
F = F + set(obj==t);

% We will nowe try to solve this using the same strategy as before
solvesdp(F,t,ops)

pause
clc

% As an alternative, we might use additional cuts using the CUT functionality.
% These additional constraints are only used in domain reduction and
% solution of lower bound programs, but not in the local solver.

% We square all linear constraints and add these as cuts
%
% Note that the triu operator leaves returns a matrix with
% zeros below the diagonal. These are however neglected by YALMIP
% in the call to SET.
%
f = sdpvar(F(find(is(F,'linear'))));
F = F + cut(triu(f*f') > 0);
pause

% ...and solve (To speed up things, we turn off LP based domain reduction)
%
% NOTE : The first ~15 iterations are rather slow, but things starts to 
%        speed up when YALMIP starts pruning redundnant linear constraints.
%
solvesdp(F,t,sdpsettings(ops,'bmibnb.lpreduce',0))
pause
clc

% The global solver in YALMIP can only solve bilinear programs, 
% but a pre-processor is capable of transforming higher order problems 
% to bilinear programs. As an example, the variable x^3y^2 will be replaced 
% with the the variable w and the constraints w == uv, u == zx, z == x^2, v == y^2. 
% The is done automatically if the global solver, or PENBMI, is called with 
% a higher order polynomial problem. Note that this conversion is rather 
% inefficient, and only very small problems can be addressed using this simple approach.
% Additionally, this module has not been tested in any serious sense.
pause

% Let us solve a small trivial polynomial program
sdpvar x y
F = set(x^3+y^5 < 5) + set(y > 0);
F = F + set(-100 < [x y] < 100) % Always bounded domain
pause

solvesdp(F,-x,ops)
pause

%***********************************************************
% The main motivation for implementing a bilinear solver in
% YALMIP is matrix-valued problems, i.e. BMI problems.
%
% Of-course, BMI problems are even harder than the nonconvex 
% quadratic problems solved above. However, in some cases, 
% it is actually possible to solve BMI problems globally.
%
% Don't expect too much though...
%***********************************************************

pause
clc

%***********************************************************
% The following is a classical problem
%***********************************************************
A0 = [-10 -0.5 -2;-0.5 4.5 0;-2 0 0];
A1 = [9 0.5 0;0.5 0 -3 ; 0 -3 -1];
A2 = [-1.8 -0.1 -0.4;-0.1 1.2 -1;-0.4 -1 0];
K12 = [0 0 2;0 -5.5 3;2 3 0];

x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);

F = set(A0+x*A1+y*A2+x*y*K12-t*eye(3)<0);
F = F + set(-0.5 < x < 2) + set(-3 < y < 7);
pause

% Specify solvers (this requires PENBMI!)
ops = sdpsettings('solver','bmibnb');                 % Global solver
ops = sdpsettings(ops,'bmibnb.uppersolver','penbmi'); % Local solver

pause

% Solve!
solvesdp(F,t,ops);
pause
clc
%***********************************************************
% The next example shows how to solve a standard LQ control
% problem with constraints on the feedback gain
%***********************************************************

A = [-1 2;-3 -4];B = [1;1];

P = sdpvar(2,2);
K = sdpvar(1,2);

F = set(P>0)+set((A+B*K)'*P+P*(A+B*K) < -eye(2)-K'*K);
F = F + set(K<0.1)+set(K>-0.1);
F = F + set(1000> P(:)>-1000);
solvesdp(F,trace(P),ops);
pause
clc
%***********************************************************
% Decay-rate estimation as a global BMI example ...
%
% Note, this is a quasi-convex problem that PENBMI actually 
% solves to global optimality, so the use of a global solver
% is pretty stupid. It is only included here to show some 
% properties of the global solver.
%
% may run for up to 100 iterations depending on the solvers
%***********************************************************

A = [-1 2;-3 -4];
t = sdpvar(1,1);
P = sdpvar(2,2);
F = set(P>eye(2))+set(A'*P+P*A < -2*t*P);
F = F + set(-1e4 < P(:) < 1e4) + set(100 > t > 0) ;
pause
solvesdp(F,-t,ops);
pause

% Now, that sucked, right! 100 iterations and still 10% gap...
%
% The way we model a problem can make a huge difference.
% The decay-rate problem can be substantially improved by 
% using a smarter model. The original problem is homogenous 
% w.r.t P, and a better way to make the feasible set bounded
% is to constrain the trace of P
A = [-1 2;-3 -4];
t = sdpvar(1,1);
P = sdpvar(2,2);
F = set(P>0)+set(A'*P+P*A < -2*t*P);
F = F + set(trace(P)==1) + set(100 > t > 0) ;
pause
solvesdp(F,-t,ops);
pause

% Nothing prevents us from adding an additional 
% constraint derived from trace(P)=1
F = F + set(trace(A'*P+P*A) < -2*t);
pause

solvesdp(F,-t,ops);
pause

% When we add redundant constraints, we improve the relaxations
% and obtain better lower bounds. For the upper bounds however,
% these cuts only add additional computational burden.
% To overcome this, the command CUT can be used.
% Constraints defined with CUT are not used when the upper 
% bound problems are solved, but are only used in relaxations
F = set(P>0)+set(A'*P+P*A < -2*t*P);
F = F + set(trace(P)==1) + set(100 > t > 0) ;
F = F + cut(trace(A'*P+P*A) < -2*t);
pause

solvesdp(F,-t,ops);
pause

% Finally, it should be mentioned that the branch & bound
% algorithm can run without any installed local BMI solver. 
% This version of the algorithms is obtained by specifying 
% the upper bound solver as 'none'
ops = sdpsettings(ops,'bmibnb.uppersolver','none');
F = set(P>0)+set(A'*P+P*A < -2*t*P);
F = F + set(trace(P)==1) + set(100 > t > 0) ;
F = F + cut(trace(A'*P+P*A) < -2*t);
pause

solvesdp(F,-t,ops);
pause

echo off



break

% ********************************
% Sahinidis (-6.666 in 1)
% ********************************
yalmip('clear');
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
F = set(0 < x1 < 6)+set(0 < x2 < 4) + set(x1*x2 < 4) 
solvesdp(F,-x1-x2,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bmibnb.maxiter',100,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','cdd','cdd.method','dual-simplex'))

% ****************************
% Decay-rate (-2.5)
% ****************************
yalmip('clear');
A = [-1 2;-3 -4];
t = sdpvar(1,1);
P = sdpvar(2,2);
F = set(P>eye(2))+set(A'*P+P*A < -2*t*P);
F = F + set(-1e4 < P(:) < 1e4) + set(100 > t > -100) ;
ops = sdpsettings('bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bmibnb.maxiter',100,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','cdd','cdd.method','dual-simplex');
solvesdp(F,-t,ops);

% ********************************
% Classical example from Safonov etc
% ********************************
yalmip('clear')
x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
A0 = [-10 -0.5 -2;-0.5 4.5 0;-2 0 0];
A1 = [9 0.5 0;0.5 0 -3 ; 0 -3 -1];
A2 = [-1.8 -0.1 -0.4;-0.1 1.2 -1;-0.4 -1 0];
K12 = [0 0 2;0 -5.5 3;2 3 0];
F = lmi(x>-0.5)+lmi(x<2) + lmi(y>-3)+lmi(y<7) + lmi(t>-1000)+lmi(t<1000);%set(x+y<2.47)+set(x+y>2.47);
F = F + lmi(A0+x*A1+y*A2+x*y*K12-t*eye(3)<0);
solvesdp(F,t,sdpsettings('bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.09e-1,'bmibnb.lpsolver','glpk'))

% ********************************
% Lassere 1
% optimal : -4, takes 13
% ********************************
clear all
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x = [x1;x2;x3];
B = [0 0 1;0 -1 0;-2 1 -1];
b = [3;0;-4];
v = [0;-1;-6];
r = [1.5;-0.5;-5];
p = -2*x1+x2-x3;
F = set(x1+x2+x3 < 4) + set(x1<2)+set(x3<3) + set(3*x2+x3 < 6);
F = F + set(x>0) + set(x'*B'*B*x-2*r'*B*x+r'*r-0.25*(b-v)'*(b-v) >0);
solvesdp(F,p,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.09e-1,'bmibnb.lpsolver','glpk'))

% ********************************
% Lassere 2 
% optimal : -310 in 3
% ********************************
clear all
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);
t = sdpvar(1,1);
p = -25*(x1-2)^2-(x2-2)^2-(x3-1)^2-(x4-4)^2-(x5-1)^2-(x6-4)^2;
F = set((x3-3)^2+x4>4)+set((x5-3)^2+x6>4);
F = F + set(x1-3*x2<2)+set(-x1+x2<2);
F = F + set(x1-3*x2 <2)+set(x1+x2>2);
F = F + set(6>x1+x2>2);
F = F + set(1<x3<5) + set(0<x4<6)+set(1<x5<5)+set(0<x6<10)+set(x1>0)+set(x2>0);

solvesdp(F,p,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.09e-1,'bmibnb.lpsolver','glpk'))
F = F + set(p<t);
solvesdp(F,t,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.09e-1,'bmibnb.lpsolver','glpk'))
% ********************************
% Lassere 2.6
% optimal : -39 in 6
% ********************************
clear all
Q = 100*eye(10);
c = [48, 42, 48, 45, 44, 41, 47, 42, 45, 46]';
b = [-4, 22,-6,-23,-12]';
A =[-2 -6 -1 0 -3 -3 -2 -6 -2 -2;
6 -5 8 -3 0 1 3 8 9 -3;
-5 6 5 3 8 -8 9 2 0 -9;
9 5 0 -9 1 -8 3 -9 -9 -3;
-8 7 -4 -5 -9 1 -7 -1 3 -2];
x = sdpvar(10,1);
t = sdpvar(1,1);
p = c'*x-0.5*x'*Q*x;
F = set(0<x<1)+set(A*x<b);
solvesdp(F,p,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))
F = F + set(p<t);
solvesdp(F,t,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))
% ********************************
% Lassere
% -0.37? >100
% ********************************
clear all
x = sdpvar(10,1);
t = sdpvar(1,1);
x(1) = 1-sum(x(2:end));
obj = 0;
for i = 1:9
    obj = obj+x(i)*x(i+1);
end
for i = 1:8
    obj = obj+x(i)*x(i+2);
end
obj = obj + x(1)*x(7)+ x(1)*x(9)+ x(2)*x(10)+ x(4)*x(7);
F = set(1>x>0);
F = F + set(-obj<t);
solvesdp(F,-obj,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))
F = F + set(-obj<t);
solvesdp(F,t,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))

% ********************************
% Is J co-positive ?
% Optimal : 0
% ********************************
clear all
t = sdpvar(1,1);
x = sdpvar(5,1);
J = [1 -1  1  1 -1;
    -1  1 -1  1  1;
    1 -1  1 -1  1;
    1  1 -1  1 -1;
    -1  1  1 -1  1];
F = set(x>0) + set(sum(x)<1);
F = F+set(-10<t<10);
solvesdp(F,x'*J*x,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))
F = F + set(x'*J*x<t);
solvesdp(F,t,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))

% ********************************
% BARON example problems
% -17 in 9
% ********************************   
clear all
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
t = sdpvar(1,1);
p = 42*x1-50*x1^2+44*x2-50*x2^2+45*x3-50*x3^2+47*x4-50*x4^2+47.5*x5-50*x5^2;
F = set([20 12 11 7 4]*[x1;x2;x3;x4;x5] < 40) + set(0<[x1;x2;x3;x4;x5]<1);
solvesdp(F,p,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','nag','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))
F = F+set(p<t);
solvesdp(F,t,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','nag','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))

% ********************************
% BARON example
% -13 in 1
% ********************************   
clear all
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
t = sdpvar(1,1);
p = x1 - x2 - x3 - x1*x3 + x1*x4 + x2*x3 - x2*x4;
F = set(x1+4*x2 <= 8) + set(4*x1+x2 <= 12) + set(3*x1+4*x2 <= 12) + set(2*x3+x4 <= 8) + set(x3+2*x4 <= 8) + set(x3+x4 <= 5)+set([x1;x2;x3;x4]>0);
F = F+set(p<t);
solvesdp(F,p,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.09e-1,'bmibnb.lpsolver','glpk'))
F = F+set(p<t);
solvesdp(F,t,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.09e-1,'bmibnb.lpsolver','glpk'))

% ********************************
% Control design
% 5.46 in 12
% ********************************   
A = [1 2;-3 0];B = [1;1];
[K0,P0] = lqr(A,B,eye(2),1);
P = sdpvar(2,2);setsdpvar(P,2*P0);K0(K0>1)=1;K0(K0<-1)=-1;
K = sdpvar(1,2);setsdpvar(K,-K0);
F = set(K<1)+set(K>-1)+set(P>0)+set((A+B*K)'*P+P*(A+B*K) < -eye(2)-K'*K);
%F = F + set([-(A+B*K)'*P-P*(A+B*K)-eye(2) K';K eye(1)] > 0)
F = F+lmi(diag(P)>0)+lmi(P(:)>-151) + lmi(P(:)<150) + lmi(P>P0)+lmi(K>-100) + lmi(K<100);
solvesdp(F,trace(P),sdpsettings('usex0',1,'bmibnb.uppersolver','penbmi','bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk','bmibnb.lpreduce',1))

% ********************************
% GAMS st_qpc_m1
% optimal -473.77 in 3
% ********************************
clear all
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
t = sdpvar(1,1);

F = set([x1;x2;x3;x4;x5]>0);
F = F + set(x1 + x2 + 2*x3 + x4 + x5 > 10);
F = F + set(2*x1 + 3*x2 + x5 > 8);
F = F + set(x2 + 4*x3 - x4 + 2*x5 > 12);
F = F + set(8*x1 - x2 - x3 + 6*x4 > 20);
F = F + set(- 2*x1 - x2 - 3*x3 - x4 - x5 > -30);
obj = -(- (10*x1 - 0.34*x1*x1 - 0.28*x1*x2 + 10*x2 - 0.22*x1*x3 + 10*x3 - 0.24*x1*x4 + 10*x4 - 0.51*x1*x5 + 10*x5 - 0.28*x2*x1 - 0.34*x2*x2 - 0.23*x2*x3 -0.24*x2*x4 - 0.45*x2*x5 - 0.22*x3*x1 - 0.23*x3*x2 - 0.35*x3*x3 - 0.22*x3*x4 - 0.34*x3*x5 - 0.24*x4*x1 - 0.24*x4*x2 - 0.22*x4*x3 - 0.2*x4*x4 - 0.38*x4*x5 - 0.51*x5*x1 - 0.45*x5*x2 - 0.34*x5*x3 - 0.38*x5*x4 - 0.99*x5*x5));

solvesdp(F,obj,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))
F = F + set(obj<t);
solvesdp(F,t,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))

% ********************************
% GAMS st_qpk1
% optimal -3 in 7
% ********************************
clear all
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
t = sdpvar(1,1);

F = set([x1;x2]>0);
F = F + set(- x1 + x2 < 1);
F = F + set(x1 - x2 < 1);
F = F + set(- x1 + 2*x2 < 3);
F = F + set(2*x1 - x2 < 3);
obj = (2*x1 - 2*x1*x1 + 2*x1*x2 + 3*x2 - 2*x2*x2);
F = F + set(obj<t);
solvesdp(F,t,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))


% ********************************
% GAMS st_robot
% optimal -3.52 in 4
% ********************************
yalmip('clear')
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);
x7 = sdpvar(1,1);
x8 = sdpvar(1,1);
x = [x1;x2;x3;x4;x5;x6;x7;x8];
F = set(-1<x<1);
F = F + set( 0.004731*x1*x3 - 0.1238*x1 - 0.3578*x2*x3 - 0.001637*x2 - 0.9338*x4 + x7 == 0.3571);
F = F + set(0.2238*x1*x3 + 0.2638*x1 + 0.7623*x2*x3 - 0.07745*x2 - 0.6734*x4 - x7 == 0.6022);
F = F + set( x6*x8 + 0.3578*x1 + 0.004731*x2 == 0);
F = F + set( - 0.7623*x1 + 0.2238*x2 == -0.3461);
F = F + set(x1^2 + x2^2 == 1);
F = F + set(x3^2 + x4^2 == 1);
F = F + set(x5^2 + x6^2 == 1);
F = F + set(x7^2 + x8^2 == 1);

solvesdp(F,sum(x),sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','pensdp','verbose',2,'solver','bmibnb','bnb.maxiter',30,'penbmi.P0',0.9e-1,'bmibnb.lpsolver','glpk'))


% ********************************
% GAMS st_robot
% optimal -382 in 5
% ********************************
yalmip('clear')
x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);
x7 = sdpvar(1,1);
x8 = sdpvar(1,1);
x9 = sdpvar(1,1);
x10 = sdpvar(1,1);
t = sdpvar(1,1);
x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10];

e1=    20*x1 + 20*x2 + 60*x3 + 60*x4 + 60*x5 + 60*x6 + 5*x7 + 45*x8 + 55*x9   + 65*x10 - 600.1;

e2=    5*x1 + 7*x2 + 3*x3 + 8*x4 + 13*x5 + 13*x6 + 2*x7 + 14*x8 + 14*x9      + 14*x10 - 310.5;

e3=    100*x1 + 130*x2 + 50*x3 + 70*x4 + 70*x5 + 70*x6 + 20*x7 + 80*x8 + 80*x9      + 80*x10 - 1800;

e4=    200*x1 + 280*x2 + 100*x3 + 200*x4 + 250*x5 + 280*x6 + 100*x7 + 180*x8      + 200*x9 + 220*x10 - 3850;

e5=    2*x1 + 2*x2 + 4*x3 + 4*x4 + 4*x5 + 4*x6 + 2*x7 + 6*x8 + 6*x9 + 6*x10      - 18.6;

e6=    4*x1 + 8*x2 + 2*x3 + 6*x4 + 10*x5 + 10*x6 + 5*x7 + 10*x8 + 10*x9      + 10*x10 - 198.7;

e7=    60*x1 + 110*x2 + 20*x3 + 40*x4 + 60*x5 + 70*x6 + 10*x7 + 40*x8 + 50*x9      + 50*x10 - 882;

e8=    150*x1 + 210*x2 + 40*x3 + 70*x4 + 90*x5 + 105*x6 + 60*x7 + 100*x8      + 140*x9 + 180*x10 - 4200;

e9=    80*x1 + 100*x2 + 6*x3 + 16*x4 + 20*x5 + 22*x6 + 20*x8 + 30*x9 + 30*x10      - 40.25;

e10=    40*x1 + 40*x2 + 12*x3 + 20*x4 + 24*x5 + 28*x6 + 40*x9 + 50*x10       - 327;

obj =   (10*x1 - 6.8*x1*x1 - 4.6*x1*x2 + 10*x2 - 7.9*x1*x3 + 10*x3 - 5.1*x1*x4       + 10*x4 - 6.9*x1*x5 + 10*x5 - 6.8*x1*x6 + 10*x6 - 4.6*x1*x7 + 10*x7 -      7.9*x1*x8 + 10*x8 - 5.1*x1*x9 + 10*x9 - 6.9*x1*x10 + 10*x10 - 4.6*x2*x1       - 5.5*x2*x2 - 5.8*x2*x3 - 4.5*x2*x4 - 6*x2*x5 - 4.6*x2*x6 - 5.5*x2*x7 -      5.8*x2*x8 - 4.5*x2*x9 - 6*x2*x10 - 7.9*x3*x1 - 5.8*x3*x2 - 13.3*x3*x3 -      6.7*x3*x4 - 8.9*x3*x5 - 7.9*x3*x6 - 5.8*x3*x7 - 13.3*x3*x8 - 6.7*x3*x9 -      8.9*x3*x10 - 5.1*x4*x1 - 4.5*x4*x2 - 6.7*x4*x3 - 6.9*x4*x4 - 5.8*x4*x5 -      5.1*x4*x6 - 4.5*x4*x7 - 6.7*x4*x8 - 6.9*x4*x9 - 5.8*x4*x10 - 6.9*x5*x1 -      6*x5*x2 - 8.9*x5*x3 - 5.8*x5*x4 - 11.9*x5*x5 - 6.9*x5*x6 - 6*x5*x7 - 8.9*      x5*x8 - 5.8*x5*x9 - 11.9*x5*x10 - 6.8*x6*x1 - 4.6*x6*x2 - 7.9*x6*x3 - 5.1      *x6*x4 - 6.9*x6*x5 - 6.8*x6*x6 - 4.6*x6*x7 - 7.9*x6*x8 - 5.1*x6*x9 - 6.9*      x6*x10 - 4.6*x7*x1 - 5.5*x7*x2 - 5.8*x7*x3 - 4.5*x7*x4 - 6*x7*x5 - 4.6*x7      *x6 - 5.5*x7*x7 - 5.8*x7*x8 - 4.5*x7*x9 - 6*x7*x10 - 7.9*x8*x1 - 5.8*x8*      x2 - 13.3*x8*x3 - 6.7*x8*x4 - 8.9*x8*x5 - 7.9*x8*x6 - 5.8*x8*x7 - 13.3*x8      *x8 - 6.7*x8*x9 - 8.9*x8*x10 - 5.1*x9*x1 - 4.5*x9*x2 - 6.7*x9*x3 - 6.9*x9      *x4 - 5.8*x9*x5 - 5.1*x9*x6 - 4.5*x9*x7 - 6.7*x9*x8 - 6.9*x9*x9 - 5.8*x9*      x10 - 6.9*x10*x1 - 6*x10*x2 - 8.9*x10*x3 - 5.8*x10*x4 - 11.9*x10*x5 - 6.9      *x10*x6 - 6*x10*x7 - 8.9*x10*x8 - 5.8*x10*x9 - 11.9*x10*x10);

F = set(0<x<10000);
F = F + set([e1;e2;e3;e4;e5;e6;e7;e8;e9;e10]<0);
%F = F + set(obj==t);
solvesdp(F,obj,sdpsettings('bmibnb.uppersolver','penbmi','bmibnb.lowersolver','glpk','solver','bmibnb','bmibnb.maxiter',15,'penbmi.P0',0.01,'bmibnb.lpsolver','glpk'))
