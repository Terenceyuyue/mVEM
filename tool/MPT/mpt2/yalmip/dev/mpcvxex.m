
figure

% Create numerics for a discretized double integrator
N = 5;
A = [2 -1;1 0];
B = [1;0];
C = [0.5 0.5];
[H,S] = create_CHS(A,B,C,N);

% Variables
x = sdpvar(2,1);
U = sdpvar(N,1);
Y = sdpvar(N,1);
% Parameterized predictions
Y = H*x+S*U; 

%Control constraints
F = set(-1 < U < 1);

% Terminal constraint
F = F+set(-1 < Y(N) < 1); 

% Exploration space
F = F + set(-5 < x < 5);

ops = sdpsettings('solver','mpcvx','mpcvx.plot',0,'mpcvx.relgaptol',0.01,'mpcvx.solver','quadprog')
[mpsol1,a1,b1,c1,d1] = solvemp(F,Y'*Y+U'*U,ops,x,U(1));


[mpsol2,a2,b2,c2,d2] = solvemp(F,Y'*Y+U'*U,sdpsettings('solver','mpt'),x,U(1));
close all
figure
plot(c1);hold on
plot(c2)
figure
plot(d1);hold on
plot(d2)

% SDP
% Compare with
% http://www.springerlink.com/content/272k25u7852702r8/fulltext.pdf

sdpvar t1 t2 t3 x1 x2 x3

A0 = [1 2 -3; 2 4 -1; -3 -1 3];
A1 = [1 -1 2; -1 1 3; 2 3 2];
A2 = [-1 1 0;1 1 2;0 2 -2];
A3 = [3 -2 4;-2 1 -2;4 -2 -2];
A4 = [-3 1 1;1 -2 -1;1 -1 1];
A5 = [5 4 2;4 1 1;2 1 -1];

F = set(-2 < [t1 t2] < 2);
F = F + set(A0 + A1*t1 + A2*t2 + A3*x1 + A4*x2 + A5*x3);
obj = x1 -2*x2 + x3;

[mpsol1,a1,b1,c1,d1] = solvemp(F,obj,sdpsettings('solver','mpcvx','mpcvx.plot',1,'mpcvx.absgaptol',0.5),[t1 t2]);

plot(c1)

% Compare costs
assign([t1;t2],[0.5;1.3])
double(c1)
solvesdp(F + set([t1;t2]==[0.5;1.3]),obj)
double(obj)



A = randn(10,2);
b = rand(10,1)*2;
W = randn(10,2);
x = sdpvar(2,1);
t = sdpvar(2,1);
solvesdp(set(A*x < b + W*t) + set(-0.1 < t < 0.1),-geomean(b-A*x))

ops = sdpsettings('solver','mpcvx','mpcvx.plot',1,'mpcvx.eps',0.050,'mpcvx.solver','')
[mpsol1,a1,b1,c1,d1] = solvemp(set(-3 < t < 3),-geomean(b+W*t-A*x),ops,t);


