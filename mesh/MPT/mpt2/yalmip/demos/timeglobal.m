t0 = cputime;

yalmip('clear')
clc
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

% Explicitly specify solvers (this is currently recommended)
ops = sdpsettings('solver','bmibnb');                 % Global solver
ops = sdpsettings(ops,'bmibnb.lowersolver','glpk');   % Lower solver
ops = sdpsettings(ops,'bmibnb.uppersolver','penbmi'); % Local solver
ops = sdpsettings(ops,'bmibnb.lpsolver','glpk');      % LP solver
ops = sdpsettings(ops,'verbose',0);      % LP solver
ops = sdpsettings(ops,'penbmi.P0',0.01);     

solvesdp(F,objective,ops)

yalmip('clear')
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


t = sdpvar(1,1);
F = F + set(objective<t);

solvesdp(F,t,ops)


yalmip('clear')
A0 = [-10 -0.5 -2;-0.5 4.5 0;-2 0 0];
A1 = [9 0.5 0;0.5 0 -3 ; 0 -3 -1];
A2 = [-1.8 -0.1 -0.4;-0.1 1.2 -1;-0.4 -1 0];
K12 = [0 0 2;0 -5.5 3;2 3 0];

x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);

F = set(A0+x*A1+y*A2+x*y*K12-t*eye(3)<0);
F = F + set(-0.5 < x < 2) + set(-3 < y < 7);

% Specify solvers
ops = sdpsettings('solver','bmibnb');                 % Global solver
ops = sdpsettings(ops,'bmibnb.lowersolver','pensdp'); % Lower solver
ops = sdpsettings(ops,'bmibnb.uppersolver','penbmi'); % Local solver
ops = sdpsettings(ops,'bmibnb.lpsolver','glpk');      % LP solver
ops = sdpsettings(ops,'penbmi.P0',0.01);              % Typically good

ops = sdpsettings(ops,'verbose',0);              % Typically good

solvesdp(F,t,ops)


A = [-1 2;-3 -4];B = [1;1];

P = sdpvar(2,2);
K = sdpvar(1,2);

F = set(P>0)+set((A+B*K)'*P+P*(A+B*K) < -eye(2)-K'*K);
F = F + set(K<0.1)+set(K>-0.1);
F = F + set(1000> P(:)>-1000);
solvesdp(F,trace(P),ops)

A = [-1 2;-3 -4];
t = sdpvar(1,1);
P = sdpvar(2,2);
F = set(P>eye(2))+set(A'*P+P*A < -2*t*P);
F = F + set(-1e4 < P(:) < 1e4) + set(100 > t > 0) ;
solvesdp(F,-t,ops)

A = [-1 2;-3 -4];
t = sdpvar(1,1);
P = sdpvar(2,2);
F = set(P>0)+set(A'*P+P*A < -2*t*P);
F = F + set(trace(P)==1) + set(100 > t > 0) ;
solvesdp(F,-t,ops)

F = F + set(trace(A'*P+P*A) < -2*t);
solvesdp(F,-t,ops)

F = set(P>0)+set(A'*P+P*A < -2*t*P);
F = F + set(trace(P)==1) + set(100 > t > 0) ;
F = F + cut(trace(A'*P+P*A) < -2*t);

solvesdp(F,-t,ops)

yalmip('clear');
A = [-1 2 ;0.5 0.3];
B = [0;1];

P = sdpvar(2,2);
K = sdpvar(1,2);

F = lmi(P>0);
F = F + lmi([-(A+B*K)'*P-P*(A+B*K)-eye(2) K';K eye(1)] > 0);
F = F+lmi(diag(P)>0)+lmi(P(:)>-150) + lmi(P(:)<150) + lmi(K>-5) + lmi(K<5);
F = F + cut((A+B*K)'*P+P*(A+B*K) < -eye(2)-K'*K);

F = F + cut(kron(P,P)>0);

ops = sdpsettings('bmibnb.lowersolver','pensdp','verbose',0,'solver','bmibnb','pensdp.P0',0.9);
ops = sdpsettings(ops,'bmibnb.vartol',1e-4,'bmibnb.lpreduce',1,'bmibnb.roottight',1,'bmibnb.maxiter',30,'penbmi.P0',0.01,'bmibnb.lpsolver','glpk');
solvesdp(F,trace(P),ops)



cputime-t0
