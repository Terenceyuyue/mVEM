function test_robust_5


% Problem with two separable polytope uncertainties
A1 = randn(8,2);
A2 = randn(8,2);
x0 = randn(2,1);
b1 = A1*x0+1;
b2 = A2*x0+1;

w2 = sdpvar(2,1);
w1 = sdpvar(2,1);
x = sdpvar(1,1);

sol1 = solvesdp([uncertain([w1;w2]), x+sum(w1) + sum(w2) < 0,A1*w1 < b1, A2*w2 < b2],-x,sdpsettings('debug',1))
obj1 = double(-x);
sol2 = solvesdp([A1*w1 < b1,A2*w2 < b2],-sum(w1)-sum(w2))
obj2 = double(-sum(w1)-sum(w2));

mbg_asserttolequal(sol1.problem, 0, 1e-5);
mbg_asserttolequal(obj1+obj2,0, 1e-5);
