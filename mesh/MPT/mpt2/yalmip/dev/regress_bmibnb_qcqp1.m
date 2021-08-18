function fail = regress_bmibnb(ops)

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

sol = solvesdp(F,objective,ops);
fail=getfail(sol.problem,double(objective),-4,checkset(F));
   