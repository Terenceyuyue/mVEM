function fail = regress_AtimesX(ops)


A = [-1 2;-3 -4];
t = sdpvar(1,1);
P = sdpvar(2,2);
F = set(P>eye(2))+set(A'*P+P*A < -2*t*P);
F = F + set(-1e4 < P(:) < 1e4) + set(100 > t > -100) ;

ops.bmibnb.lowersolver = 'pensdp';

obj = -t;

sol = solvesdp(F,obj,ops);

fail=getfail(sol.problem,double(obj),-2.5,checkset(F));
   