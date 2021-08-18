DELTA = 1e-4;

y = sdpvar(4,1);
Z = [ -y(2)         1/2+0.5*y(1)  -y(3);
   1/2+0.5*y(1) DELTA         -y(4);
   -y(3)          -y(4)       DELTA ];

F = set(Z > 0)
sol = solvesdp(F,-y(1) - 2*DELTA*y(2),sdpsettings('sdpt3.gaptol',1e-8,'solver','vsdp','savesolveroutput',1,'vsdp.solver','sdpt3'))



DELTA = 1e-4;

y = sdpvar(4,1);
Z = [ -y(2)         1/2+0.5*y(1)  -y(3);
   1/2+0.5*y(1) DELTA         -y(4);
   -y(3)          -y(4)       DELTA ];

F = set(Z > 0)
sol = solvesdp(F,-y(1) - 2*DELTA*y(2),sdpsettings('sdpt3.gaptol',1e-8,'solver','vsdp','savesolveroutput',1))




double(y)
double(y(1) + 2*DELTA*y(2))
sol.solveroutput 


s = midrad(1,0.2);
sdpvar x
sol = solvesdp(set([s + x 1;1 s + x] > 0),x,sdpsettings('solver','vsdp','savesolveroutput',1))


s = midrad(0,0.02);
A = [-1 1+s;0.1 -2];
P = sdpvar(2,2);
F = set(A'*P + P*A <= -eye(2));
solvesdp(F,trace(P),sdpsettings('solver','vsdp','savesolveroutput',1));


s = sdpvar(4,1);
A = [-1+s(1) 1+s(2);0.1+s(3) -2+s(4)];
P = sdpvar(2,2);
F = set(A'*P + P*A <= -eye(2));
F = F + set(-0.02 < s < 0.02) + set(uncertain(s));
solvesdp(F,trace(P),sdpsettings('solver','sedumi'))


s = midrad([0;0;0;0],0.02);
A = [-1+s(1) 1+s(2);0.1+s(3) -2+s(4)];
P = sdpvar(2,2);
F = set(A'*P + P*A <= -eye(2));
solvesdp(F,trace(P),sdpsettings('solver','vsdp','savesolveroutput',1))


solvesdp(F,trace(P),sdpsettings('solver','sedumi'))
