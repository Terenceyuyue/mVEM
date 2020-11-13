sdpvar x1 x2 x3 x4 x5 x6 t1 t2

lhs = [-x1+x5;-x1-x5;-x2+x6;-x2-x6;-x3;-x3-x5;-x3;-x3+x5;-x4-x5-x6;-x4+x5;-x4+x5;-x4+x5+x6;x5;-x5;x6;-x6]


rhs = [0 0 0 0 t1+t2 t2 -t1-t2 -t2 t1+2*t2 t2 -t1-2*t2 -t2 1 1 1 1]';


F = set(lhs < rhs);
x = [x1 x2 x3 x4 x5];
sol = solvemp(F+set(-4< [t1 t2]<4),x1+x2+x3+x4,[],[t1;t2])

sdpvar x1 x2 y1 y2 t1 t2

F = set(x1+3*x2-y1 < 9-2*t1+t2);
F = F + set(2*x1+x2 < 8 + t1-2*t2);
F = F + set(x1<4+t1+t2);
F = F + set(-x1 < 0);
F = F + set(-x2 < 0);
F = F + set(binary([y1 y2])) + set(-10 <[t1 t2] < 10);


[sol,dgn,Z,J] = solvemp(F,-2*x1-x2+y1+y2,sdpsettings('mp.al',3),[t1;t2]);sol = rmovlps(sol);


sdpvar x1 x2 x3 x4 x5 x6 t1 t2

lhs = [-x1+x5;-x1-x5;-x2+x6;-x2-x6;-x3;-x3-x5;-x3;-x3+x5;-x4-x5-x6;-x4+x5;-x4+x5;-x4+x5+x6;x5;-x5;x6;-x6]


rhs = [0 0 0 0 t1+t2 t2 -t1-t2 -t2 t1+2*t2 t2 -t1-2*t2 -t2 1 1 1 1]';


F = set(lhs < rhs);
x = [x1 x2 x3 x4 x5];
sol = solvemp(F+set(-4< [t1 t2]<4),x1+x2+x3+x4,[],[t1;t2])



sdpvar x1 x2 t1 t2
binvar y1 y2

F = set(x1+3*x2-y1 < 9-2*t1+t2);
F = F + set(2*x1+x2 < 8 + t1-2*t2);
F = F + set(x1<4+t1+t2);
F = F + set(-x1 < 0);
F = F + set(-x2 < 0);
F = F + set(-10 <[t1 t2] < 10);

y = [y1;y2];
t = [t1;t2];
obj = -2*x1-x2+y1+y2;
sol = solvemp(F,obj,[],t);
sol = rmovlps(sol);

sold = dua_test(F,obj,sdpsettings('relax',1,'mp.presolve',1),t,y)
sold = rmovlps(sold);



sol = rmovlps(sol);
plot(sol.Pn)

sol = mimpmilp(F,obj,[],t)

sol = mimpmilp(F,-2*x1-x2+y1+y2,t,y);






sol = mimpmilp(F,obj,x{k},dd);



% Dua example 1
sdpvar x1 x2  t1 t2
binvar y1 y2

obj = -3*x1-8*x2+4*y1+2*y2

F = set([x1+x2 < 13+t1;5*x1-4*x2 < 20;-8*x1+22*x2 < 121+t2;4*x1+x2>8;x1-10*y1<0;x2-15*y2 < 0]);
F = F + set([x1 x2]>0) + set(0<[t1 t2]<10)
sol = solvemp(F,obj,[],[t1;t2]);



% Dua example 2
sdpvar x1 x2  t1 t2 t3
binvar y1 y2

obj = -3*x1-2*x2+10*y1+5*y2;

F=set([x1 < 10+t1+2*t2;
x2 < 10-t1 +t2;
x1+x2<20-t2;
x1+2*x2 < 12+t1-t3;
x1-20*y1<0;
x2-20*y2<0;
-x1+x2>4-t3;
y1+y2>1]);

F = F + set([x1 x2]>0) + set(0<[t1 t2 t3]<5)

sol = solvemp(F,obj,[],[t1;t2;t3]);

sol = rmovlps(sol)


% Dua 3
sdpvar x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23
sdpvar t1 t2 t3 t4
binvar y1 y2 y3 y4 y5 y6 y7 y8

x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23]

% x19 = 90000*y1+9600*x1;
% x20 = 40000*y2+8500*x2;
% x21 = 45000*(y3+y4) + 25*x17;
% x22 = 25000*(y5+y6) + 14.5*x18;
% x23 = 20000*(y7+y8)

F = set([x19 == 90000*y1+9600*x1;
x20 == 40000*y2+8500*x2;
x21 == 45000*(y3+y4) + 25*x17;
x22 == 25000*(y5+y6) + 14.5*x18;
x23 == 20000*(y7+y8);
x1  == x3+x4+x5+x6+x7;
x12 == x2+x3+x4+x5-x8-x9-x10;
x11 == x6+x7+x8+x9+x10;
x17 == x13 + x14;
x18 == x15+x16;
x13 == 42.6*x3+121*x6;
x14 == 42.6*x4+121*x7;
x15 == 78.4*x8;
x16 == 78.4*x9;
x12 > 20+10*t1;
x13+x15 > 7000+1000*t2;
x14+x16 > 4000+1000*t3;
x11 > 80+10*t4;
x6 - 200*y7 < 0;
x7 - 200*y8 < 0;
x3+x6-200*y3 < 0;
x4+x7-200*y4 < 0;
x8-200*y5 < 0;
x9-200*y6 < 0;
y3-y1<0;
y4-y1<0;
y7-y3<0;
y8-y4<0;
y3+y5 == 1;
y4+y6 == 1;
x > 0;
[t1 t2 t3 t4] <1;
[t1 t2 t3 t4] > 0])
obj = x19+x20+x21+x22+x23;
sol = solvemp(F,obj,[],[t1;t2;t3;t4]);
sol = rmovlps(sol);

sold = dua_test(F,obj,[],[t1;t2;t3;t4],[y1 y2 y3 y4 y5 y6 y7 y8])
sold = rmovlps(sold);

sol = solvesdp(F,obj)

y=[y1 y2 y3 y4 y5 y6 y7 y8]
for i = 1:23
solvesdp(F+set(x<1e7)+set(y==[1 0 0 1 1 0 0 0]),x(i),sdpsettings('relax',1));xl(i) = double(x(i));
solvesdp(F+set(x<1e7)+set(y==[1 0 0 1 1 0 0 0]),-x(i),sdpsettings('relax',1));xu(i) = double(x(i));
end








































[sol,dgn,Z,J] = solvemp(F,-2*x1-x2+y1+y2,[],t,y);


sol = rmovlps(sol);

binary = [y1;y2];

% Solve relaxed problem
[sol,dgn,Z,J] = solvemp(F,-2*x1-x2+y1+y2,sdpsettings('relax',1),[t1;t2],y);


not_integer = find((sum(abs([sol{1}.Fi{:}]),2) ~= 0) | sum(abs([sol{1}.Gi{:}]-round([sol{1}.Gi{:}])),2)~=0);

if ~isempty(not_integer)
    j = not_integer(1);
    [sol_down,dgn,Z,J] = solvemp(F+set(y(j)==0),-2*x1-x2+y1+y2,sdpsettings('relax',1),[t1;t2],y);
    [sol_up,dgn,Z,J] =   solvemp(F+set(y(j)==1),-2*x1-x2+y1+y2,sdpsettings('relax',1),[t1;t2],y);
    sol = rmovlps({sol_down{1},sol_up{1}})
end


