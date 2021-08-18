function regress_sos

ops1 = sdpsettings('solver','sdpt3','sos.cong',0,'sos.model',1,'verbose',0);
ops2 = sdpsettings('solver','sdpt3','sos.cong',1,'sos.model',2,'verbose',0);
ops3 = sdpsettings('solver','sdpt3','sos.cong',0,'sos.newton',0,'verbose',0,'sos.extlp',0);

if 1
sdpvar x s t;
F = set(sos(1+x+(1-s)*x^2-s))+set(sos(2+2*x+x^4-8*t))+set(s>0.49)+set(s+t<0.55);
obj= -s-t;

fail = regresstest(F,obj,ops1);
regressreport('Test 1',fail)
fail = regresstest(F,obj,ops2);
regressreport('Test 2',fail)
fail = regresstest(F,obj,ops3);
regressreport('Test 3',fail)


% Disjoint variables
x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
F = set(sos(1+y^2-t))+set(sos(1+x^2-t));
obj = -t;
fail = regresstest(F,obj,ops1);
regressreport('Test 4',fail)
fail = regresstest(F,obj,ops2);
regressreport('Test 5',fail)
fail = regresstest(F,obj,ops3);
regressreport('Test 6',fail)

% Disjoint variables and parameters
x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
F = set(sos(1+y^2-t))+set(sos(1+x^2-s));
obj = -s-t;
fail = regresstest(F,obj,ops1);
regressreport('Test 7',fail)
fail = regresstest(F,obj,ops2);
regressreport('Test 8',fail)
fail = regresstest(F,obj,ops3);
regressreport('Test 9',fail)


% Disjoint variables and parameters
x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
F = set(sos(1+y^2-t))+set(sos(1+x^2-s))+set(t>0);
obj = -s-t;
fail = regresstest(F,obj,ops1);
regressreport('Test 10',fail)
fail = regresstest(F,obj,ops2);
regressreport('Test 11',fail)
fail = regresstest(F,obj,ops3);
regressreport('Test 12',fail)

% Disjoint variables and parameters
x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
F = set(sos(1+y^2-t))+set(sos(1+x^2-s))+set(t>0)+set(s>0);
obj = -s-t;
fail = regresstest(F,obj,ops1);
regressreport('Test 13',fail)
fail = regresstest(F,obj,ops2);
regressreport('Test 14',fail)
fail = regresstest(F,obj,ops3);
regressreport('Test 15',fail)

x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
F = set(sos(1+y^2-t))+set(sos(1+x^2-s))+set(t>0)+set(s>0.5);
obj = -s-t;
fail = regresstest(F,obj,ops1);
regressreport('Test 16',fail)
fail = regresstest(F,obj,ops2);
regressreport('Test 17',fail)
fail = regresstest(F,obj,ops3);
regressreport('Test 18',fail)

x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
F = set(sos(1+y^2-t-s))+set(sos(1+x^2-s))+set(t>0)+set(s>-0.5);
obj = -s-t;
fail = regresstest(F,obj,ops1);
regressreport('Test 19',fail)
fail = regresstest(F,obj,ops2);
regressreport('Test 20',fail)
fail = regresstest(F,obj,ops3);
regressreport('Test 21',fail)


x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
F = set(sos(1+y^2-t-s))+set(sos(1+x^2-s))+set(t>0)+set(s>-0.5)+set(s>0);
obj = -s-t;
fail = regresstest(F,obj,ops1);
regressreport('Test 19',fail)
fail = regresstest(F,obj,ops2);
regressreport('Test 20',fail)
fail = regresstest(F,obj,ops3);
regressreport('Test 21',fail)

x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
F = set(sos(t*x^4+s*y^2-t))+set(sos(1+x^2-s));
obj = -t;
fail = regresstest(F,obj,ops1,[s t]);
regressreport('Test 22',fail)
fail = regresstest(F,obj,ops2,[s t]);
regressreport('Test 23',fail)
fail = regresstest(F,obj,ops3,[s t]);
regressreport('Test 24',fail)

x = sdpvar(1,1);
y = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
sdpvar u
F = set(sos(t*x^4+s*y^2-t))+set(sos(1+x^2-s));
obj = -t;
fail = regresstest(F,obj,ops1,[s t u]);
regressreport('Test 25',fail)
fail = regresstest(F,obj,ops2,[s t u]);
regressreport('Test 26',fail)
fail = regresstest(F,obj,ops3,[s t u]);
regressreport('Test 27',fail)


x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
t = sdpvar(1,1);
s = sdpvar(1,1);
w = sdpvar(1,1);
F = set(sos(x^4+z^6))+set(sos(x^2+(t+s+w-6)*x*z)) + set([w;s]>0)+set(t>3);
obj = w;
fail = regresstest(F,obj,ops1,t+s+w);
regressreport('Test 28',fail)
fail = regresstest(F,obj,ops2,t+s+w);
regressreport('Test 29',fail)
fail = regresstest(F,obj,ops3,t+s+w);
regressreport('Test 30',fail)

sdpvar x s t u
F = set(sos(1+x+16*s*x^2+13*u+t))+set(sos(2+2*x+(-8+u)*x^4+5-pi*t))+set(t>0)+set(t+u>0);
obj = t;
fail = regresstest(F,obj,ops1,t+s+u);
regressreport('Test 31',fail)
fail = regresstest(F,obj,ops2,t+s+u);
regressreport('Test 32',fail)
fail = regresstest(F,obj,ops3,t+s+u);
regressreport('Test 33',fail)

sdpvar x s t u
F = set(sos(1+x+16*s*x^2+13*u+t))+set(sos(2+2*x+(-8+u)*x^4+5-pi*t))+set(t>0)+set(t+u>0);
obj = t;
fail = regresstest(F,obj,ops1,t+s+u);
regressreport('Test 34',fail)
fail = regresstest(F,obj,ops2,t+s+u);
regressreport('Test 35',fail)
fail = regresstest(F,obj,ops3,t+s+u);
regressreport('Test 36',fail)

sdpvar x s t u
F = set(sos(1+x+16*s*x^2+13*u+t))+set(sos(2+2*x+(-8+u)*x^4+5-pi*t))+set(t>0)+set(t+u>0);
obj = s;
fail = regresstest(F,obj,ops1,t+s+u);
regressreport('Test 37',fail)
fail = regresstest(F,obj,ops2,t+s+u);
regressreport('Test 38',fail)
fail = regresstest(F,obj,ops3,t+s+u);
regressreport('Test 39',fail)

sdpvar x s u
F = set(sos(x+(t-100)*x^2+4*s-2+u))+set(s>0) + set(s<3) + set(u<4) + set(u>3) + set(t<102)+set(t>101.8)+set(t>0);
obj = -s-t;
fail = regresstest(F,obj,ops1,t+s+u);
regressreport('Test 40',fail)
fail = regresstest(F,obj,ops2,t+s+u);
regressreport('Test 41',fail)
fail = regresstest(F,obj,ops3,t+s+u);
regressreport('Test 42',fail)


sdpvar x s t u
F = set(sos(1+x+s*x^2))+set(sos(2+2*s*x+u*x^4));
obj = s;
fail = regresstest(F,obj,ops1,t+s+u);
regressreport('Test 43',fail)
fail = regresstest(F,obj,ops2,t+s+u);
regressreport('Test 44',fail)
fail = regresstest(F,obj,ops3,t+s+u);
regressreport('Test 45',fail)

end
sdpvar x y a b c
V = x^2*a+b*y^2+x*y+3;
F=set(sos(V))+set([a-4 b c]);
obj=a;
fail = regresstest(F,obj,ops1,a+b+c);
regressreport('Test 46',fail)
fail = regresstest(F,obj,ops2,a+b+c);
regressreport('Test 47',fail)
fail = regresstest(F,obj,ops3,a+b+c);
regressreport('Test 48',fail)


sdpvar x s t
F = set(sos(x^4+(s+t-2)*x^3+s+t));
obj = max(s,t);
fail = regresstest(F,obj,ops1,s+t);
regressreport('Test 49',fail)
fail = regresstest(F,obj,ops2,s+t);
regressreport('Test 50',fail)
fail = regresstest(F,obj,ops3,s+t);
regressreport('Test 51',fail)


% Markus tentacles problem
sdpvar x y a
f = x^4 * y^2 + x^2 * y^4 - 3 * x^2 * y^2 + 1, k = 0;
df = jacobian(f, [x y]);
g = 1 - (df(1)^2 + df(2)^2) * (x^2 + y^2);
if k >= 0
    v = monolist([x; y], 2*k);
    coeffVec = sdpvar(length(v), 1);
    t = coeffVec' * v;
    constraints = set(sos(f - a - t * g)) + set(sos(t));
else
    coeffVec = [];
    constraints = set(sos(f - a));
end
F = constraints;
obj = -a;
fail = regresstest(F,obj,ops1,[a; coeffVec]);
regressreport('Test 52',fail)
fail = regresstest(F,obj,ops2,[a; coeffVec]);
regressreport('Test 53',fail)
fail = regresstest(F,obj,ops3,[a; coeffVec]);
regressreport('Test 54',fail)




return



F = set(sos(x^4+z^6))+set(sos(x^2+(t+s+w-6)*x*z)) + set([w;s]>0)+set(t>3);
solvesos(F,[],sdpsettings('sos.traceobj',1,'sos.model',1,'solver','sedumi'),s+t+w)
double([s t w])
solvesos(F,[],sdpsettings('sos.traceobj',1,'sos.model',0,'solver','sedumi'),s+t+w)
double([s t w])

% Automatic detection of non-negative cones
sdpvar x s t u
F = set(sos(1+x+16*x^2+1));
solvesos(F,[],sdpsettings('sos.traceobj',1,'sos.solvedual',0,'solver','sedumi'))


% Failure before since no free
sdpvar x s t u
F = set(sos(x+(t-100)*x^2+4*s-2+u))+set(s>0) + set(s<3) + set(u<4) + set(u>3) + set(t<102)+set(t>101.8)+set(t>0);
solvesos(F,-s-t,sdpsettings('sos.traceobj',0,'sos.solvedual',1),s+t)
double([s t])

% Mixed constraint, all free
sdpvar x s t
F = set(sos(x+(t-100)*x^2+4*s-2))+set(s>1) + set(s<3) + set(t<102)+set(t>101.8)+set(2*s+t>104);
solvesos(F,s+t,sdpsettings('sos.traceobj',0,'sos.solvedual',1),s+t)
double([s t])

% Mixed constraint, one free
sdpvar x s t
F = set(sos(x+(t-100)*x^2+4*s-2))+set(s>1) + set(s>0) + set(s<3) + set(t<102)+set(t>101.8)+set(2*s+t>104);
solvesos(F,s+t,sdpsettings('sos.traceobj',0,'sos.solvedual',1),s+t)
double([s t])


% Mixed constraint, all free
sdpvar x s t
F = set(sos(x+(t-100)*x^2+4*s-2))+set(s>1) + set(t>0) + set(s<3) + set(t<102)+set(t>101.8)+set(2*s+t>104);
solvesos(F,s+t,sdpsettings('sos.traceobj',0,'sos.solvedual',1),s+t)
double([s t])



% Mixed constraint, all free
sdpvar x s t
F = set(sos(x+(t-100)*x^2+4*s-2))+set(1.001*s>1) + set(1.001*t>1) + set(s<3) + set(t<102)+set(1.01*t>101.8)+set(2*s+t>104);
solvesos(F,s+t,sdpsettings('sos.traceobj',0,'sos.solvedual',1),s+t)
double([s t])


% Mixed constraint, negative variables
sdpvar x s u t
F = set(sos(x+(t-100)*x^2+pi*s+6+u)) + set(t>1) + set(-12.34567<s<-1) + set(t<102)+set(t>101.8);
solvesos(F,s+t+u,sdpsettings('sos.traceobj',0,'sos.solvedual',1),s+t+u)
double([s t u])

% Mixed constraint, all free
sdpvar x s t u
F = set(sos(x+t*x^2+s+u)) + set(t>1.2) + set(s+t==4) + set(u<0.01) + set(u+t==0.1)
solvesos(F,s+t,sdpsettings('sos.traceobj',0,'sos.solvedual',1),s+t+u)
double([s t u])


% Mixed constraint, all constrained
sdpvar x s t u
F = set(sos(x+t*x^2+s+u)) + +set(u>0) + set(t==9) + set(s+t==11) + set(u==0.01)
solvesos(F,s+t,sdpsettings('sos.traceobj',0,'sos.solvedual',1),s+t+u)
double([s t u])

% 
% 
% % Failure : t is not used
sdpvar x s t u
F = set(sos(1+x+s*x^2))+set(sos(2+2*s*x+u*x^4));
solvesos(F,[],sdpsettings('sos.traceobj',1,'sos.solvedual',1,'solver','sedumi'),s+t+u)
double([s t u])





% Automatic detection of diagional constraints
sdpvar x s t u
F = set(sos(s+x+x^2))+set(sos(s+2*s*x+u*x^4));
solvesos(F,[],sdpsettings('sos.traceobj',1,'sos.solvedual',1,'solver','sedumi'),s+u)
double([s t u])



% Automatic detection of diagional constraints
sdpvar x s t u
F = set(sos(1+x^2+x^4));
solvesos(F,[],sdpsettings('sos.traceobj',1,'sos.solvedual',1,'solver','sedumi'))
double([s t u])


% Automatic detection of diagional constraints
sdpvar x y s t u
F = set(sos(1+x^2+x^4+x^2*y^2+y^4))+set(sos((7-t)*y^8+y^12));
solvesos(F,[],sdpsettings('sos.traceobj',1,'sos.solvedual',1,'solver','sedumi'),t)
double([s t u])


x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
p = x^4*y^2+x^2*y^4+z^6-3*x*x*y*y*z*z;
solvesos(set(sos(p)))
solvesos(set(sos(p)),[],sdpsettings('sos.model',2))

% PP, p 41 
x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
p = 2*x^4+2*x^3*y-x^2*y^2+5*y^4;
solvesos(set(sos(p)))
solvesos(set(sos(p)),[],sdpsettings('sos.model',2))

% PP, p 42
x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
p = x^4-(2*y*z+1)*x^2+(y^2*z^2+2*y*z+2)
solvesos(set(sos(p)))
solvesos(set(sos(p)),[],sdpsettings('sos.model',2))


% PP, p 43
x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
p = x^4+x^2+z^6-3*x*x*z*z;
solvesos(set(sos(p)))
solvesos(set(sos(p)),[],sdpsettings('sos.model',2))

x = sdpvar(1,1);
y = sdpvar(1,1);
z = sdpvar(1,1);
p = x^2+y^2+x^4+z^6+(x*x-y^3+z*z*y)^2
solvesos(set(sos(p)))
solvesos(set(sos(p)),[],sdpsettings('sos.model',2))


% robinson
x = sdpvar(1,1);
y = sdpvar(1,1);
p = x^6+y^6-x^4*y^2-y^4*x^2-x^4-y^4-x^2-y^2+3*x^2*y^2+10
solvesos(set(sos(p)))
solvesos(set(sos(p)),[],sdpsettings('sos.model',2))









function regressreport(text,fail)

switch fail
    case 0
        disp(['No problems in ' text]);
    case 1
        disp(['Objective differ in ' text]);
    case 2
        disp(['Infeasible solution in ' text]);
    otherwise
end

    
function fail  = regresstest(F,obj,ops,pv);

if nargin==3
    pv = [];
end

ops.sos.model = 1;
solvesos(F,obj,ops,pv);
obj1 = double(obj);
p1s = checkset(F(find(is(F,'sos'))));
p1e = checkset(F(find(~is(F,'sos'))));

ops.sos.model = 2;
solvesos(F,obj,ops,pv);
obj2 = double(obj);
p2s = checkset(F(find(is(F,'sos'))));
p2e = checkset(F(find(~is(F,'sos'))));


fail = 0;

if abs(obj1-obj2) > 1e-4
    fail = 1;
end

if any(p1s>1e-4)
   fail = 2;
   p1s
end
if any(p2s>1e-4)
   fail = 2;
   p2s
end
if any(p1e<-1e-4)
   fail = 2;
   p1e
end
if any(p2e<-1e-4)
   fail = 2;
   p2e
end
if fail==0
    disp('Correct solution');
end

