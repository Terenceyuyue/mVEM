function regress_dualize


ops = sdpsettings('verbose',0);

% TEST 1
yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
t = sdpvar(1,1);
y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(P > A*A') + set(P(3,3)>0) + set(t+y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>12);
obj = trace(P)+y;

fail  = regresstest(F,obj,ops);
regressreport('Test 1',fail)

% TEST 2
yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
t = sdpvar(1,1);
y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(P > 0) + set(P(3,3)>0) + set(t+y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>12);
obj = trace(P)+y;

fail  = regresstest(F,obj,ops);
regressreport('Test 2',fail)

yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
t = sdpvar(1,1);
y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(P > 0) + set(P(3,3)>0) + set(t+y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>0);
obj = trace(P)+y;

fail  = regresstest(F,obj,ops);
regressreport('Test 3',fail)

yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
t = sdpvar(1,1);
y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(P > 0) + set(P(3,3)>0) + set(t-y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>0);
obj = trace(P);

fail  = regresstest(F,obj,ops);
regressreport('Test 4',fail)

yalmip('clear')
X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = set(X>0) + set(cone(x(2:end),x(1))) + set(trace(X)==x(1)+2*x(2)+3*x(3)+4)+set(X(1,3)==8);
fail  = regresstest(F,obj,ops);
regressreport('Test 5',fail)

yalmip('clear')
X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = set(X>0) + set(cone(x(2:end),1+x(1))) + set(trace(X)==x(1)+2*x(2)+3*x(3)+4)+set(X(1,3)==8);
fail  = regresstest(F,obj,ops);
regressreport('Test 6',fail)

yalmip('clear')
X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = set(X>0) + set(cone(x(2:end),1+x(1))) + set(trace(X)==x(1)+2*x(2)+3*x(3)+4)+set(X(1,3)==8);
fail  = regresstest(F,obj,ops);
regressreport('Test 7',fail)

yalmip('clear')
X = sdpvar(3,3);
x = sdpvar(3,1);
obj = trace(X)+sum(x);
F = set(X>0) + set(cone(1-x(2:end),1+x(1))) + set(trace(X)==x(1)+2*x(2)+3*x(3)+4)+set(X(1,3)==8);
fail  = regresstest(F,obj,ops);
regressreport('Test 8',fail)

yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
t = sdpvar(1,1);
y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(P > A*A') + set(P(3,3)>0) + set(t+y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>12)+set(t>-12);
obj = trace(P)+y;
fail  = regresstest(F,obj,ops);
regressreport('Test 9',fail)

%yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
%t = sdpvar(1,1);
%y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(P > A*A') + set(P(3,3)>0) + set(t+y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>12)+set(t>-12);
obj = trace(P)+y+t;
fail  = regresstest(F,obj,ops);
regressreport('Test 10',fail)

yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
t = sdpvar(1,1);
y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(P > A*A') + set(P(3,3)>0) + set(P>0) + set(t+y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>12)+set(t>-12);
obj = trace(P)+y+t;
fail  = regresstest(F,obj,ops);
regressreport('Test 11',fail)

%yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
%t = sdpvar(1,1);
y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(P > A*A') + set(P(3,3)>0) + set(P>0) + set(t+y > 7) + set(t+y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>12)+set(t>-12);
obj = trace(P)+y+t;
fail  = regresstest(F,obj,ops);
regressreport('Test 12',fail)

yalmip('clear')
A = randn(3,3);A = -A*A';
P = sdpvar(3,3);
t = sdpvar(1,1);
y = sdpvar(1,1);
F = set(A'*P+P*A < -eye(3));
F = F + set(2*P > A*A') + set(P(3,3)>0) + set(P>0) + set(t+y > 7) + set(t+y > 7) + set(P(2,2)>4)+set(P(1,1:2)>t) + set(t>12)+set(t>-12);
obj = trace(P)+y+t;
fail  = regresstest(F,obj,ops);
regressreport('Test 13',fail)




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

    
function fail  = regresstest(F,obj,ops);
solvesdp(F,obj,ops);
obj1 = double(obj);
p1 = checkset(F);

x = recover(getvariables(F));
setsdpvar(x,0*double(x));

[Fdual,objdual,X,free,err]  = dualize(F,obj);
solvesdp(Fdual,-objdual,ops);

if length(X)>0
    for i = 1:length(X);setsdpvar(X{i},dual(Fdual(i)));end
end
if ~isempty(free)
    setsdpvar(free,dual(Fdual(end)));
end

obj2 = double(obj);
p2 = checkset(Fdual);

fail = 0;

if abs(obj1-obj2) > 1e-4
    fail = 1;
    [obj1 obj2]
end
if any(p2<-1e-5)
    fail = 2;
end

