function fail = regress_bmibnb_gamsrobot(ops)

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
obj = t;

sol = solvesdp(F,obj,ops);

fail=getfail(sol.problem,double(obj),-0.375,checkset(F));
   
function fail =  getfail(problem,obj,objgoal,infeas)
fail = 0;
if problem == 0
    if abs(obj-objgoal)>1e-3
        fail = 1;
    elseif max(infeas)<-1e-6
        fail = 2;
    end
else
    fail = 3;
end

