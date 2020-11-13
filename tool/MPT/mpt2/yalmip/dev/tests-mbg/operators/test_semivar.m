function test_abs


A = magic(10);
b = A(:,end);
A = A(:,1:5);
x = 0.5*semivar(5,1);
e = b-A*x;
obj = norm(e,1);
sol = solvesdp(set(-100 < x < 100),obj);

mbg_asserttrue(sol.problem == 0)
mbg_asserttolequal(double(obj), 52.5784, 1e-4);

% Still not working
obj = e'*e
sol = solvesdp(set(-100 < x < 100),obj);
mbg_asserttrue(sol.problem == 0)
