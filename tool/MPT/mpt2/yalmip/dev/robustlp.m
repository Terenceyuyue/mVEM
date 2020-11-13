yalmip('clear')
sdpvar u v x(2,1)

% With decomposed problem
% min_x max_w c(w)'*x 
%
% s.t A(w)*x <= b(w)   for all E*w <=f
A = [-9 -2+5*u+v;3 -2;u 1-v;0 1]
c = [6-u;-7-u+v*8];
b = [4;1-v;9;6];
E = [eye(2);-eye(2)];
f = 0.3*[1;1;1;1];
solve_robust_lp(c,A,b,E,f,x,[u;v],set([]),sdpsettings)
double(x)

% Automatic decomposition!
w = [u;v];
F = set(A(1:end-1,:)*x < b(1:end-1)) + set(E*w < f) + set(x(2) < 6);
solverobust(F,c'*x,[],w);
double(x)

% Automatic decomposition!
w = [u;v];
F = set(A*x < b) + set(E*w < f);
solverobust(F,c'*x,sdpsettings('solver','sedumi'),w);
double(x)

% Automatic decomposition!
w = [u;v];
F = set(A*x < b) + set(E*w < f)
solverobust(F,c'*x,sdpsettings('solver','sedumi'),w);
double(x)



% By hand
sdpvar t
uu{1} = 0.3;
vv{1} = 0.3;
uu{2} = -0.3;
vv{2} = 0.3;
uu{3} = 0.3;
vv{3} = -0.3;
uu{4} = -0.3;
vv{4} = -0.3;
F = set([]);
for i = 1:4
    F = F + set([-9 -2+5*uu{i}+vv{i};3 -2;uu{i} 1-vv{i};0 1]*x < [4;1-vv{i};9;6]);
    F = F + set([6-uu{1} -7-uu{1}+vv{1}*8]*x < t);
end
solvesdp(F,t)
double(x)

