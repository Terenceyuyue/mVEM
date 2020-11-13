%test for mpmiqp with media workload
%4 proc,12 tasks,25 subtasks,n=4,m=12,nu=12,nx=4
%allocation matrix,n*m

%matrix initialization

n=3;
m=6;
allo_m = [0.1167 0 0.1625 0.05833 0 0.2 0 0.09 0.09 0 0 0;
    0.1833 0.09 0 0.0417 0.075 0.1375 0 0.15 0 0.075 0 0;
    0.1833 0 0.1125 0 0.065 0 0.175 0 0 0 0.1125 0;
    0.15 0.11 0 0.075 0.115 0 0.1083 0 0 0 0 0.0692];
allo_m = allo_m(1:n,1:m);
tp=pascal(n);

z = sdpvar(m,1);
x = sdpvar(n,1);

obj = (allo_m*z - x)'*(allo_m*z - x)

F = set(binary(z)) + set(allo_m*z <= x) + set(x <= 2*tp(:,1));

solvemp(F,obj,[],x)


n=4;
m=12;
allo_m = [0.1167 0 0.1625 0.05833 0 0.2 0 0.09 0.09 0 0 0;
    0.1833 0.09 0 0.0417 0.075 0.1375 0 0.15 0 0.075 0 0;
    0.1833 0 0.1125 0 0.065 0 0.175 0 0 0 0.1125 0;
    0.15 0.11 0 0.075 0.115 0 0.1083 0 0 0 0 0.0692];
%F=-2*allo_m
Matrices.F = -2*allo_m;
%Y=I,n*n
Matrices.Y = eye(n);
%G=allo_m
Matrices.G = allo_m;
%H=2*allo_m'*allo_m
Matrices.H = 2*(Matrices.G)'*Matrices.G;
%W=0,n*1
Matrices.W = 0*eye(n,1);
%E=I,n*n
Matrices.E = eye(n);
%Cf=Q2,assume zero here,1*m
Matrices.Cf = 0*eye(1,m);
%Cx=0,1*n
Matrices.Cx = 0*eye(1,n);
%Cc=0,1*1
Matrices.Cc = [0];
%bndA=I,n*n
Matrices.bndA = eye(n);
%bndb=2*[1;1;1...]
tp=pascal(n);
Matrices.bndb = 2*tp(:,1);

for i=1:m,
    Matrices.vartype(i)='B';
end

%call solver
Options.verbose = 1;
sol = mpt_mpmiqp(Matrices,Options);

%test
Pn=sol.Pn;
D0=[0.8;0.6;0.7;0.8];
reg=0;
mincost=inf;
count=0;
for i=1:length(Pn),
    [H{i},K{i}]=double(Pn(i));
    if H{i}*D0-K{i}<0
    	count=count+1;
        cost=D0'*sol.Ai{i}*D0 + sol.Bi{i}*D0 + sol.Ci{i};
        if cost<mincost
            mincost = cost;
            reg = i;
        end
    end
end
u = sol.Fi{reg}*D0+sol.Gi{reg};

