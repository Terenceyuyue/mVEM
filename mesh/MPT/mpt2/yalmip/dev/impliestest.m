
A1 = randn(8,2);
b1 = rand(8,1)*2-A1*[3;3];
A2 = randn(8,2);
b2 = rand(8,1)*2-A2*[-3;3];
A3 = randn(8,2);
b3 = rand(8,1)*2-A3*[3;-3];
A4 = randn(8,2);
b4 = rand(8,1)*2-A4*[-3;-3];

plot(polytope(A1,b1),polytope(A2,b2),polytope(A3,b3),polytope(A4,b4))

yalmip('clear')
binvar inp1 inp2 inp3 inp4

x = sdpvar(2,1);bounds(x,-100,100)
y = sdpvar(1);  bounds(y,-100,100);

F = set([]);
F = F + set(iff(inp1,A1*x < b1));
F = F + set(iff(inp2,A2*x < b2));
F = F + set(iff(inp3,A3*x < b3));
F = F + set(iff(inp4,A4*x < b4))
F = F + set(iff(inp4,y==pi));
F = F + set(inp1 | inp2 | inp3 | inp4);
%F = F + set(inp1 + inp2 + inp3 + inp4 == 1);
F=F+set(x==[-3;-4])

solvesdp(F,x(2))

F = set([]);
F = F + set(inp1 == (A1*x < b1));
F = F + set(inp2 == (A2*x < b2));
F = F + set(inp3 == (A3*x < b3));
F = F + set(inp4 == (A4*x < b4))
F = F + set(inp2 == (y==pi));
F = F + set(inp1 | inp2 | inp3 | inp4);
solvesdp(F,x(2))


binvar inp1 inp2 inp3 inp4

x = sdpvar(2,1);

F = set([]);
F = F + set(iff(A1*x < b1,inp1));
F = F + set(iff(A2*x < b2,inp2));
F = F + set(iff(A3*x < b3,inp3));
F = F + set(iff(A4*x < b4,inp4))

F = F + set(inp1 | inp2 | inp3 | inp4);


solvesdp(F,x(2))



x = sdpvar(2,1);bounds(x,-1000,1000)
F = set( (A1*x < b1) | (A2*x < b2) | (A3*x < b3) | (A4*x < b4));
solvesdp(F,x(2));


clear all
pwa_di
yalmip('clear')
N=3;
for i=1:N
    x{i} = sdpvar(2,1);bounds(x{i},-[5;5],[5;5]);
    u{i} = sdpvar(1,1);bounds(u{i},-1,1);
    y{i} = sdpvar(2,1);bounds(y{i},-[5;5],[5;5]);
end
delta = binvar(N,4);

F = set(x{1} == [0.1;0.11]);

for i = 1:N
    F = F + set(sysStruct.umin < u{i} < sysStruct.umax);
    F = F + set(sysStruct.ymin < y{i} < sysStruct.ymax);
    
    for j = 1:4
        F = F + set(iff(sysStruct.guardX{j}*x{i} < sysStruct.guardC{j},delta(i,j)));
        F = F + set(iff(delta(i,j),y{i} == [sysStruct.C{j} sysStruct.D{j} sysStruct.g{j}]*[x{i};u{i};1]));
    end
end
for i = 1:N-1
    for j = 1:4
        F = F + set(iff(delta(i,j),x{i+1} == [sysStruct.A{j} sysStruct.B{j} sysStruct.f{j}]*[x{i};u{i};1]));
    end
end

tu = sdpvar(N,1)
ty = sdpvar(N,1)
tx = sdpvar(N,1)
for i = 1:N
    F = F + set(-tu(i) < u{i} < tu(i))+set(-ty(i) < y{i} < ty(i))+set(-tx(i) < 8*(x{i}-0*[0.1;0.11]) < tx(i));
end
F=F+set(sum(delta,2)==1);
solvesdp(F,sum(2*tu)+sum(ty(1:end))+sum(tx(1:end)))





