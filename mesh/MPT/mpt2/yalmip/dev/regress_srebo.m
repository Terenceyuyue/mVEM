randn('state',8342182);
n = 5;
m = 10;

S = randn(n,m)>0;
y = sign(randn(n,m)).*S;

% 1) Standard primal model, very inefficient
yalmip('clear')
X = sdpvar(n,m);
F = set(X(S).*y(S)>1);
obj = tracenorm(X);
solvesdp(F,obj);
X1 = double(X);

% 2) Try to dualize, will be better, but still not optimal
[Fd,objd] = dualize(F,obj);
assign(X,0*double(X))
solvesdp(Fd,-objd)
X2 = double(X);

% 3) Reduce problem after dualization (correspponds to imagemodel)
assign(X,0*double(X))
solvesdp(Fd,-objd,sdpsettings('remove',2));
X3 = double(X);

% 4) Manual
[i,a] = find(S);
q = sdpvar(1,length(i));
Q = sparse(i,a,q,n,m);
spQ = set([eye(n) (Q.*y);(Q.*y)' eye(m)]>=0); 
solvesdp(spQ + set(Q>=0),-sum(Q(S))); 
optAXXB = dual(spQ);
X4 = -optAXXB(1:n,(n+1):end)*2;




