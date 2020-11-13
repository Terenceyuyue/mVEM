sdpvar x y z
F = set(max(1,x)+max(y^2,z) < 3)+set(max(1,-min(x,y)) < 5)+set(norm([x;y],2) < z);
sol = solvesdp(F,max(x,z)-min(y,z)-z);

if ~(sol.problem == 0 & abs(double(max(x,z)-min(y,z)-z) - -sqrt(2))<1e-3)    
    error('Failed')    
end

sdpvar x y z
F = set(max(1,x)+max(y^2,z) < 3)+set(max(1,min(x,y)) < 5)+set(norm([x;y],2) < z);
sol = solvesdp(F,max(x,z)-min(y,z)-z);
if sol.problem ~= 14    
    error('Failed')    
end

randn('seed',456);
rand('seed',456);
A = randn(15,2);
b = rand(15,1)*5;
x = sdpvar(2,1);
sol = solvesdp([],-geomean(b-A*x)); 
if ~(sol.problem == 0 & abs(double(-geomean(b-A*x)) - -2.7797)<1e-3)    
    error('Failed')    
end

randn('seed',456);
rand('seed',456);
A = randn(15,2);
b = rand(15,1)*5;
x = sdpvar(2,1);
sol = solvesdp([],-geomean([b-A*x;min(x)])); 
if ~(sol.problem == 0 & abs(double(-geomean([b-A*x;min(x)])) - -2.4936)<1e-3)    
    error('Failed')    
end

n = 500;
randn('seed',456);
A = randn(2*n,n);
b = randn(2*n,1);
x=sdpvar(n,1);
sol = solvesdp([],norm(A*x-b));
if ~(sol.problem == 0 & abs(double(norm(A*x-b)) - 24.126)<1e-3)    
    error('Failed')    
end




