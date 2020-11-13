yalmip('clear')
clear all

% Prediction horizon
N = 5;

double_integrator

probStruct.N=N-1;
probStruct.R = 100;
probStruct.subopt_lev=0;
probStruct.y0bounds=1;
%probStruct.Tconstraint=0;
probStruct.P_N = zeros(2);
sysStruct.ymin = [-100;-100];
sysStruct.ymax = [100;100];

ctrl=mpt_control(sysStruct,probStruct);

M=mpt_constructMatrices(sysStruct,probStruct)

[H,K] = double(M.Pinvset);

nx = 2;
nu = 1;

% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Value functions
J = cell(1,N);

% Initialize value function at stage N
J{N} = 0;
sysStruct.xmin = sysStruct.ymin;
sysStruct.xmax = sysStruct.ymax;
F = set([]);
obj = 0;
for k = N-1:-1:1
    F = F + set(sysStruct.umin < u{k}     < sysStruct.umax);
    F = F + set(sysStruct.ymin < sysStruct.C*x{k}   < sysStruct.ymax);
    F = F + set(sysStruct.ymin < sysStruct.C*x{k+1} < sysStruct.ymax);  
    F = F + set(x{k+1} == sysStruct.A*x{k}+sysStruct.B*u{k});    
    obj = obj + norm([x{k};u{k}],1);%x{k}'*probStruct.Q*x{k}+u{k}*probStruct.R*u{k};    
end

[mpsol{k},sol{k},Uz{k},J{k}] = solvemp(F,obj,[],x{k},u{k});

%[pass,tol] = mpt_isPWAbigger(mpsol{1}{1},ctrl)