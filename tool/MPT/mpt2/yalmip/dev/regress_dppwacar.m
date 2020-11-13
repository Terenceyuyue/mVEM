yalmip('clear')
clear all

% Prediction horizon
N = 4;

pwa_car
%sysStruct.xmin = sysStruct.ymin;
%sysStruct.xmax = sysStruct.ymax;
probStruct.R = 1;
probStruct.Q = eye(2);
probStruct.N=N-1;
probStruct.norm = 1;
probStruct.subopt_lev=0;
probStruct.P_N = zeros(2);
probStruct.y0bounds = 1;

if 0
    ctrl=mpt_control(sysStruct,probStruct)
    mpt_plotpwa(ctrl.Pn,ctrl.Bi,ctrl.Ci)
end
nx = 2;
nu = 1;

% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Binary for PWA selection
d = binvar(4,1);

% Value functions
J = cell(1,N);

% Initialize value function at stage N
J{N} = 0;
sysStruct.xmin = sysStruct.ymin;
sysStruct.xmax = sysStruct.ymax;
for k = N-1:-1:1
    % Feasible region
    t = sdpvar(nx+nu,1);
    bounds(x{k},sysStruct.xmin,sysStruct.xmax);
    bounds(u{k},sysStruct.umin,sysStruct.umax);
    bounds(x{k+1},sysStruct.xmin,sysStruct.ymax);
    bounds(t,0,40*2+5+1);

    F =     set(sysStruct.umin < u{k}     < sysStruct.umax);
    F = F + set(sysStruct.xmin < x{k}     < sysStruct.xmax);
    F = F + set(sysStruct.xmin < x{k+1}   < sysStruct.xmax);
    F = F + set(sysStruct.ymin < sysStruct.C{1}*x{k}   < sysStruct.ymax);
    F = F + set(sysStruct.ymin < sysStruct.C{1}*x{k+1} < sysStruct.ymax);

    F = F + set(-t < [x{k};u{k}] < t) ;

    % PWA Dynamics
    for i = 1:length(sysStruct.A)
        F = F + set(implies(d(i),x{k+1} == sysStruct.A{i}*x{k}+sysStruct.B{i}*u{k}+sysStruct.f{i}));
        F = F + set(implies(d(i),sysStruct.guardX{i}*x{k} <= sysStruct.guardC{i}));
    end
    F = F + set(sum(d) == 1);

    % Compute value function for one step backwards
    [mpsol{k},sol{k},Uz{k},J{k}] = solvemp(F,sum(t) + J{k+1},[],x{k},u{k});
end
break
mpsol{1} = rmovlps(mpsol{1})

%[pass,tol] = mpt_isPWAbigger(mpsol{1},ctrl)



break

% On-line solution
J{N} = 0;
sysStruct.xmin = sysStruct.ymin;
sysStruct.xmax = sysStruct.ymax;
F = set([]);
obj = 0;
nx = 2;
nu = 1;

% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Binary for PWA selection
d = binvar(4,1);
for k = N-1:-1:1
    % Feasible region
    t = sdpvar(nx+nu,1);
    d = binvar(4,1);
    obj = obj + sum(t);
    bounds(x{k},sysStruct.xmin,sysStruct.xmax);
    bounds(u{k},sysStruct.umin,sysStruct.umax);
    bounds(x{k+1},sysStruct.xmin,sysStruct.ymax);
    bounds(t,0,600);

    F = F + set(sysStruct.umin < u{k}     < sysStruct.umax);
    F = F + set(sysStruct.xmin < x{k}     < sysStruct.xmax);
    F = F + set(sysStruct.xmin < x{k+1}   < sysStruct.xmax);
    F = F + set(sysStruct.ymin < sysStruct.C{1}*x{k}   < sysStruct.ymax);
    F = F + set(sysStruct.ymin < sysStruct.C{1}*x{k+1} < sysStruct.ymax);

    F = F + set(-t < [x{k};u{k}] < t) ;

    % PWA Dynamics
    for i = 1:length(sysStruct.A)
        F = F + set(implies(d(i),x{k+1} == sysStruct.A{i}*x{k}+sysStruct.B{i}*u{k}+sysStruct.f{i}));
        F = F + set(implies(d(i),sysStruct.guardX{i}*x{k} <= sysStruct.guardC{i}));
    end
    F = F + set(sum(d) == 1);
    obj = obj + sum(t);
   
end
[mpsol{k},sol{k},Uz{k}] = solvemp(F,obj,[],x{k},u{k});
sol = solvesdp(F+set(x{k}==[-6;20]),obj)




