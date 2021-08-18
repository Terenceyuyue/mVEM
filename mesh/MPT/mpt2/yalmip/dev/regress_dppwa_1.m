yalmip('clear')
%clear all

% Data
A = [2 -1;1 0];nx = 2;
B1 = [1;0];
B2 = pi*B1;
nu = 1;
C = [0.5 0.5];

% Prediction horizon
N = 12;

% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Binary for PWA selection
d = binvar(repmat(2,1,N),repmat(1,1,N));

% Value functions
J = cell(1,N);

% Initialize value function at stage N
J{N} = 0;

% DP iteration from final state
obj = 0;
F = set([])
for k = N-1:-1:1

    % Control constraints
   % F = set(-1 < u{k} < 1);
    
    % We are in a feasible region...
    F = F + set(-1 < C*x{k}   < 1);
    F = F + set(-5 < x{k}     < 5);
    
    % ...and go to a feasible region
    F = F + set(-1 < C*x{k+1} < 1);
    F = F + set(-5 < x{k+1}   < 5);
 F = F + set(-5 < u{k}   < 5);
    % PWA Dynamics
    F = F + set(implies(d{k}(1),[x{k+1} == A*x{k}+B1*u{k}, x{k}(1) > 0]));
    F = F + set(implies(d{k}(2),[x{k+1} == A*x{k}+B2*u{k}, x{k}(1) < 0]));
    
%     F = F + set([1 0]*(A*x{k}+B1*u{k}) >= 100*(d{k}(1) + d{k+1}(1)-2));
%     F = F + set([-1 0]*(A*x{k}+B1*u{k}) >= 100*(d{k}(1) + d{k+1}(2)-2));
%     F = F + set([1 0]*(A*x{k}+B2*u{k}) >= 100*(d{k}(2) + d{k+1}(1)-2));
%     F = F + set([-1 0]*(A*x{k}+B2*u{k}) >= 100*(d{k}(2) + d{k+1}(2)-2));

    F = F + set(sum(d{k}) == 1);
    
    % L1 objective
  %  obj = norm(x{k},1) + norm(u{k},1);
  
    obj = obj + x{k}'*x{k} + u{k}^2;%norm(x{k},1) + norm(u{k},1),
 %  [mpsol{k},sol{k},aux{k},J{k},U{k}] = solvemp(F,obj + J{k+1},[],x{k},u{k});

end





break



break
mpsol{1} = mpt_removeOverlaps(mpsol{1})

% Compare
sysStruct.A{1} = A;
sysStruct.B{1} = B1;
sysStruct.C{1} = C;
sysStruct.D{1} = [0];
sysStruct.A{2} = A;
sysStruct.B{2} = B2;
sysStruct.C{2} = C;
sysStruct.D{2} = [0];
sysStruct.guardX{1} = [-1 0];
sysStruct.guardU{1} = [0];
sysStruct.guardC{1} = [0];
sysStruct.guardX{2} = [1 0];
sysStruct.guardU{2} = [0];
sysStruct.guardC{2} = [0];

%set constraints on output
sysStruct.ymin    =   -1;
sysStruct.ymax    =    1;

%set constraints on input
sysStruct.umin    =   -1;
sysStruct.umax    =   1;

sysStruct.xmin    =   [-5;-5];
sysStruct.xmax    =   [5;5];

probStruct.norm=1;
probStruct.Q=eye(2);
probStruct.R=1;
probStruct.N=N-1;
probStruct.P_N=zeros(2);
probStruct.subopt_lev=0;
probStruct.y0bounds=1;
probStruct.Tconstraint=0;
ctrl=mpt_control(sysStruct,probStruct)

mpt_isPWAbigger(ctrl,mpsol{1})
break
%[ii,jj] = isinside(ctrl.Pn,[1.2;0.8]);
%ctrl.Bi{jj}*[1.2;0.8]+ctrl.Ci{jj}

% 
% 
% Online
obj = 0;
F = set([]);
dd = [];
for k = N-1:-1:1
    
    % Feasible region
    F = F + set(-1 < u{k}     < 1);
    F = F + set(-1 < C*x{k}   < 1);
    F = F + set(-5 < x{k}     < 5);
    F = F + set(-1 < C*x{k+1} < 1);
    F = F + set(-5 < x{k+1}   < 5);

    % PWA Dynamics
    d = binvar(2,1);dd = [dd;d];
    F = F + set(implies(d(1),x{k+1} == (A*x{k}+B1*u{k})));
    F = F + set(implies(d(2),x{k+1} == (A*x{k}+B2*u{k})));
    F = F + set(implies(d(1),x{k}(1) > 0));
    F = F + set(implies(d(2),x{k}(1) < 0));
    F = F + set(sum(d) == 1);
    
    F1 = set(x{k+1} == (A*x{k}+B1*u{k})) + set(x{k}(1) > 0);
    F2 = set(x{k+1} == (A*x{k}+B2*u{k})) + set(x{k}(1) < 0);
    
    %F = F + hull(F1,F2);
    
   % obj = obj + norm([x{k};u{k}],1);
    obj = obj + norm([x{k};u{k}],1);
%    obj = obj + norm([x{k};u{k}],1);
    % Compute value function for one step backwards     
end
mpsol2{k} = solvemp(F,obj,[],x{k},u{k});
mpsol3{k} = dua_test(F,obj,[],x{k});
solvesdp(F+set(x{k}==[0.5;0.5]),obj)
solvesdp(F+set(x{k}==[1.2;0.8]),obj)
mpsol{k} = solvemp(F,obj,[],x{k},u);
mpsol2{1} = rmovlps(mpsol2{1});

mpsol1{k} = solvemp(F,obj,[],x{k},u{k});
% 
