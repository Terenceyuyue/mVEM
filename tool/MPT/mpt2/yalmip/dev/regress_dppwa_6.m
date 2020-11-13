yalmip('clear')
%clear all

% Data
A = [2 -1;1 0];nx = 2;
B = [1;0];nu = 1;
C = [0.5 0.5];

% Prediction horizon
N = 4;

% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Binary for PWA selection
d = binvar(2,1);

% Value functions
J = cell(1,N);

% Initialize value function at stage N
J{N} = 0;
t = sdpvar(nx+nu,1);
bounds(t,0,600);
for k = N-1:-1:1
    
    bounds(x{k},-5,5);
    bounds(u{k},-1,1);
    bounds(x{k+1},-5,5);
    
    % Feasible region
    F =     set(-1 < u{k}     < 1);
    F = F + set(-1 < C*x{k}   < 1);
    F = F + set(-5 < x{k}     < 5);
    F = F + set(-1 < C*x{k+1} < 1);
    F = F + set(-5 < x{k+1}   < 5);

    % PWA Dynamics
    F = F + set(implies(d(1),x{k+1} == (A*x{k}+B*u{k})));
    F = F + set(implies(d(2),x{k+1} == (pi*A*x{k}+B*u{k})));
    F = F + set(implies(d(1),x{k}(1) > 0));
    F = F + set(implies(d(2),x{k}(1) < 0));
    F = F + set(sum(d) == 1);
    sdpvar v
    F = F + set(J{k+1}<v);
    
    F = F + set(-t < [x{k};u{k}] < t) ;

   [mpsol{k},sol{k},Uz{k},J{k}] = solvemp(F,sum(t) + v,[],x{k},u{k});
  % plot(mpsol{k}{1}.Pn)
  % pause
end

mpsol{1} = rmovlps(mpsol{1})

% Compare
sysStruct.A{1} = A;
sysStruct.B{1} = B;
sysStruct.C{1} = C;
sysStruct.D{1} = [0];
sysStruct.A{2} = A;
sysStruct.B{2} = B*pi;
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
% % Online
% obj = 0;
% F = set([]);
% dd = [];
% for k = N-1:-1:1
%     
%     bounds(x{k},-5,5);
%     bounds(u{k},-1,1);
%     bounds(x{k+1},-5,5);
%     
%     % Feasible region
%     F = F + set(-1 < u{k}     < 1);
%     F = F + set(-1 < C*x{k}   < 1);
%     F = F + set(-5 < x{k}     < 5);
%     F = F + set(-1 < C*x{k+1} < 1);
%     F = F + set(-5 < x{k+1}   < 5);
% 
%     % PWA Dynamics
%     d = binvar(2,1);dd = [dd;d];
%     F = F + set(implies(d(1),x{k+1} == (A*x{k}+B*u{k})));
%     F = F + set(implies(d(2),x{k+1} == (A*x{k}+pi*B*u{k})));
%     F = F + set(implies(d(1),x{k}(1) > 0));
%     F = F + set(implies(d(2),x{k}(1) < 0));
%     F = F + set(sum(d) == 1);
%    % F = F + set(-0.1 < u{k}-u{k+1} < 0.1);
%     obj = obj + norm([x{k};u{k}],1);
%     %obj = obj + x{k}'*x{k}+u{k}'*u{k};%norm([x{k};u{k}],1);
%     % Compute value function for one step backwards     
% end
% [mpsol2{k},sol{k},Uz{k},J2{k},U{k}] = solvemp(F,obj,[],x{k},u{k});
% solvesdp(F+set(x{k}==[0.5;0.5]),obj)
% solvesdp(F+set(x{k}==[1.2;0.8]),obj)
% mpsol{k} = solvemp(F,obj,[],x{k},u);
% mpsol{1} = rmovlps(mpsol{1});
% 
% 
