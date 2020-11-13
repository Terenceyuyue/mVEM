yalmip('clear')


% Data
A = [2 -1;1 0];nx = 2;
B = [1;0];nu = 1;
C = [0.5 0.5];

% Prediction horizon
N = 3;

% Future state
% Now for two different noises
x1  = sdpvar(repmat(nx,1,N),repmat(1,1,N));
x2 = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Current state
x  = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
v = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Binary for PWA selection
d = binvar(2,1);

% Value functions
J = cell(1,N);

% Initialize value function at stage N
J{N} = 0;
J1{N} = pwa(norm(x1{N},1),set(-10<x1{N}(1)<10));
J2{N} = pwa(norm(x2{N},1),set(-10<x2{N}(1)<10));

t = sdpvar(nx+nu,1);
bounds(t,0,600);
k = N-1
for k = N-1:-1:1
    
    bounds(x{k},-5,5);
    bounds(u{k},-1,1);
    bounds(x1{k+1},-5,5);
    bounds(x2{k+1},-5,5);
    
    % Feasible region
    F =     set(-1 < u{k}      < 1);
    F =     set(-1 < u{k}+v{k} < 1);
    F = F + set(-1 < C*x{k}    < 1);
    F = F + set(-5 < x{k}      < 5);
    
    F = F + set(-1 < C*x1{k+1} < 1);
    F = F + set(-1 < C*x2{k+1} < 1);
    F = F + set(-5 < x1{k}     < 5);    
    F = F + set(-5 < x2{k}     < 5);    

    % Two possible extreme predictions
    F = F + set(x1{k+1} ==  A*x{k}+B*u{k});
    F = F + set(x2{k+1} ==  pi*A*x{k}+B*(u{k}+v{k}));
       
    F = F + set(-t < [x{k};u{k}] < t) ;
    
    if k<N-1
        % Create two value functions, minimize worst case
        J1{k+1} = pwf(mpsol{k+1},x1{k+1},'convex');
        J2{k+1} = pwf(mpsol{k+1},x2{k+1},'convex');
        sdpvar w
        F = F + set(J1{k+1} < w) + set(J2{k+1} < w);
        obj = sum(t) + w;
    else
       % J1{N} = 0;
       % J2{N} = 0;
        sdpvar w
        F = F + set(J1{k+1} < w) + set(J2{k+1} < w);
        obj = sum(t)+w;        
    end
   [mpsol{k},sol{k},Uz{k},J{k}] = solvemp(F,obj,[],x{k},u{k});
end


mpsol{k} = rmovlps(mpsol{k})





J{N} = pwa(norm(x{N},1),set(-10<x{N}(1)<10));
t = sdpvar(nx+nu,1);
bounds(t,0,600);
k = N-1
for k = N-1:-1:1
  
    bounds(x{k},-5,5);
    bounds(u{k},-1,1);
    bounds(x{k+1},-5,5);    
    % Feasible region
    F =     set(-1 < u{k}     < 1);
    F = F + set(-1 < C*x{k}   < 1);
    F = F + set(-5 < x{k}     < 5);    
    % Two possible extreme predictions
    F = F + set(x{k+1} ==  A*x{k}+B*u{k});
    F = F + set(-t < [x{k};u{k}] < t) ;  
    obj = sum(t)+J{k+1};        
    [mpsol1{k}] = solvemp(F,obj,[],x{k},u{k});
end
    
J{N} = pwa(norm(x{N},1),set(-10<x{N}(1)<10));
t = sdpvar(nx+nu,1);
bounds(t,0,600);
k = N-1
for k = N-1:-1:1
  
    bounds(x{k},-5,5);
    bounds(u{k},-1,1);
    bounds(x{k+1},-5,5);    
    % Feasible region
    F =     set(-1 < u{k}     < 1);
    F = F + set(-1 < C*x{k}   < 1);
    F = F + set(-5 < x{k}     < 5);    
    % Two possible extreme predictions
    F = F + set(x{k+1} ==  pi*A*x{k}+B*u{k});
    F = F + set(-t < [x{k};u{k}] < t) ;  
    obj = sum(t)+J{k+1};        
    [mpsol2{k}] = solvemp(F,obj,[],x{k},u{k});
end
    
    
   
   


  J{N} = pwa(norm(x{N},1),set(-10<x{N}(1)<10));
    bounds(x{k},-5,5);
    bounds(u{k},-1,1);
    bounds(x{k+1},-5,5);    
    % Feasible region
    F =     set(-1 < u{k}     < 1);
    F = F + set(-1 < C*x{k}   < 1);
    F = F + set(-5 < x{k}     < 5);    
    % Two possible extreme predictions
    F = F + set(x{k+1} ==  pi*A*x{k}+B*u{k});
    F = F + set(-t < [x{k};u{k}] < t) ;
    
        obj = sum(t)+J{k+1};        
   [mpsol2{k}] = solvemp(F,obj,[],x{k},u{k});
   































break


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
