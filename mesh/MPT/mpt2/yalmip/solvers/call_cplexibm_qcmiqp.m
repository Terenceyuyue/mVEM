function output = call_cplexibm_miqcp(interfacedata)

% Author Johan Löfberg
% $Id: call_cplexibm_miqcp.m,v 1.21 2009-11-03 11:08:47 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;

UB      = interfacedata.ub;
LB      = interfacedata.lb;

showprogress('Calling CPLEX',options.showprogress);

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [LB,UB,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,LB,UB);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(K.f+cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
end

n_original = length(c);
if K.q(1)>0
    % To simplify code, we currently normalize everything to z'*z<x0^2
    nNEW = sum(K.q);
    if ~isempty(x0)
        x0 = [x0;full(F_struc(1+K.f+K.l:end,:))*[1;x0]];
    end
        
    F_strucSOCP = [F_struc(1+K.f+K.l:end,:) -speye(nNEW)];
    F_struc = [F_struc(1:K.f+K.l,:) spalloc(K.f+K.l,nNEW,0)];
    UB = [UB;inf(nNEW,1)];
    c = [c;spalloc(nNEW,1,0)];
    Q = blkdiag(Q,spalloc(nNEW,nNEW,0));
    
    iCone = n_original+1;
    ri = zeros(1,length(K.q));
    Li = spalloc(n_original+nNEW,length(K.q),0);
    for i = 1:length(K.q);
        Qi{i} = sparse(iCone:iCone+K.q(i)-1,iCone:iCone+K.q(i)-1,[-1 ones(1,K.q(i)-1)],n_original+nNEW,n_original+nNEW);
        LB = [LB;0;-inf(K.q(i)-1,1)];
        iCone = iCone + K.q(i);
    end
    if i == 1
        Qi = Qi{1};
    end
    F_struc = [F_strucSOCP;F_struc];
    K.f = K.f + nNEW;
else
    Qi = [];
    ri = [];
    Li = [];
end

if K.l+K.f>0
    A =-(F_struc(1:K.f+K.l,2:end));
    B = full(F_struc(1:K.f+K.l,1));
end

if K.f > 0
    Aeq = A(1:K.f,:);
    beq = B(1:K.f);
else
    Aeq = [];
    beq = [];
end
if K.l > 0
    Aineq = A(1+K.f:K.f+K.l,:);
    bineq = B(1+K.f:K.f+K.l);
else
    Aineq = [];
    bineq = [];
end

ctype = char(ones(length(c),1)*67);
ctype(setdiff(integer_variables,semicont_variables)) = 'I';
ctype(binary_variables)  = 'B';  % Should not happen except from bmibnb
ctype(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
ctype(intersect(semicont_variables,integer_variables)) = 'N';

if nnz(Q)==0
    H = [];
else
    H = full(2*Q);
end

options.cplex.simplex.display = options.verbose;
options.cplex.mip.display = options.verbose;
options.cplex.barrier.display = options.verbose;

if options.savedebug
    save cplexintdebug H c Aineq bineq Aeq beq Li Qi ri LB UB ctype
end

% Call mex-interface
solvertime = clock;
if isempty(integer_variables) & isempty(binary_variables) & isempty(semicont_variables) & isempty(K.sos.type)
    [x,fval,exitflag,output] = cplexqcp(H, c(:), Aineq,bineq,Aeq,beq,Li,Qi,ri,LB,UB,x0,options.cplex);
else
    [x,fval,exitflag,output] = cplexmiqcp(H, c(:), Aineq,bineq,Aeq,beq,Li,Qi,ri,K.sos.type,K.sos.variables,K.sos.weight,LB,UB,ctype',x0,options.cplex);
end
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

if isempty(x)
    x = zeros(n_original,1);
else
    x = x(1:n_original);
end

problem = 0;
D_struc = [];

% Check, currently not exhaustive...
switch output.cplexstatus
    case {1,101,102}
        problem = 0;
    case {3,103,106}
        problem = 1; % Infeasible
    case {2,20,21,118}
        problem = 2; % Unbounded
    case 4
        problem = 1;
    case {10,11,104,105,107,108,111,112}
        problem = 3; % Iteration/time
    case {5,6,109,110}
        problem = 4; % Numerics
    case 119
        problem = 15;        
    otherwise
        problem = -1;
end

infostr = yalmiperror(problem,'CPLEX-IBM');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.H = H;
    solverinput.f = c(:);
    solverinput.Aineq = Aineq;
    solverinput.bineq = bineq;
    solverinput.Aeq = Aeq;
    solverinput.beq = beq;
    solverinput.ctype = ctype;
    solverinput.LB = LB;
    solverinput.UB = UB;
    solverinput.x0 = [];
    solverinput.options = options.cplex;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.FMIN = FMIN;
    solveroutput.SOLSTAT = SOLSTAT;
    solveroutput.DETAILS=DETAILS;
else
    solveroutput = [];
end


% Standard interface
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;

function nSOCP=countSOCP(F_struc,K)
nSOCP = 0;
top = K.f+K.l + 1;
ri = zeros(1,length(K.q));
Li = [];
for i = 1:length(K.q)
    % [cx+d;Ax+b]   |Ax+b|<cx+d, originally a QCQP
    m = K.q(i);
    ci = F_struc(top,2:end)';
    di = F_struc(top,1);
    Ai = F_struc(top+1:top+m-1,2:end);
    if(min(eig(full(Ai'*Ai - ci*ci')))<0)
        nSOCP = nSOCP + 1;
    end
    top = top+m;
end
