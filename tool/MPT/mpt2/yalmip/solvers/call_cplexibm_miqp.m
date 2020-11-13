function output = call_cplexibm_miqp(interfacedata)

% Author Johan Löfberg
% $Id: call_cplexibm_miqp.m,v 1.4 2010-03-24 19:49:16 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;

ub      = interfacedata.ub;
lb      = interfacedata.lb;

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [lb,ub,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,lb,ub);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(K.f+cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
end

showprogress('Calling CPLEX-IBM',options.showprogress);

% Notation used
H = 2*interfacedata.Q;
f = c;
if ~isempty(F_struc)
    Aineq = -F_struc(:,2:end);
    bineq = F_struc(:,1);
else
    Aineq = [];
    bineq = [];
end

if K.f > 0
    Aeq = Aineq(1:K.f,:);
    beq = bineq(1:K.f);
    Aineq(1:K.f,:)=[];
    bineq(1:K.f)=[];
else
    Aeq = [];
    beq = [];
end

if all(isinf(lb))
    lb = [];
end
if all(isinf(ub))
    ub = [];
end

ctype = char(ones(length(f),1)*67);
ctype(setdiff(integer_variables,semicont_variables)) = 'I';
ctype(binary_variables)  = 'B';  % Should not happen except from bmibnb
ctype(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
ctype(intersect(semicont_variables,integer_variables)) = 'N';

solvertime = clock;
options.cplex.simplex.display = options.verbose;
options.cplex.mip.display = options.verbose;
options.cplex.barrier.display = options.verbose;


if options.savedebug
    save cplexdebug
end

if isempty(integer_variables) & isempty(binary_variables) & isempty(semicont_variables) & isempty(K.sos.type)
    if isempty(Aineq) & isempty(Aeq) & isempty(lb) & isempty(ub)
        [x,fval,exitflag,output,lambda] = cplexqp(H,f,zeros(1,length(f)),1,Aeq,beq,lb,ub,x0,options.cplex);
        if ~isempty(lambda)
            lambda.ineqlin = [];
        end
    else
        [x,fval,exitflag,output,lambda] = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);
    end
else
    if isempty(Aineq) & isempty(Aeq) & isempty(lb) & isempty(ub)
        [x,fval,exitflag,output] = cplexmiqp(H,f,zeros(1,length(f)),1,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype',x0,options.cplex);
    else
        [x,fval,exitflag,output] = cplexmiqp(H,f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype',x0,options.cplex);
    end
    lambda = [];
end
solvertime = etime(clock,solvertime);

if ~isempty(lambda)
    D_struc = [-lambda.eqlin;-lambda.ineqlin];
else
    D_struc = [];
end

% Check, currently not exhaustive...
switch output.cplexstatus
    case {1,101,102}
        problem = 0;
    case {3,103,106}
        problem = 1; % Infeasible
    case {2,20,21,118}        
        problem = 2; % Unbounded
    case 4
        
        if isempty(integer_variables) & isempty(binary_variables)
    if isempty(Aineq) & isempty(Aeq) & isempty(lb) & isempty(ub)
        [x,fval,exitflag,output,lambda] = cplexqp(H*0,f*0,zeros(1,length(f)),1,Aeq,beq,lb,ub,[],options.cplex);
        if ~isempty(lambda)
            lambda.ineqlin = [];
        end
    else
        [x,fval,exitflag,output,lambda] = cplexqp(H*0,f*0,Aineq,bineq,Aeq,beq,lb,ub,[],options.cplex);
    end
else
    if isempty(Aineq) & isempty(Aeq) & isempty(lb) & isempty(ub)
        [x,fval,exitflag,output] = cplexmiqp(H*0,f*0,zeros(1,length(f)),1,Aeq,beq,[],[],[],lb,ub,ctype',[],options.cplex);
    else
        [x,fval,exitflag,output] = cplexmiqp(H*0,f*0,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype',[],options.cplex);
    end
    lambda = [];
        end
if output.cplexstatus == 3
    problem = 1;
else
    problem = 2;
end
   %     problem = 15;
    case {10,11,104,105,107,108,111,112}
        problem = 3; % Iteration/time
    case 119
        problem = 15;
    case {5,6,109,110}
        problem = 4; % Numerics
    otherwise
        problem = -1;
end

infostr = yalmiperror(problem,'CPLEX-IBM');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.H = H;
    solverinput.f = f;
    solverinput.Aineq = Aineq;
    solverinput.Aineq = Aeq;
    solverinput.bineq = bineq;
    solverinput.beq = beq;    
    solverinput.LB = lb;
    solverinput.UB = ub;
    solverinput.X0 = [];
    solverinput.ctype = ctype;
    solverinput.options = options.cplex;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.lambda = lambda;
else
    solveroutput = [];
end

if isempty(x)
    x = zeros(length(f),1);
end
% Standard interface
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;