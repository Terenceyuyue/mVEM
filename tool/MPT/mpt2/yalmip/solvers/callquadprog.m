function output = callquadprog(interfacedata)

% Author Johan Löfberg
% $Id: callquadprog.m,v 1.17 2007-08-02 11:39:36 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

[Q,c,A,b,Aeq,beq,lb,ub,ops] = yalmip2quadprog(interfacedata);

if options.savedebug
    %  ops = options.quadprog;
    save quadprogdebug Q c A b Aeq beq lb ub x0 ops
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock;
if nnz(Q) == 0
    % BUG in LIN/QUADPROG, computation of lambda crashes in some rare
    % cases. To avoid seeing this when we don't want the lambdas anyway, we
    % don't ask for it
    if options.saveduals
        [x,fmin,flag,output,lambda] = linprog(c, A, b, Aeq, beq, lb, ub, x0,ops);
    else
        lambda = [];
        [x,fmin,flag,output] = linprog(c, A, b, Aeq, beq, lb, ub, x0,ops);
    end
else
    if options.saveduals
        [x,fmin,flag,output,lambda] = quadprog(2*Q, c, A, b, Aeq, beq, lb, ub, x0,ops);
    else
        lambda = [];
        [x,fmin,flag,output] = quadprog(2*Q, c, A, b, Aeq, beq, lb, ub, x0,ops);
    end
    if flag==5
        [x,fmin,flag,output,lambda] = quadprog(2*Q, c, A, b, Aeq, beq, lb, ub, [],ops);
    end
end
%etime(clock,solvertime)
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

problem = 0;

% Internal format for duals
if ~isempty(lambda)
    D_struc = [lambda.eqlin;lambda.ineqlin];
else
    D_struc = [];
end

% Check, currently not exhaustive...
if flag==0
    problem = 3;
else
    if flag==-2
        problem = 1;
    else
        if flag>0
            problem = 0;
        else
            if isempty(x)
                x = repmat(nan,length(c),1);
            end
            if any((A*x-b)>sqrt(eps)) | any( abs(Aeq*x-beq)>sqrt(eps))
                problem = 1; % Likely to be infeasible
            else
                if c'*x<-1e10 % Likely unbounded
                    problem = 2;
                else          % Probably convergence issues
                    problem = 5;
                end
            end
        end
    end
end
infostr = yalmiperror(problem,'QUADPROG');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aeq;
    solverinput.beq = beq;
    solverinput.c = c;
    solverinput.H = Q;
    solverinput.options = options.quadprog;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
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