function output = callqpip(interfacedata)

% Author Johan Löfberg
% $Id: callqpip.m,v 1.4 2007-05-22 13:40:36 joloef Exp $

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
    ops = options.qpip;
    save qpipdebug Q c A b Aeq beq lb ub x0 ops
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock;
[x,flag,lm] = qpip(2*Q, c, A, b, Aeq, beq, lb, ub, options.verbose,options.qpip.mu,options.qpip.method);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% Internal format for duals
if ~isempty(lm)
    D_struc = [lm.equality;lm.inequality];
else
    D_struc = [];
end

% Check, currently not exhaustive...
switch flag
    case 0
        problem = 0;
    case 2
        problem = 1;
    otherwise
        problem = 1;
end

infostr = yalmiperror(problem,'QPIP');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aq;
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
    solveroutput.flag = flag;
    solveroutput.lm = lm;
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