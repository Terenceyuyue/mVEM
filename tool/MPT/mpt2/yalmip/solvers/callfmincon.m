function output = callfmincon(model)

% Author Johan Löfberg
% $Id: callfmincon.m,v 1.61 2010-01-20 10:20:57 joloef Exp $


model = yalmip2nonlinearsolver(model);

switch model.options.verbose
    case 0
        model.options.fmincon.Display = 'off';
    case 1
        model.options.fmincon.Display = 'final';
    otherwise
        model.options.fmincon.Display = 'iter';
end

if isfield(model.options.fmincon,'LargeScale')
    if isequal(model.options.fmincon.LargeScale,'off')
        model.A = full(model.A);
        model.b = full(model.b);
        model.Aeq = full(model.Aeq);
        model.beq = full(model.beq);
    end
end

if model.derivative_available
    model.options.fmincon.GradObj = 'on';    
end
if model.derivative_available
    model.options.fmincon.GradConstr = 'on';   
end

if model.options.savedebug
    ops = model.options.fmincon;
    save fmincondebug model %A b Aeq beq x0 lb ub ops
end

showprogress('Calling FMINCON',model.options.showprogress);

warning('off','optim:fmincon:NLPAlgLargeScaleConflict')
solvertime = clock;
[xout,fmin,flag,output,lambda] = fmincon('fmincon_fun',model.x0,model.A,model.b,model.Aeq,model.beq,model.lb,model.ub,'fmincon_con',model.options.fmincon,model);
solvertime = etime(clock,solvertime);
warning('on','optim:fmincon:NLPAlgLargeScaleConflict')

x = RecoverNonlinearSolverSolution(model,xout);

problem = 0;

% Internal format for duals
D_struc = [];

% Check, currently not exhaustive...
if flag==0
    problem = 3;
else
    if flag>0
        problem = 0;
    else
        if isempty(x)
            x = repmat(nan,length(model.c),1);
        end
        if model.c'*x<-1e10 % Likely unbounded
            problem = 2;
        else          % Probably convergence issues
            problem = 5;
        end
    end
end

% Save all data sent to solver?
if model.options.savesolverinput
    solverinput.model = model;
    solverinput.A = model.A;
    solverinput.b = model.b;
    solverinput.Aeq = model.Aeq;
    solverinput.beq = model.beq;
    solverinput.options = model.options.fmincon;
else
    solverinput = [];
end

% Save all data from the solver?
if model.options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'FMINCON',solverinput,solveroutput,solvertime);