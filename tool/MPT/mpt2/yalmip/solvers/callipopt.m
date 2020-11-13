function output = callipopt(model)

% Author Johan Löfberg
% $Id: callipopt.m,v 1.6 2009-09-29 10:30:15 joloef Exp $

options = model.options;
model = yalmip2nonlinearsolver(model);

if ~model.derivative_available
    disp('Derivate-free call to ipopt not yet implemented')
    error('Derivate-free call to ipopt not yet implemented')
end
if model.options.savedebug
    save ipoptdebug model
end
showprogress('Calling IPOPT',model.options.showprogress);

Fupp = [ repmat(0,length(model.bnonlinineq),1);
    repmat(0,length(model.bnonlineq),1);
    repmat(0,length(model.b),1);
    repmat(0,length(model.beq),1)];

Flow = [ repmat(-inf,length(model.bnonlinineq),1);
    repmat(0,length(model.bnonlineq),1);
    repmat(-inf,length(model.b),1);
    repmat(0,length(model.beq),1)];

usrf = 'ipopt_callback_f';
usrdf = 'ipopt_callback_df';
usrg = 'ipopt_callback_g';
usrdg = 'ipopt_callback_dg';

if length(model.bnonlineq) == 0
    jac_c_constant = 'yes';
else
    jac_c_constant = 'no';
end
if length(model.bnonlinineq) == 0
    jac_d_constant = 'yes';
else
    jac_d_constant = 'no';
end

fields = fieldnames(options.ipopt);
values = struct2cell(options.ipopt);
ops = cell(2*length(fields),1);
for i = 1:length(fields)
    ops{2*i-1} = fields{i};
    ops{2*i} = values{i};
end
ops{end+1} = 'print_level';
ops{end+1} = 3*options.verbose;

% Since ipopt react strangely on lb>ub, we should bail if that is detected
% (ipopt creates an exception)
if ~isempty(model.lb)
    if any(model.lb<model.ub)
        problem = 1;   
        solverinput = [];
        solveroutput = [];  
        output = createoutput(model.lb*0,[],[],problem,'IPOPT',solverinput,solveroutput,0);
    end
end

% These are needed to avoid recomputation due to ipopts double call to get
% f and df, and g and dg
global latest_x_f
global latest_x_g
global latest_df
global latest_f
global latest_G
global latest_g
latest_G= [];
latest_g = [];
latest_x_f = [];
latest_x_g = [];


solvertime = clock;

funcs.objective = @(x)ipopt_callback_f(x,model);
funcs.gradient = @(x)ipopt_callback_df(x,model);
funcs.constraints = @(x)ipopt_callback_g(x,model);
funcs.jacobian  = @(x)ipopt_callback_dg(x,model);

options.lb = model.lb;
options.ub = model.ub;
options.cl = Flow;
options.cu = Fupp;
options.ipopt.print_level = options.verbose;

n = length(Flow);
m = length(model.lb);
funcs.jacobianstructure = @() sparse(ones(n,m));
%eval(['@(x)sparse(ones(' num2str(n) ',' num2str(m) '));']);
solvertime = clock;
[xout,info] = ipopt(model.x0,funcs,options);
solvertime = etime(clock,solvertime);

% Duals currently not supported
lambda = [];

x = RecoverNonlinearSolverSolution(model,xout);

switch info.status
    case {0,1}
        problem = 0;
    case {2}
        problem = 1;
    case {-1}
        problem = 3;
    case {3,4,-2,-3}
        problem = 4;
    case {-11,-12,-13}
        problem = 7;
    case {-10,-100,-101,-102,-199}
        problem = 11;
    otherwise
        problem = -1;
end

% Internal format for duals
D_struc = [];

% Save all data sent to solver?
if options.savesolverinput
    solverinput.model = model;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = xout;
    solveroutput.lambda = lambda;
    solveroutput.iters = iters;
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'IPOPT',solverinput,solveroutput,solvertime);