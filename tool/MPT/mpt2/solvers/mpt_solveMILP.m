function [xmin,fmin,how,exitflag]=mpt_solveMILP(f,A,B,Aeq,Beq,lb,ub,vartype,param,options,solver)
%MPT_SOLVEMILP Interface to various MILP solvers
%
% [xmin,fmin,how,exitflag]=mpt_solveMILP(f,A,B,Aeq,Beq,lb,ub,vartype,param,options,solver)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Solves an MILP problem:
%
%     min  f'x
%     s.t. A*x  <= B
%          Aeq*x = Beq
%          some 'x' integer/boolean
%
% by using the method specified in 'solver'
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% f          An (n x 1) vector containing the linear objective function coefficients.
%            REQUIRED INPUT ARGUMENT.
% 
% A          An (m x n) matrix (in full or sparse format) containing the constraint
%            coefficients. REQUIRED INPUT ARGUMENT.
% 
% b          An (m x 1) vector containing the right-hand side value for each
%            constraint in the constraint matrix. REQUIRED INPUT ARGUMENT.
% 
% Aeq        An (k x n) matrix (in full or sparse format) containing the constraint
%            coefficients for equality constraints (i.e. Aeq*x = Beq)
% 
% beq        An (k x 1) vector containing the right-hand side value for each
%            constraint in the constraint matrix for equality constraints
%
% LB         An (n x 1) vector containing the lower bound on each of the variables.
%            Any lower bound that is set to a value less than or equal to that of
%            the constant -CPX_INFBOUND will be treated as negative \infty.
%            CPX_INFBOUND is defined in the header file cplex.h.
%            Default: [], (lower bound of all variables set to -CPX_INFBOUND).
% 
% UB         An (n x 1) vector containing the upper bound on each of the variables.
%            Any upper bound that is set to a value greater than or equal to that of
%            the constant CPX_INFBOUND will be treated as \infty.
%            CPX_INFBOUND is defined in the header file cplex.h.
%            Default: [], (upper bound of all variables set to CPX_INFBOUND).
%
% VARTYPE    An (n x 1) vector containing the types of the variables
%            VARTYPE(i) = 'C' Continuous variable
%            VARTYPE(i) = 'B' Binary(0/1) variable
%            VARTYPE(i) = 'I' Integer variable
%            VARTYPE(i) = 'S' Semi-continuous variable
%            VARTYPE(i) = 'N' Semi-integer variable
%            (This is case sensitive).
%            Default: [], (all variables are continuous).
%
% solver     which solver to use:
%              0 - CPLEX 9 (cplexint)
%              1 - YALMIP
%              2 - GLPK (glpkmex)
%              3 - XPRESS
%              4 - MOSEK
%              5 - bintprog
%              6 - CPLEX 8 (cplexmex)
%              7 - CPLEXMEX (by Nicolo Giorgetti)
%
% Note: if 'solver' is not specified, mptOptions.milpsolver will be used instead
%       (see help mpt_init)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% xopt      - The optimizer
% fmin      - Value of the objective
% exitflag  - An integer specifying result of the optimization
%             (1 - optimal solution found, -1 - problem is infeasible)
% how       - States the result of optimization ('ok', 'unbounded', 'infeasible')
%
% see also MPT_SOLVELP, MPT_SOLVEMIQP

% Copyright is with the following author(s):
%
%(C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%              kvasnica@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This program is free software; you can redistribute it and/or
%          modify it under the terms of the GNU General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This program is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          General Public License for more details.
%
%          You should have received a copy of the GNU General Public
%          License along with this library; if not, write to the
%          Free Software Foundation, Inc.,
%          59 Temple Place, Suite 330,
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

global mptOptions

error(nargchk(3,11,nargin));

if ~isstruct(mptOptions)
    mpt_error;
end

% set default values for parameters which have not been specified
if nargin <= 3,
    Aeq = [];
    Beq = [];
    lb = [];
    ub = [];
    vartype = [];
    param = [];
    options = [];
    solver = mptOptions.milpsolver;
elseif nargin <= 4,
    error('mpt_solveMILP: you must specify ''Beq'' as well!');
elseif nargin <= 5,
    lb = [];
    ub = [];
    vartype = [];
    param = [];
    options = [];
    solver = mptOptions.milpsolver;
elseif nargin <= 6,
    error('mpt_solveMILP: you must specify ''ub'' as well!');    
elseif nargin <= 7,
    vartype = [];
    param = [];
    options = [];
    solver = mptOptions.milpsolver;
elseif nargin <= 8,
    param = [];
    options = [];
    solver = mptOptions.milpsolver;
elseif nargin <= 9,
    options = [];
    solver = mptOptions.milpsolver;
elseif nargin <= 10,
    solver = mptOptions.milpsolver;
end
f = f(:);
lb = lb(:);
ub = ub(:);
vartype = vartype(:);

% CPLEXINT doesn't like options.verbose=-1
if isfield(options, 'verbose'),
    options.verbose = max(options.verbose, 0);
end

if solver==0
    % use cplex 9
    
    % use initial guess for xmin if provided.
    % note that CPLEXINT uses a special syntax - x0 must be a 2 column matrix
    % where the first column represents indices of elements for which an
    % appropriate guess is provided in the second column
    if isfield(options, 'usex0'),
        options.x0 = [[1:length(options.usex0)]' options.usex0];
    end
    
    indeq = [];
    if ~isempty(Aeq),
        nc = size(A,1);
        A = [A; Aeq];
        B = [B; Beq];
        indeq = (nc+1:size(A,1))';
    end
    [xmin,fmin,status,details]=cplexint([], f, A, B, indeq, [], lb, ub, vartype, param, options);
    how = lower(details.statstring);
    exitflag = -1;
    if strcmp(how, 'optimal') | strcmp(how, 'optimaltol') | ...
            strcmp(how, 'integer optimal solution') | strcmp(how, 'integer optimal, tolerance'),
        how = 'ok';
        exitflag = 1;
    end
    
elseif solver==1
    % YALMIP / BNB
    [xmin,fmin,how,exitflag]=yalmipMILP(f,A,B,Aeq,Beq,lb,ub,vartype,param,options,'bnb');
    
elseif solver==2
    % use glpk
    [nc, nx] = size(A);
    indeq = [];
    if ~isempty(Aeq),
        A = [A; Aeq];
        B = [B; Beq];
        indeq = (nc+1:size(A,1))';
        nc = size(A,1);
    end
    if isempty(lb),
        lb = -1e9*ones(nx, 1);
    end
    if isempty(ub),
        ub = 1e9*ones(nx, 1);
    end
    if length(lb)~=nx | length(ub)~=nx,
        error('Wrong size of LB and/or UB.');
    end
    % change 'B' (boolean) vartype to 'I' (integer)
    % glpkmex does not support 'B' type variables
    bpos = find(vartype=='B');
    vartype(bpos) = 'I';
    lb(bpos) = 0;
    ub(bpos) = 1;
    
    SENSE = 1; % minimize
    CTYPE = repmat('U',nc,1); % all constraints <=
    CTYPE(indeq) = 'S';       % change requested constraints to equality
    param.msglev = 0;
    param.itlim = -1;
    if isfield(options,'glpk_solver'),
        lpsolver = options.glpk_solver;
        [xmin,fmin,status,details]=glpkmex(SENSE, f, A, B, CTYPE, lb, ub, ...
            vartype, param, lpsolver);
    else
        [xmin,fmin,status,details]=glpkmex(SENSE, f, A, B, CTYPE, lb, ub, ...
            vartype, param);
    end
    
    switch status
        case {171, 180, 181, 151, 200}
            how = 'ok';
            exitflag = 1;
        otherwise
            how = 'infeasible';
            exitflag = -1;
    end

elseif solver==8
    % use glpkcc
    [nc, nx] = size(A);
    indeq = [];
    if ~isempty(Aeq),
        A = [A; Aeq];
        B = [B; Beq];
        indeq = (nc+1:size(A,1))';
        nc = size(A,1);
    end
    if isempty(lb),
        lb = -1e9*ones(nx, 1);
    end
    if isempty(ub),
        ub = 1e9*ones(nx, 1);
    end
    if length(lb)~=nx | length(ub)~=nx,
        error('Wrong size of LB and/or UB.');
    end
    % change 'B' (boolean) vartype to 'I' (integer)
    % glpkmex does not support 'B' type variables
    bpos = find(vartype=='B');
    vartype(bpos) = 'I';
    lb(bpos) = 0;
    ub(bpos) = 1;
    
    lb(isinf(lb)) = -1e9;
    ub(isinf(ub)) = 1e9;
    B(B==Inf) = 1e9;
    B(B==-Inf) = -1e9;
    
    SENSE = 1; % minimize
    CTYPE = repmat('U',nc,1); % all constraints <=
    CTYPE(indeq) = 'S';       % change requested constraints to equality
    param.msglev = 0;
    param.itlim = -1;
    if isfield(options,'glpk_solver'),
        param.lpsolver = options.glpk_solver;
    end
    propagate_fields = {'save', 'savefilename', 'savefiletype', ...
        'itlim', 'presolve', 'itcnt', 'tmlim', 'msglev', ...
        'tolint', 'tolobj', 'msglev', 'tolbnd', 'branch', 'presol', ...
        'scale', 'dual', 'round', 'btrack'};
    for i = 1:length(propagate_fields)
        vv = propagate_fields{i};
        if isfield(options, vv)
            param = setfield(param, vv, getfield(options, vv));
        end
    end

    [xmin,fmin,status,details]=glpkcc(f, A, B, lb, ub, ...
        CTYPE, vartype, SENSE, param);
    
    switch status
        case {5, 2}
            how = 'ok';
            exitflag = 1;
        otherwise
%             xmin = zeros(length(f), 1);
            how = 'infeasible';
            exitflag = -1;
    end

elseif solver==3
    % YALMIP / XPress
    [xmin,fmin,how,exitflag]=yalmipMILP(f,A,B,Aeq,Beq,lb,ub,vartype,param,options,'xpress');

elseif solver==4
    % YALMIP / Mosek
    [xmin,fmin,how,exitflag]=yalmipMILP(f,A,B,Aeq,Beq,lb,ub,vartype,param,options,'mosek');
    
elseif solver==5
    % YALMIP / bintprog.m
    [xmin,fmin,how,exitflag]=yalmipMILP(f,A,B,Aeq,Beq,lb,ub,vartype,param,options,'bintprog');
    
elseif solver==6,
    % YALMIP / CPLEX8 (cplexmex)
    [xmin,fmin,how,exitflag]=yalmipMILP(f,A,B,Aeq,Beq,lb,ub,vartype,param,options,'cplex-milp-cplexint');

elseif solver==7,
    % CPLEX interfaced with CPLEXMEX (by Nicolo Giorgetti)

    % all constraints are default of the form A <= b
    ctype = repmat('L', size(A,1), 1);
    if ~isempty(Aeq),
        nc = size(A,1);
        A = [A; Aeq];
        B = [B; Beq];
        ctype = [ctype; repmat('E', size(Aeq, 1), 1)];
    end
    if isfield(options, 'usex0'),
        x0 = options.usex0;
    else
        x0 = [];
    end
    if isfield(options, 'save_prob'),
        % note from cplexmex help file:
        % The file name can not be specified and defaults to "cplexpb.lp".
        save = 1;
    else
        save = 0;
    end

    sense = 1; % minimization
    [xmin,fmin,status,details]=cplexmex(sense, [], f, A, B, ctype, lb, ub, vartype, x0, param, save);
    switch status
        case {1,101,102}
            exitflag = 1;
            how = 'ok';
        case {3,103}
            % infeasible
            exitflag = -1;
            how = 'infeasible';
        case {2,118}
            % unbounded
            exitflag = -1;
            how = 'unbounded';
        case {4,119}
            % infeasible or unbounded
            exitflag = -1;
            how = 'InfOrUnb';
        otherwise
            exitflag = -1;
            how = 'other error';
    end

elseif solver==-1,
    % no solver available
    error('mpt_solveMILP: no MILP solver available!');
    
else
    error('mpt_solveMILP: unknown solver');
end


%------------------------------------------------------------------------------
function [xmin,fmin,how,exitflag]=yalmipMILP(f,A,B,Aeq,Beq,lb,ub,vartype,param,options,solver)

global mptOptions

if isfield(mptOptions, 'yalmipdata'),
    % use fast problem setup by exploiting a dummy problem structure stored in
    % mptOptions.yalmipdata

    % use hashtable indexing with a string
    
    % NOTE! we observed problems with certain MILPs being solved with GLPK. It
    % is therefore possible to switch to some other solver in branch&bound code
    % by setting modelprefix = 'miqp:'
    % NOTE! this change will require you to have some QP solver installed!
    
    modelprefix = 'milp:';
    
    model = mptOptions.yalmipdata([modelprefix solver]);
    
    if isempty(model),
        % this should never occur, since we prepare a model structure for every
        % solver in mpt_init. but in case we forget to add some solver in
        % mpt_init, we compute the dummy model once more here
        model = sub_getyalmipdata(solver, 'milp');
    end

    if isempty(model),
        % throw an error if model is still empty
        error(sprintf('Solver "%s" not available.', solver));
    end
    
    % use initial guess of optimizer if available
    if isfield(options, 'usex0'),
        x0 = options.usex0;
    else
        x0 = [];
    end
    
    % update the model structure with current data
    model = sub_setyalmipdata(model,[],f,A,B,Aeq,Beq,lb,ub,vartype,x0);
    
    % call an appropriate solver
    solution = eval([model.interfacedata.solver.call '(model.interfacedata);']);

    xmin = solution.Primal;
    fmin = f'*xmin;

    % analyze the solution
    if solution.problem==0,
        % solution is optimal, obtain optimizer and compute objective value
        how = 'ok';
        exitflag = 1;
    else
        % a problem occured
        if solution.problem < 0,
            how = 'nosolver';
            exitflag = -3;
        else
            how = 'infeasible';
            exitflag = -1;
        end
    end
    return
end

nx = size(A,2);

f = f;
x = sdpvar(nx,1);
F = set(A*x <= B);
if ~isempty(Aeq),
    F = F + set(Aeq*x == Beq);
end
if ~isempty(vartype),
    indbin = find(vartype=='B');
    if ~isempty(indbin),
        F = F + set(binary(x(indbin)));
    end
    
    indint = find(vartype=='I');
    if ~isempty(indint),
        F = F + set(integer(x(indint)));
    end
end
if ~isempty(lb),
    F = F + set(lb <= x);
end
if ~isempty(ub)
    F = F + set(x <= ub);
end

if ~isempty(options),
    optfields = fields(options);
    if length(optfields)==1,
        % if only "usex0" is present in "options", we can use the global
        % "sdpsettings" structure stored in mptOptions, just update this one
        % field. this will greatly reduce the setup time!
        if strcmpi(optfields{1}, 'usex0'),
            x0 = options.usex0;
            options = mptOptions.sdpsettings;
            options.solver = solver;
            options.usex0 = 1;  % tell YALMIP to use x0
            assign(x, x0);
        else
            options=sdpsettings(options,'Verbose',0,'warning',0,'solver',solver);
        end
    else
        options=sdpsettings(options,'Verbose',0,'warning',0,'solver',solver);
    end
else
    options = mptOptions.sdpsettings;
    options.solver = solver;
end

%solution = solvesdp(F, f'*x, sdpsettings(options, 'bnb.solver','nag'));
solution = solvesdp(F, f'*x, options);
xmin = double(x);
fmin = f'*xmin;

if solution.problem==0,
    how = 'ok';
    exitflag = 1;
elseif solution.problem < 0
    how = 'nosolver';
    exitflag = -3;
else
    how = 'infeasible';
    exitflag = -1;
end
