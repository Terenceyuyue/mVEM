function output = callgurobi(interfacedata)

% Author Johan Löfberg 
% $Id: callgurobimex.m,v 1.2 2010-03-23 11:11:46 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
interfacedata.gettime = 0;
n = length(c);

if ~isempty(ub)
    LB = lb;
    UB = ub;
    LB(binary_variables)  = round(LB(binary_variables));
    UB(binary_variables)  = round(UB(binary_variables));
    LB(integer_variables) = round(LB(integer_variables));
    UB(integer_variables) = round(UB(integer_variables));
end

if options.showprogress;showprogress('Calling GUROBI',options.showprogress);end

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [LB,UB,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,LB,UB);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(K.f+cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
end

SENSE = 1;     % Minimize
C = full(c);   % Must be full
B = full(F_struc(:,1));         % Must be full
A =-F_struc(:,2:end);

% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
CTYPE = [char(ones(K.f,1)*61); char(ones(K.l,1)*60)];
VARTYPE = char(ones(length(c),1)*67);
VARTYPE(setdiff(integer_variables,semicont_variables)) = 'I';
VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb
VARTYPE(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
VARTYPE(intersect(semicont_variables,integer_variables)) = 'N';


if options.savedebug
    save gurobidebug
end

if ~options.verbose
    options.gurobi.DisplayInterval = 1e12;
    options.gurobi.Display = 0;
end
% Call mex-interface
solvertime = clock; 
if isempty(binary_variables) & isempty(integer_variables)
    [x,val,flag,output,lambda] = gurobi_mex(C,SENSE,sparse(A),B,CTYPE,LB,UB,VARTYPE,options.gurobi);
else
    [x,val,flag,output] = gurobi_mex(C,SENSE,sparse(A),B,CTYPE,LB,UB,VARTYPE,options.gurobi);
    lambda = [];
end

if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

problem = 0;
D_struc = -lambda;

% Check, currently not exhaustive...
switch flag
    case {2}
        problem = 0;
    case {3}
        x = zeros(length(C),1);
        problem = 1;
    case 4
        x = zeros(length(C),1);
        problem = 12;
    case 5
        problem = 2;
    case {7,8,9}
        problem = 3;
    case {12}
        problem = 4;
    case {1,6,10,11}
        problem = 11;
        
    otherwise
        problem = -1;
end

% Save all data sent to solver?
if options.savesolverinput
	solverinput.A = A;
	solverinput.C = C;
	solverinput.B = B;
	solverinput.CTYPE = CTYPE;
	solverinput.LB = LB;
	solverinput.UB = UB;
    solverinput.param = options.glpk;
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput    
	solveroutput.x = x;
    solveroutput.val = val;
    solveroutput.flag = flag;
    solveroutput.lambda=lambda;
    solveroutput.output = output; 
else
	solveroutput = [];
end

% Standard interface 
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;