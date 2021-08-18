function model = sub_setyalmipdata(model,H,f,A,B,Aeq,Beq,lb,ub,vartype,x0)

% Copyright is with the following author(s)
%
%(C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch

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

% equality constraints go first, i.e. A = [Aeq; A]; B = [Beq; B];

% model.interfacedata.F_struc = [b -A]
% model.interfacedata.K.f - number of equality constraints
% model.interfacedata.K.l - number of inequality constraints
% model.interfacedata.c   - linear part of the cost function
% model.interfacedata.Q   - quadratic term of the cost function (don't forget to multiply it by 2!!!)
% model.interfacedata.f   - constant part of the cost function
% model.interfacedata.x0  - initial guess
% model.interfacedata.lb  - lower bound
% model.interfacedata.ub  - upper bound
% model.interfacedata.binary_variables   - indices of binary variables
% model.interfacedata.integer_variables  - indices of integer variables
% model.interfacedata.monomtable         - identity matrix (with "nx" elements)

[nr, nx] = size(A);

interfacedata = model.interfacedata;

interfacedata.K.f = size(Aeq, 1);
interfacedata.K.l = nr;

if ~isempty(Aeq),
    interfacedata.problemclass.constraint.equalities.linear = 1;
end

A = [Aeq; A];
B = [Beq; B];

interfacedata.F_struc = sparse([B -A]);
interfacedata.c = f;
interfacedata.x0 = x0;
interfacedata.lb = lb;
interfacedata.ub = ub;
interfacedata.monomtable = eye(nx);

interfacedata.problemclass.constraint.binary = 0;
interfacedata.problemclass.constraint.integer = 0;
if ~isempty(vartype),
    interfacedata.binary_variables = find(vartype=='B')';
    interfacedata.integer_variables = find(vartype=='I')';
    if interfacedata.solver.constraint.binary==0 & ~isempty(interfacedata.binary_variables),
        % some solvers (e.g. MOSEK) only support integer variables. therefore we
        % have to denote binary variables as integer with 0 <= x <= 1 bounds
        interfacedata.integer_variables = union(interfacedata.binary_variables, ...
            interfacedata.integer_variables);
        if isempty(lb),
            lb = repmat(-Inf, size(A, 2), 1);
        end
        if isempty(ub),
            ub = repmat(Inf, size(A, 2), 1);
        end
        lb(interfacedata.binary_variables) = 0;
        ub(interfacedata.binary_variables) = 1;
        interfacedata.binary_variables = [];
    end
        
    if ~isempty(interfacedata.binary_variables),
        interfacedata.problemclass.constraint.binary = 1;
    end
    if ~isempty(interfacedata.integer_variables),
        interfacedata.problemclass.constraint.integer = 1;
    end
else
    interfacedata.binary_variables = [];
    interfacedata.integer_variables = [];
end

if isempty(H),
    % we have an LP / MILP
    interfacedata.Q = spalloc(nx, nx, 0);
    interfacedata.problemclass.objective.linear = 1;
    interfacedata.problemclass.objective.quadratic = struct('convex', 0, 'nonconvex', 0);
else
    % we have an QP / MIQP
    interfacedata.Q = H;
    interfacedata.problemclass.objective.linear = 0;
    interfacedata.problemclass.objective.quadratic = struct('convex', 1, 'nonconvex', 0);
end

% according to bnb.m this should be a speed hack
interfacedata.getsolvertime = 0;

model.interfacedata = interfacedata;
