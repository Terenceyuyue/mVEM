function [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,Aeq,Beq,x0,lpsolver,lb,ub)
%MPT_SOLVELP Interface to various LP solvers
%
% [xopt,fval,lambda,exitflag,how]=mpt_solveLP(f,A,B,Aeq,Beq,x0,lpsolver,lb,ub)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Solves an LP problem:
%
%     min  f'x
%     s.t.   A*x <= B
%          Aeq*x = Beq
%          lb <= x <= ub
%
% by using the method specified in 'solver'
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% f        - Optimization objective
% A,B      - Matrices defining inequality constraints
% Aeq,Beq  - Matrices defining equality constraints (optional)
% x0       - Initial value                          (optional)
% lpsolver - Which LP solver to use:
%              lpsolver=0: uses NAG (E04NAG.M)
%              lpsolver=9: uses NAG (E04MBF.M)
%              lpsolver=1: uses linprog.m
%              lpsolver=2: uses CPLEX9 (cplexint)
%              lpsolver=3: uses CDD Criss-Cross (cddmex)
%              lpsolver=4: uses GLPK (glpkmex)
%              lpsolver=5: uses CDD Dual Simplex (cddmex)
%              lpsolver=6: uses SeDuMi
%              lpsolver=7: uses QSopt (qsoptmex)
%              lpsolver=8: uses CPLEX8 (lp_cplex)
%              lpsolver=10: uses XPRESS
%              lpsolver=11: uses MOSEK
%              lpsolver=12: uses OOQP
%              lpsolver=13: uses CLP
%              lpsolver=14: uses BPMPD (bpmpd_mex)
%              lpsolver=15: uses CPLEX (cplexmex)
%              lpsolver=16: uses PDCO
%
% lb, ub   - lower and upper bounds on variables (optional)
%            Note! set "lpsolver=[]" to select the default solver.
%
% Note: if 'lpsolver' is not specified, mptOptions.lpsolver will be used instead
%       (see help mpt_init)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% xopt      - The optimizer
% fval      - Value of the objective
% lambda    - Vector of Lagrangian multipliers
% exitflag  - An integer specifying result of the optimization:
%                1 - feasible optimal solution found
%                0 - unbounded or undecided problem
%               -1 - infeasible problem
% how       - States the result of optimization ('ok', 'unbounded', 'infeasible')
%
% see also MPT_SOLVEQP, MPT_MPLP, MPT_MPQP

% Copyright is with the following author(s):
%
%(C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch
%(C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%         baotic@control.ee.ethz.ch

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
if ~isstruct(mptOptions)
    mpt_error;
end

if nargin==3,
    Aeq = [];
    Beq = [];
    x0 = [];
    lpsolver = mptOptions.lpsolver;
elseif nargin==4,
    Beq = 0;
    x0 = [];
    lpsolver = mptOptions.lpsolver;
elseif nargin==5
    x0 = [];
    lpsolver = mptOptions.lpsolver;
elseif nargin==6
    lpsolver = mptOptions.lpsolver;
elseif nargin<3,
    error('mpt_solveLP: Not enough input arguments!');
end

if isempty(lpsolver),
    lpsolver = mptOptions.lpsolver;
end

% handle lower and upper bounds on optimization variables
if nargin >= 8,
    if nargin == 8,
        error('mpt_solveLP: "ub" must be provided if "lb" is given.');
    end
    if ~isempty(lb) | ~isempty(ub),
        nx = max(size(A, 2), size(Aeq, 2));
        lb = lb(:);
        ub = ub(:);
        if length(lb) ~= length(ub),
            error('mpt_solveLP: incompatible dimensions of "lb" and "ub".');
        end
        if length(lb) ~= nx,
            error('mpt_solveLP: "lb" has wrong dimension.');
        end
        
        % include "lb" and "ub" into constraints
        A = [A; eye(nx); -eye(nx)];
        B = [B; ub; -lb];
    end
end
    
rescue = mptOptions.rescueLP;

% some solvers don't like +/- Inf terms in constraints
minfb = find(B == -Inf);
if ~isempty(minfb),
    % trivially infeasible problem
    xopt = repmat(NaN, length(f), 1);
    fval = NaN;
    lambda = [];
    exitflag = -1;
    how = 'infeasible';
    return
end
pinfb = find(B == Inf);
A(pinfb, :) = [];
B(pinfb) = [];

% NAG and CDD do not like inputs which are sparse matrices
no_sparse_solvers = [0 9 3 5];
yalmip_solvers = [6 10 11 12];

if ismember(lpsolver, no_sparse_solvers),
    % convert sparse matrices to full format for solvers which have problems
    % with sparse inputs
    if issparse(A),
        A = full(A);
    end
    if issparse(B),
        B = full(B);
    end
    if issparse(Aeq),
        Aeq = full(Aeq);
    end
    if issparse(Beq),
        Beq = full(Beq);
    end
    if issparse(f),
        f = full(f);
    end
end

% do all the computation in mpt_solveLPi
if mptOptions.rescueLP,
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPs(f(:)',A,B,Aeq,Beq,x0,lpsolver);
else
    [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f(:)',A,B,Aeq,Beq,x0,lpsolver);
end

if ~isempty(A),
    if max(A*xopt - B) > mptOptions.abs_tol,
        % CDD sometimes reports a feasible solution in spite of infeasible
        % constraints, catch this case
        how = 'infeasible';
        exitflag = -1;
    end
end
