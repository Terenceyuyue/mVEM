function [xopt,fval,lambda,exitflag,how]=mpt_solveLPs(f,A,B,Aeq,Beq,x0,lpsolver)
%MPT_SOLVELPS Interface to various LP solvers ("safe" version)
%
% [xopt,fval,lambda,exitflag,how]=mpt_solveLPs(f,A,B,Aeq,Beq,x0,lpsolver);
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Solves an LP problem:
%
%     min  f'x
%     s.t. A*x<=B,
%          Aeq*x=Beq
%
% by using the method specified in 'lpsolver'
%
% WARNING: This is the "safe" version, which means that if the problem is
% infeasible with the default solver, another available solver will be used
% until: 
%   - a feasible solution is found, or
%   - all available solvers are tried
%
% WARNING: This is a "fast" version of mpt_solveLP:
%   - objective vector "f" must be a row vector !!!
%   - no (or at least very few) error checks are being performed
%   - all 7 input arguments must be specified
%
% NOTE! NOTE! NOTE! This function is used only for internal purposes, ordinary
% users should always use "mpt_solveLP" instead.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% f        - Optimization objective (must be a row vector!!!!!!!!)
% A,B      - Matrices defining inequality constraints
% Aeq,Beq  - Matrices defining equality constraints
% x0       - Initial value                         
% lpsolver - Which LP solver to use
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
% see also MPT_SOLVELP
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

[xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,A,B,Aeq,Beq,x0,lpsolver);
if exitflag ~= 1,
    % solution is not optimal, try another solver

    global mptOptions
    allsolvers = mptOptions.solvers.lp;

    for solverpos = 1:length(allsolvers),
        if allsolvers(solverpos) == lpsolver,
            % skip solver which was already used
            continue
        else
            % try next solver
            nextsolver = allsolvers(solverpos);
            [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(f,A,B,Aeq,Beq,x0,nextsolver);
            if exitflag == 1,
                % optimal solution found
                return
            end
            % otherwise continue with the loop which will try another solver
        end
    end
end
