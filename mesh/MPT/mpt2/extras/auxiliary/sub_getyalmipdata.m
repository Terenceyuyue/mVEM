function model = sub_getyalmipdata(solver, solvertype)
% SUB_GETYALMIPDATA Internal function

% The purpose of this function is to obtain model of internal YALMIP structure
% which will be used to quickly initialize a given optimization problem without
% creating it using sdpvar / set / solvesdp

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

global mptOptions

if ~isstruct(mptOptions)
    mpt_error;
end

% set up a dummy mixed-integer problem
f = 1;
if ~isempty(findstr(solvertype, 'lp')),
    H = 0;
else
    H = 1;
end
A = 1;
B = 1;
Aeq = [];
Beq = [];
lb = [];
ub = [];
if ~isempty(findstr(solvertype, 'mi')),
    vartype = 'B';
else
    vartype = [];
end

[nc,nx] = size(A);
f = f(:);
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

options = mptOptions.sdpsettings;
options.solver = solver;
options.verbose = 0;
options.bnb.maxiter = 5000; % maximum number of iterations in branch&bound code

try
    % export the dummy model
    [dummy1, dummy2, dummy3, model.interfacedata] = export(F, 0.5*x'*H*x + f'*x, options);
catch
    model = [];
end
