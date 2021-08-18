function sol = mpt_mpmilp(Matrices, Options)
%MPT_MPMILP Multi-Parametric Mixed-Integer LP solver
%
% sol = mpt_mpmilp(Matrices, yalmipOptions)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Solves a multi-parametric Mixed Integer LP (mpMILP):
%
%  min    H U + F x
%   U
%  s.t.   G U <= W + E x
%         bndA*x <= bndb
%
% where certain elements of the U vector are known to be either continuous, or
% binary, or they belong to a finite alphabet.
%
% As a solution we get n regions
%   sol.Pn(i) = {x : H x <= K}
% 
% with the optimizer
%   U = sol.Fi{i}*x + sol.Gi{i}
%
% and the corresponding cost function expression
%   V(x) = sol.Bi{i}*x + sol.Ci{i}
%  
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Matrices - a struct with all the parameters which are needed.
%            See description above for explanation.
%   Matrices.H, Matrices.G, Matrices.W, Matrices.E
%   Matrices.bndA=bndA - Limits on exploration space, i.e. bndA*x<=bndb
%   Matrices.bndb=bndb
%   Matrices.vartype   - Vector of {'C', 'B', 'A'} indicies.
%                          vartype(i)='C' denotes a continuous variable
%                          vartype(i)='B' denotes a binary variable (0/1)
%                          vartype(i)='A' denotes a variable which can take
%                                          values from a finite alphabet
%   Matrices.alphabet  - Matrices.alphabet{i} must contain a list of values
%                        which the variable U(i) can take. For instance:
%                          Matrices.vartype(2)  = 'A';
%                          Matrices.alphabet{2} = [-0.5 0 3.4];
%                        and U(2) can only take values -0.5, 0, or 3.4
%
% Options - additional options to pass to YALMIP:
%    .verbose       - level of verbosity
%    .mp.algorithm  - which enumeration algorithm to use
%                     (Options.mp_algorithm=3 uses different enumeration)
%    .mp.presolve   - perform pre-solving if this option is true (default)
%                 
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% "sol" structure with following fields:
%   Pn,Fi,Gi    - for region Pn(i).H*x <= Pn(i).K computed input is
%                 U=Fi{i}*x+Gi{i}   
%   Pfinal       - Defines the feasible state space partition (i.e. union of
%                  all regions) as Phard.H*x<=Phard.K
%   Bi, Ci      - cost associate to every region
%
% see also MPT_MPLP. SOLVEMP

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch
% (C) 2005 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%     loefberg@control.ee.ethz.ch

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

error(nargchk(1,2,nargin));

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end
if nargin<2,
    Options = [];
end

%===========================================================================
% copy all fields from Options to respective fields of yalmipOptions
Options = mpt_defaultOptions(Options, ...
    'verbose', mptOptions.verbose);
yalmipOptions = mptOptions.sdpsettings;
f = fields(Options);
for ii = 1:length(f),
    yalmipOptions = setfield(yalmipOptions, f{ii}, getfield(Options, f{ii}));
end


if ~isfield(Matrices, 'vartype'),
    error('mpt_mpmilp: Matrices.vartype must be given.');
end
nu = size(Matrices.G, 2);
nx = size(Matrices.E, 2);
if nu~=length(Matrices.vartype),
    error('mpt_mpmilp: Matrices.vartype has wrong dimension.');
end
Matrices.vartype = upper(Matrices.vartype);

U = sdpvar(nu, 1);
X = sdpvar(nx, 1);
F = set([]);

% G*U <= W+E*X
F = F + set(Matrices.G*U <= Matrices.W + Matrices.E*X);

% bndA*x <= bndb
if ~isempty(Matrices.bndA),
    F = F + set(Matrices.bndA*X <= Matrices.bndb);
end

% derive bounds on variables
[ll,uu] = find_variable_bounds([Matrices.G -Matrices.E],Matrices.W);
ll = ll(1:length(U));
uu = uu(1:length(U));
bin = find(Matrices.vartype=='B');
ll(bin) = max(ll(bin),0);
uu(bin) = min(uu(bin),1);
constrained = find(~any(isinf([ll uu]),2));
bounds(U(constrained),ll(constrained),uu(constrained));

% constraint given variables to be binary
bin_vars = find(Matrices.vartype=='B');
if ~isempty(bin_vars),
    F = F + set(binary(U(bin_vars)));
end

% constraint given variables to be from finite alphabet
alph_vars = find(Matrices.vartype=='A');
if ~isempty(alph_vars),
    if ~isfield(Matrices, 'alphabet'),
        error('mpt_mpmilp: Matrices.alphabet must be given.');
    elseif ~iscell(Matrices.alphabet),
        error('mpt_mpmilp: Matrices.alphabet must be a cell array.');
    end
end
for ii = 1:length(alph_vars),
    ivar = alph_vars(ii);
    if ~isempty(Matrices.alphabet{ivar}),
        F = F + set(ismember(U(ivar), Matrices.alphabet{ivar}));
    else
        error(sprintf('mpt_mpmilp: No alphabet defined in Matrices.alphabet{%d}.', ivar));
    end
end

% objective
obj = Matrices.H*U + Matrices.F*X;

% solve the mp-MILP
yalmipsol = solvemp(F, obj, yalmipOptions, X);

% remove overlaps
try
    sol = mpt_removeOverlaps(yalmipsol, Options);
    sol = rmfield(sol, 'sysStruct');
    sol = rmfield(sol, 'probStruct');
    sol = rmfield(sol, 'details');
    sol = rmfield(sol, 'dynamics');
catch
    sol = [];
end
