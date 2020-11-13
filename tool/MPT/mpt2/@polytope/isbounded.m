function status = isbounded(P,Options)
%ISBOUNDED Checks if a polytope is bounded
%
% status = isbounded(P,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = ISBOUNDED(P) returns TRUE (1) if P is a bounded polytope
%
% Check is performed using Minkowski theorem on polytopes.
%
% USAGE:
%   isbounded(P)
%   isbounded(P,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                - polytope
% Options.abs_tol  - absolute tolerance
% Options.lpsolver - LP solver to use
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status           - Logical statement
%
% see also ISFULLDIM
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Miroslav Baric, Automatic Control Laboratory, ETH Zurich,
%               baric@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich
%          baotic@control.ee.ethz.ch

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


global mptOptions;

if nargin<2
    Options=[];
end

if ~isstruct(mptOptions),
    mpt_error;
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;      % absolute tolerance
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver = mptOptions.lpsolver;  % LP solver to use
end

lenP = length(P.Array);
if lenP>1,
    for ii=1:lenP,
        Q = P.Array{ii};
        status = isbounded(Q,Options);
        if status==0,
            return
        end
    end
else
    [ChebyC,ChebyR] = chebyball(P,Options);
    
    % polytope must be in minimal representation
    %
    if ( ~isminrep(P) ),
        P = reduce(P,Options);
    end
    
    if ChebyR == Inf,   % don't trust this one too much
        status = 0;
        return;
    elseif ChebyR <= Options.abs_tol,
        status = 1;
        return;
    end;
    
    A = P.H;
    b = P.K;
    
    % check the rank of matrix of vector normals
    %
    [n,m] = size(A);
    
    if m >= n,     % we have more variables than constraints: unbounded
        status = 0;
        return;
    end;
    
    nullA = null(A');
    colsNullA = size(nullA,2);
    rankA = n - colsNullA;
    
    if rankA < m,  % linearly dependant constraints
        status = 0;
        return;
    end;
    
    % we're using the Minkowski's theorem on polytopes:
    % Given a_1, ..., a_m unit vectors, and x_1, ... x_m > 0, there
    % exists a polytope having a_1,...,a_m as facets and x_1, ... x_m as
    % facet areas iff:
    %
    %  a_1 x_1 + ... + a_m x_m = 0
    %
    % Hence, checking boundedness of a polytope boils down to feasibility
    % LP.
    
    Aineq = -nullA;
    bineq = -ones(n,1);
    f = zeros(1,colsNullA);
    [xopt,fval,lambda,exitflag,how] = mpt_solveLPi(f,Aineq,bineq,[],[],[],Options.lpsolver);
    
    if strcmpi(how, 'ok')    % problem is feasible (unconstrained case)
        status = 1;          % cannot actually happen (a_i are
                             % not co-planar and sum(a_i x_i) =
                             % 0 with x_i > 0)
                             %  => polyhedron is bounded
    else
        % problem is infeasible => polyhedron is unbounded
        status = 0;
    end
end
