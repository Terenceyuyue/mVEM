function status = ge(P,Q,Options)
%GE Checks if polytope P is a superset of polytope Q
%
% status = ge(P,Q,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = GE(P,Q) returns TRUE (1) if P>=Q (i.e. P is a superset of Q)
%
% USAGE:
%   P >= Q
%   ge(P,Q)
%   ge(P,Q,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P,Q              - Polytopes
% Options.rel_tol  - relative tolerance
% Options.abs_tol  - absolute tolerance
% Options.lpsolver - LP solver to use (see help mpt_solveLP)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status           - Logical statement (true if P>=Q, false otherwise)
%
% see also EQ, GT, LE, LT, NE
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
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

if ~isa(P, 'polytope') | ~isa(Q, 'polytope')
  error('GE: Argument MUST be a polytope object');
end

if nargin<3
    Options=[];
end

status = le(Q,P,Options)';