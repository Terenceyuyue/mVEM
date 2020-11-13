function status = gt(P,Q,Options)
%GT Checks if polytope P is a strict superset of polytope Q
%
% status = gt(P,Q,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = GT(P,Q) returns TRUE (1) if P>Q (i.e. P is a strict superset of Q)
%
% USAGE:
%   P>Q
%   gt(P,Q)
%   gt(P,Q,Options)
%
% NOTE:
%   comparing two polyarrays involves making a minkowski sum of one of them
%   with an epsilon-box. This operation is computationally very expensive, so
%   try to use no-strict operator (>=) if possible.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P,Q                 - Polytopes
% Options.rel_tol     - relative tolerance
% Options.abs_tol     - absolute tolerance
% Options.lpsolver    - LP solver to use (see help mpt_solveLP)
% Options.verbose     - level of verbosity
% Options.elementwise - compares two polyarrays elementwise
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status           - Logical statement (true if P>Q, false otherwise)
%
% see also GE, LT, LE, EQ, NE
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
  error('GT: Argument MUST be a polytope object');
end

if nargin<3
    Options=[];
end

status = lt(Q,P,Options)';
