function status = ne(P,Q,Options)
%NE Checks if two polytopes are not equal
%
% status = ne(P,Q,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = NE(P,Q) returns TRUE (1) if P~=Q
%
% USAGE:
%   P~=Q
%   ne(P,Q)
%   ne(P,Q,Options)
%
% NOTE:
%  if P and/or Q are arrays of polytopes, STATUS is a matrix(vector) of statements
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P,Q              - Polytopes
% Options.rel_tol  - relative tolerance
% Options.abs_tol  - absolute tolerance
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status           - Logical statement
%
% see also EQ, LE, GE, LT, GT
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
  error('NE: Argument MUST be a polytope object');
end

global mptOptions;

if ~exist('mptOptions','var'),
    mpt_error;
end

if nargin<3
    Options=[];
end

if ~isfield(Options,'rel_tol')
    Options.rel_tol=mptOptions.rel_tol;    % relative tolerance
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end

status = ~eq(P,Q,Options);
