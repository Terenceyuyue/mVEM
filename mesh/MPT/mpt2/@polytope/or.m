function [Pu,how] = or(varargin)
%OR Convex union of n polytopes
%
% [Pu, how] = or(varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function tries to compute convex union of polytopes passed as input
% arguments delimited by comma
%
% USAGE:
%   U = P1 | P2 | P3
%   [U,how] = or(P1,P2,P3)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Polytopes
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Pu   - polytope describing convex union in case such one exists, return array of input arguments otherwise
% how = 1  if union is convex
%     = 0  union is not convex (in this case Pu a polytope array of the input arguments)
%
% see also UNION
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

global mptOptions

Pn = mptOptions.emptypoly;
for ii = 1:nargin,
    if ~isa(varargin{ii}, 'polytope')
        error('OR: argument MUST be a polytope object!');
    end
    Pn = [Pn varargin{ii}];
end
[Pu,how] = union(Pn);
