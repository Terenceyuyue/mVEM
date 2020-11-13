function value = get(P,name)
%GET Get polytope properties
%
% value = get(P,name)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Get polytope properties (field values) of the given polytope
%
% USAGE:
%   VALUE = GET(P,'PropertyName')
%   i.e. value = get(P,'H')
%
% List of properties that can be accessed:
%   H, K, minrep, normal, xCheb, RCheb, vertices
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P     - Polytopes
% name  - String giving name of the property
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% value - Value of the given property
%
% see also DOUBLE, SET
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

% if ~isa(P, 'polytope')
%   error('GET: First argument MUST be a polytope object');
% end

P=struct(P);

if ~isfield(P,name)
    error('Non-existing property (field)');
end
value=getfield(P,name);