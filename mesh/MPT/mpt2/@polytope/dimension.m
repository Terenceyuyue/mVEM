function nx = dimension(P)
%DIMENSION Returns dimension of the given polytope
%
% nx = dimension(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns dimension of the polytope P
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P   - Polytope
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% nx  - dimension of the given polytope P (or polyarray)
%
% see also NCONSTR, LENGTH
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
%   error('DIMENSION: Argument MUST be a polytope object');
% end

% dimension is simply number of columns of the H matrix
if length(P.Array)>0,
    nx=size(P.Array{1}.H,2);
else
    nx = size(P.H,2);
end