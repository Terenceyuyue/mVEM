function status = isminrep(P)
%ISMINREP Checks if a given polytope is in minimal representation
%
% status = isminrep(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = ISMINREP(P) returns TRUE (1) if P is in minimal representation
%
% NOTE:
%   if P is a polyarray, STATUS will be a vector of statements
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P        - polytope
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status   - Logical statement
%
% see also ISNORMAL
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
%   error('ISMINREP: First argument MUST be a polytope object');
% end

lenP = length(P.Array);
if lenP>0,
    status=zeros(lenP,1);
    for ii=1:lenP,
        status(ii)= P.Array{ii}.minrep;
    end
else
    status= P.minrep;
end
