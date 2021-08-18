function status = isfulldim(P,Options)
%ISFULLDIM Checks if a polytope is full dimensional
%
% status = isfulldim(P,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = ISFULLDIM(P) returns TRUE (1) if P is a full dimensional polytope
%
% Note that full dimensionality of a polytope does not imply that it is also
% bounded!
%
% USAGE:
%   isfulldim(P)
%   isfulldim(P,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                - polytope
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
% see also ISBOUNDED
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


global mptOptions;

abs_tol = mptOptions.abs_tol;

lenP = length(P.Array);
if (lenP>0),
    status=zeros(lenP,1);
    for ii=1:lenP,
        if isempty(P.Array{ii}.RCheb)
            [P.Array{ii}.xCheb,P.Array{ii}.RCheb]=chebyball(P.Array{ii});
        end
        status(ii)= (P.Array{ii}.RCheb>=abs_tol); % polytope is fully dimensional if chebyshev's radius is positive
    end
    status = all(status);
else
    RCheb = P.RCheb;
    if isempty(RCheb)
        [xCheb,RCheb]=chebyball(P);
    end
    status= (RCheb>=abs_tol);    % polytope is fully dimensional if chebyshev's radius is positive
end
