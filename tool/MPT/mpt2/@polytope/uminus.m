function P=uminus(P)
%UMINUS Unary minus for a polytope
%
% R = uminus(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% Mirrors polytope around the origin
%
% USAGE:
%   R = -P
%   R=uminus(P)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P   - polytope
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R   - polytope mirrored around the origin
%
% see also UPLUS
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

if ~isa(P, 'polytope')
  error('UMINUS: Argument MUST be a polytope object');
end

if ~isempty(P.Array)
    U = polytope;
    for ii=1:length(P.Array),
        Q=P.Array{ii};
        if ~Q.minrep
            Q=reduce(Q);
        end
        Q.H = -Q.H;
        Q.xCheb = -Q.xCheb;
        Q.vertices=-Q.vertices;
        bbox = Q.bbox;
        if ~isempty(bbox),
            Q.bbox = [-bbox(:, 2) -bbox(:, 1)];
        end
        U = [U Q];
    end
    P = U;
    return
end

if ~P.minrep
    P=reduce(P);
end

P.H=-P.H;
P.xCheb=-P.xCheb;
P.vertices=-P.vertices;
bbox = P.bbox;
if ~isempty(bbox),
    P.bbox = [-bbox(:, 2), -bbox(:, 1)];
end
