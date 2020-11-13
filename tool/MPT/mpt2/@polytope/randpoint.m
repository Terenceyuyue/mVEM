function x = randpoint(P)
%RANDPOINT Returns coordinates of a random point x \in P
%
% x = randpoint(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns coordinates of a random point x \in P
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P - polytope object
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% x - coordinates of the random point
%

% Copyright is with the following author(s):
%
% (C) 2008 Michal Kvasnica, Slovak University of Technology in Bratislava
%          michal.kvasnica@stuba.sk

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

if length(P) > 1
    error('Polytope arrays are not supported.');
end
if ~isfulldim(P)
    x = [];
    return
end

% compute extremal vertices of P
E = extreme(P);
ne = size(E, 1);

% create a random convex combination with sum(L)=1
L = rand(ne-1, 1);
L = L./(ne-1);
L = [L; 1-sum(L)];
L = L(randperm(ne));

% the random point is a random convex combination of the extremal vertices
x = E'*L;
