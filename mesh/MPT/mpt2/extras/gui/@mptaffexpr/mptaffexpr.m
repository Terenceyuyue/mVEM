function expr = mptaffexpr(nx, nu, gX, gU, gC)
%MPTAFFEXPR Constructor of MPTAFFEXPR objects

% Copyright is with the following author(s):
%
%(C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch

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

if nargin<3,
    gX = zeros(1, nx);
    gU = zeros(1, nu);
    gC = 0;
end
expr.gX = gX;
expr.gU = gU;
expr.gC = gC;

expr = class(expr, 'mptaffexpr');
