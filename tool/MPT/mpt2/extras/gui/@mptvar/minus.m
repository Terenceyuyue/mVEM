function expr = minus(var1, var2)
%GE Difference of two MPTVAR objects

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
if ~(isa(var1, 'mptvar') | isa(var1, 'double'))
    error(sprintf('Variables of class "%s" not supported!', class(var1)));
end
if ~(isa(var2, 'mptvar') | isa(var2, 'double'))
    error(sprintf('Variables of class "%s" not supported!', class(var2)));
end

var2 = -var2;

expr = plus(var1, var2);
