function ind = uniquerows1(M, nc)
%UNIQUEROWS1 Fast method to obtain indices of unique rows of a given matrix
%
% ind = uniquerows1(M,nc)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns indices of unique rows of a given matrix
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% M   - matrix
% nc  - number of columns of M (i.e. size(M,2))
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% ind - indices of unique rows in M
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

% ---------------------------------------------------------------------------
% Legal note:
%          This library is free software; you can redistribute it and/or
%          modify it under the terms of the GNU Lesser General Public
%          License as published by the Free Software Foundation; either
%          version 2.1 of the License, or (at your option) any later version.
%
%          This library is distributed in the hope that it will be useful,
%          but WITHOUT ANY WARRANTY; without even the implied warranty of
%          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%          Lesser General Public License for more details.
% 
%          You should have received a copy of the GNU Lesser General Public
%          License along with this library; if not, write to the 
%          Free Software Foundation, Inc., 
%          59 Temple Place, Suite 330, 
%          Boston, MA  02111-1307  USA
%
% ---------------------------------------------------------------------------

rand1 = rand(nc,1);

hash1 = M*rand1;

%hash1 = round(hash1*1e12)/1e12;

fpos = numel(hash1);

[c,ind] = sort(hash1);

dc = diff(c);

d = (dc ~= 0);

d(fpos) = 1;

ind = ind(d);