function boxPoints=bboxvertices(BoxMin, BoxMax, n)
%BBOXVERTICES Computes all vertices of a bounding box
%
% V = bboxvertices(BoxMin, BoxMax, n)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes all vertices of a bounding box defined by lower and upper points.
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% BoxMin   - lower point of the bounding box
% BoxMax   - upper point of the bounding box
% n        - dimension of the bounding box
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% V        - vertices of the bounding box, stored column-wise!!!
%
% see also BOUNDING_BOX

% ---------------------------------------------------------------------------
% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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
% ---------------------------------------------------------------------------

binaryOne = '1';
boxPoint = zeros(n, 2^n);
for j=1:2^n
    index=dec2bin(j-1,n);
    for k=1:n
        if(index(k)==binaryOne)
            boxPoints(k, j)=BoxMax(k);
        else
            boxPoints(k, j)=BoxMin(k);
        end
    end
end
