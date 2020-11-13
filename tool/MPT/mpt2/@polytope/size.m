function [nx,ny] = size(P,a)
%SIZE Returns size of the given polytope object
%
% [nx,ny]=size(P,a)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% [NX,NY] = SIZE(P)
%
% If P is a single polytope, NX and NY will be both 1
%
% If P is an array of polytopes, NY contains number of element in the
%   array. NX will be always 1.
%
% NOTE: In order to access number of constraints and dimension of P,
%       please use NCONSTR and DIMENSION respectively.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P      - polytope (or a polytope array)
% a      - index
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% nx     - is 1 always
% ny     - number of elements in P (i.e. length(P))
%
% see also LENGTH, DIMENSION, NCONSTR
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
%   error('SIZE: Argument MUST be a polytope object');
% end

if ~isempty(P.Array)
    ny = length(P.Array);
    nx = 1;
else
    nx = 1;
    ny = 1;
end
if nargin<2,
    a=1;
end
if a==2,
    nx=ny;
end
if nargout==1,
    nx=[nx ny];
end