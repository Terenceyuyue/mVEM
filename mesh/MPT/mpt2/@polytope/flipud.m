function R = flipud(P)
%FLIPUD Flips array of polytopes
%
% R = flipud(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% If the input argument is a polytope array (i.e. (U Pi)),
% returns a polytope array R where R(1)=P(n), R(2)=P(n-1), ..., R(n)=P(1)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P   - Polytope array
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R   - Flipped polytope array
%
% see also FLIPLR, SUBSREF, LENGTH
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

% if ~isa(P,'polytope')
%     error('FLIPUD: Argument MUST be a polytope object');
% end

R = fliplr(P);