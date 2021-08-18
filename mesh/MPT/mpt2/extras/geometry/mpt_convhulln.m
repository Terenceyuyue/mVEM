function varargout = mpt_convhulln(V)
%MPT_CONVHULL Internal helper for convhulln
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Since qhull has many problems under Matlab R14, this helper is used to
% call convhulln() with additional options. However under R13 options
% cannot be passed to convhulln()...
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% V      - set of points
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% K      - indices of points forming the convex hull
% V      - volume of the convex hull (optional)
%

% Copyright is with the following author(s):
%
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich
%     kvasnica@control.ee.ethz.ch

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

% additional options for Matlab R14's version of qhull
%
% notice that following options hold only if the input should be fully
% dimensional
options = {'QJ', 'Qs', 'QR0'};

varargout = cell(1, max(nargout, 1));

if mpt_matlabrelease >= 14,
    % matlab R14 should use additional options
    [varargout{:}] = convhulln(V, options);

else
    % however matlab R13 does not support callinf delaunayn() with two
    % input arguments
    [varargout{:}] = convhulln(V);

end
