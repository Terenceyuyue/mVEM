function plotj(varargin)
%PLOTJ Plots value function associated to a given controller
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots the value function defined over regions of a polyhedral partition
%
% USAGE:
%   plotj(ctrl)
%   plotj(ctrl, gridpoints)
%   plotj(ctrl, Options)
%   plotj(ctrl, gridpoints, Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl          - Explicit controller (MPTCTRL object)
% gridpoints    - Number of grid points in one axis (for PWQ case)
% Options       - Auxiliary options
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
%
% see also MPT_PLOTJ, MPT_PLOTPWA, MPT_PLOTPWQ
%

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
%
% ---------------------------------------------------------------------------

error(nargchk(1,3,nargin));

try
    mpt_plotJ(varargin{:});
catch
    rethrow(lasterror);
end
