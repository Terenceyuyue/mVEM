function handle=plotu(varargin)
%PLOT Plots value of the control action
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots value of the PWA control law with respect to a polyhedral partition
%
% NOTE: this function only supports explicit controllers!
%
% USAGE:
%   plotu(ctrl)
%   plotu(ctrl, uind)
%   plotu(ctrl, Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl                    - Explicit controller (MPTCTRL object)
% uind                    - Which output to plot (for MIMO systems)
% Options.extreme_solver  - Which method to use for vertex enumeration 
%                           (see 'help extreme')
% Options.lpsolver        - LP solver to be used
% Options.abs_tol         - absolute tolerance
% Options.newfigure       - If set to 1, opens a new figure window
% Options.verbose         - Level of verbosity
% Options.showPn          - If set to 1, plots on polyhedral sets Pn,
%                           (default: 0)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% handle       - handle to the graph
%
% see also MPT_PLOTU, MPT_PLOTPWA
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
    handle = mpt_plotU(varargin{:});
catch
    rethrow(lasterror);
end
if nargout < 1,
    clear handle
end
