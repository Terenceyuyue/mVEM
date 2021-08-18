function handle=plot(ctrl,Options)
%PLOT Plots regions of the explicit controller
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots the polyhedral parition obtained as a solution of the 
% optimal control problem obtained by mpt_control.
%
% NOTE: this function only supports explicit controllers!
%
% USAGE:
%   plot(ctrl)
%   plot(ctrl, Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl             - Explicit controller (MPTCTRL object)
% Options.noInfCol - Supress green-shade coloring for time-optimal solutions
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% handle       - handle to the graph
%
% see also MPT_PLOTTRAJECTORY, MPT_PLOTPWA, MPT_PLOTPWQ, MPT_CONTROL
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

error(nargchk(1,2,nargin));

if nargin<2,
    Options = [];
end
Options = mpt_defaultOptions(Options, ...
    'sameUcolors', 1 );

if ~isexplicit(ctrl)
    error('MPTCTRL/PLOT: first input must be an explicit controller!');
end

handle = mpt_plotPartition(struct(ctrl), Options);

if nargout==0
    clear handle
end
