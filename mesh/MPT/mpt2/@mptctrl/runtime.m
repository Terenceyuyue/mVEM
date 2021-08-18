function secs = runtime(ctrl)
%RUNTIME returns runtime needed to compute a given explicit controller
%
%   secs = runtime(ctrl)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns runtime in seconds which was needed to compute a given explicit
% controller
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl - explicit controller (an MPTCTRL object)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% secs - runtime in seconds. returns an empty vector if the controller is an
%        on-line controller.
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

if isfield(ctrl.details, 'runTime')
    secs = ctrl.details.runTime;
else
    secs = [];
end