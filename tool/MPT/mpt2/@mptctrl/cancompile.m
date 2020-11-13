function yesno = cancompile(ctrl)
%CANCOMPILE Returns true if a given controller can be compiled
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns true if a given controller can be compiled into an executable form to
% be used on target platforms.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl   - MPTCTRL object
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% yesno  - Boolean flag
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

yesno = 1;
if ~isexplicit(ctrl),
    % cannot compile on-line controllers
    yesno = 0;
    
elseif ctrl.overlaps & ctrl.probStruct.norm == 2 & ...
        ctrl.probStruct.subopt_lev == 0,
    % currently cannot compile optimal controllers for PWA systems with
    % quadratic cost
    yesno = 0;
    
elseif ctrl.probStruct.feedback ~= 0,
    % cannot export controllers based on feedback pre-stabilization
    yesno = 0;
    
end
