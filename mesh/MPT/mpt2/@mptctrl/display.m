function sys = display(ctrl)
%DISPLAY Displays details about a given MPT controller
%
% sys = display(ctrl)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Displays details about a given MPT controller
% 
% This method is executed automatically if one types:
%   >> ctrl
% without the semicolumn (;) at the end
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl - MPT controller object
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sys  - String containing the information
%
% see also MPTCTRL
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

sys = char(ctrl);

fprintf('\n');
disp([inputname(1) '=']);
fprintf('\n');

disp(sys);

if isexplicit(ctrl)
    fprintf('\n Type ''struct(%s)'' to see all fields.\n', inputname(1))
    fprintf(' Fields can be accessed directly, e.g. type ''%s.Pn'' to get the controller partition.\n', inputname(1));
end
    
fprintf('\n');
return
