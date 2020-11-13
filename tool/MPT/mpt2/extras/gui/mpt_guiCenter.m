function mpt_guiCenter(hObject)
%Computes coordinates of a GUI object such that it is in the center

% $Revision: 2.1 $ $Date: 2005/02/27 19:03:05 $
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

resolution = get(0, 'ScreenSize');
res_width = resolution(3);
res_height = resolution(4);
set(hObject, 'Units', 'pixels');
position = get(hObject, 'Position');
height = position(4);
width = position(3);
xpos = (round(res_width/2) - round(width/2));
ypos = (round(res_height/2) - round(height/2));
xpos = xpos + resolution(1);
ypos = ypos + resolution(2);
set(hObject, 'Position', [max(0,xpos) max(0, ypos) position(3:4)]);