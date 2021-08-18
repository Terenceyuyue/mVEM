function statbar = mpt_statusbar(statbar, progress, smin, smax)
%Status bar function
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% status bar function
%
% Internal function.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
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

if nargin == 0 & nargout == 0,
    o='ShowHiddenHandles';
    t=get(0,o);
    set(0,o,'on');
    f=findobj(get(0,'Children'),'flat','Tag','StatusBar');
    set(0,o,t);
    if ~isempty(f),
        delete(f);
    end
    return
end

if nargin == 1 & nargout == 1,
    statbar = statusbar(statbar);
    return
end

progress = min(1, progress);

if nargin<3,
    smin = 0;
    smax = 1;
elseif nargin<4,
    smax = 1;
end

progress = progress*(smax-smin)+smin;
statbar = statusbar(min(1,progress), statbar);
