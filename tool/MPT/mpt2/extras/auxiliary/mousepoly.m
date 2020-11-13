function P = mousepoly(box)
%MOUSEPOLY Allows to specify polytope by mouse-clicks
%
% P = mousepoly(box)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Click on the figure to pick up points in 2D. When a right-mouse button is
% pressed, convex hull of given points is calculated and returned as a
% polytope object
% 
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% box    - range of the axis (10 by default)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% P      - polytope object
%
% see also HULL

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(0,1,nargin));

if nargin<1,
    box=10;
end
P = polytope;
figure; axis([-box box -box box]);hold on; grid on
title('Left-click to select points. Right-click to exit editing mode...');
X=[];
while 1
    [x1,x2,button]   =   ginput(1);   %graphically enter one point
    if button==3,
        P=polytope(X);
        opt.newfigure=0;
        hold on
        plot(P,opt);
        axis([-box box -box box]);
        hold off
        return
    end
    X=[X; x1 x2];
    plot(x1,x2,'kx','LineWidth',3);
end
hold off