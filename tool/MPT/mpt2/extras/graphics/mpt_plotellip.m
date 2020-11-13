function [xe,ye]=mpt_plotellip(A,b,c,linetype)
%mpt_PLOTELLIP function to plot polytopes
%
% [xe,ye]=mpt_plotellip(A,b,c,linetype)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% There are two ways to call the function:
% 
% 
% plotellip(A,xc)           plots ellipsoid (x-xc)*A*(x-xc) = 1 in R^2, or
% OR
% plotellip(A,b,c)          plots ellipsoid x'*A*x + b'*x + c = 0 in R^2
%
% linetype                  optional field used for plotting

% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%  xe,ye                    coordinates of points on ellipsoid   
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
%   S. Boyd and L. Vanderberghe, Convex Optimization
%
% see also MPT_GETINNERELLIPSOID, MPT_GETOUTTERELLIPSOID


% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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


if nargin==2
  type1=1;  xc=b;  linetype='-';
elseif nargin==3 & ischar(c)
  type1=1;  xc=b;  linetype=c;
elseif nargin==3
  type1=0;  linetype='-';
else
  type1=0;
end

A=.5*(A+A');
if min(eig(A))<=0
  disp('Invalid ellipsoid -- ''A'' not positive definite'), xe=[];, ye=[];, return
end
if ~type1
  xc=A\(-b./2);
  k=xc'*A*xc-c;
  if k<=0
    disp('Invalid ellipsoid'), xe=[];, ye=[];, return
  end
  A=A/k;
end

x0 = xc(1); y0 = xc(2);
R = chol(A);
theta = linspace(-pi,pi,100);
xy_tilde = [cos(theta); sin(theta)];
invR = inv(R);
xy_bar = invR*xy_tilde;
xe = xy_bar(1,:)+x0;
ye = xy_bar(2,:)+y0;
plot(xe,ye,linetype,'LineWidth',2);
hold on
plot(xc(1),xc(2),'b*','LineWidth',4);

axis equal
axis tight

if nargout==0
    clear xe ye
end
    
