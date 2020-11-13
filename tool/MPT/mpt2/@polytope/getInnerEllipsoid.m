function E = getInnerEllipsoid(Pset,E,Options)
%GETINNERELLIPSOID Computes the largest ellipsoid inscribed in a polytope
%
% E = getInnerEllipsoid(Pset,x0,E,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function computes he largest ellipsoid inscribed in a polytope. 
% It is also possible to pass an ellipsoid (x-x0) E (x - x0) <= rho and 
% to compute the maximum rho such that the scaled ellipsoid is still contained
% in the polytope.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------  
%   Pset    -   Polytope object constraining the ellipsoid
%   E,x0    -   Optional: input ellipsoid with center x0, i.e.
%                           (x-x0) E^(-1) (x - x0) <= rho 
%               given as an ELLIPSOID object
%               The function then computes the maximum rho such that the 
%               ellipsoid is still contained in Pset.
%   Options
%     .plotresult    - If problem is in 2D and flag is set to 1, the result 
%                      will be plotted (Default: 0).
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%  E       -   Maximum volume ellipsoid,  (x-x0) E (x - x0) <= 1, returned as
%              an ELLIPSOID object
%
% ---------------------------------------------------------------------------
% LITERATURE
% ---------------------------------------------------------------------------
%   S. Boyd and L. Vanderberghe, Convex Optimization
%
% see also MPT_PLOTELLIP, GETOUTTERELLIPSOID


% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
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

global mptOptions

error(nargchk(1,3,nargin));

if nargin < 3,
    Options = [];
end
if ~isfield(Options,'plotresult')
    Options.plotresult=0;
end

[H,K]=double(Pset);
nx=size(H,2);

if(nargin>1)
    %compute the scaled ellipsoid
    [x0, E] = parameters(E);
    E=(E+E')/2;
    iL=inv(E);
    aux=zeros(size(K));
    for i=1:length(K);
        Ai=H(i,:);
        aux(i)=(K(i)-Ai*x0)^2/(Ai*iL*Ai');
    end
    rho=min(aux);
    E=E/rho';
else
    m = size(H,1);
    yalmip('clear' ) ;

    B = sdpvar(nx,nx);
    d = sdpvar(nx,1);
    
    myprog = lmi;
    for i=1:m
        myprog = myprog + set('||B*H(i,:)''||+H(i,:)*d<K(i,:)') ;
    end
    
    options = mptOptions.sdpsettings;
    myprog=myprog+set('B>0');
    
    % Ellipsoid E = {Bu+d | ||u||_2 <=1}
    solution = solvesdp(myprog,-logdet(B),options);
    E   = double(B);
    x0  = double(d);
    E=inv(E)'*inv(E);
end

if nx==2 & Options.plotresult
    plot(Pset)
    hold on
    [xe,ye]=mpt_plotellip(E,x0);
end

E = inv(E);
E = 0.5*(E + E');
E = ellipsoid(x0, E);
