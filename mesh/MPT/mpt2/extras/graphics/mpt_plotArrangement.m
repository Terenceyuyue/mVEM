function P=mpt_plotArrangement(H,K,BBox,Options)
%MPT_PLOTARRANGEMENT Plots hyperplane arrangement of a polytope in H-representation
%
% P=mpt_plotArrangement(H,K,BBox,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% P=mpt_plotArrangement(H,K,BBox,Options) plots hyperplane arrangement if
% 1<=dimension(P)<=3 ; Otherwise a projection onto 3D space is plotted. The
% initial hyperplanes are given as Hx<=K ; 
%
%   EXAMPLE:
%       P=mpt_plotArrangement([eye(3);-eye(3)],ones(6,1));
%
%       Options.randcolors=0;
%       P=mpt_plotArrangement([eye(2);-eye(2)],ones(4,1),unitbox(2,4),Options);
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% H,K           -   Hyperplanes Hx <= K
%                   (It is also possible to pass a polytope object H)
% BBox          -   Polytope object defining the outter box
% Options       -   used for plotting ("help plot")
% Options.randcolors    if set to 1 random colors will be used for the arrangement.
%                       Otherwise the initial polytope will be black and everything 
%                       else white.
%   
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% P             -   Polyarray containing all polytopes of the arrangement
%
% see also PLOT
%

% Copyright is with the following author(s):
%
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

error(nargchk(1,4,nargin));
if(isa(H,'polytope'))
    if(nargin>1)
       Options=K;
    end
    [H,K]=double(H);
end

binOne=dec2bin(1);
[m,n]=size(H);

if(nargin>2)
    [Hbox,Kbox]=double(BBox);
else
    Hbox=[eye(n);-eye(n)];
    Kbox=ones(2*n,1)*10;
end
if(nargin<4 | ~isfield(Options, 'randcolors'))
    Options.randcolors=1;
end
if(nargin<4 | ~isfield(Options, 'linewidth'))
    Options.linewidth=2;
end


P=polytope;
for i=1:2^m  %compute all possible intersections of hyperplanes
    binVec=dec2bin(i-1,m);
    ind=find(binVec==binOne);
    Ht=H; Kt=K;
    Ht(ind,:)=-Ht(ind,:);
    Kt(ind,:)=-Kt(ind,:);
    Pt=polytope([Ht; Hbox],[Kt; Kbox]);
    P=[P Pt];
end

if(Options.randcolors)
    Options.color=ones(length(P),3);
    Options.color(1,:)=[1 1 0.7];
    if(isfield(Options,'wire') & Options.wire==1)
        plot(P,'k',Options);
    else
        plot(P,Options);
    end
else
    plot(P);
end
axis tight

if nargout==0
    clear P
end