function drawdual2(cp,ce,pv,ev,varargin)
%DRAWDUAL2 draw a dual complex.
%   DRAWDUAL2(CP,CE,PV,EV) draws the dual complex created by
%   MAKEDUAL2.

%   Darren Engwirda : 2014 --
%   Email           : engwirda@mit.edu
%   Last updated    : 10/12/2014

    fc = [.95,.95,.55];
    ec = [.15,.15,.15];
%----------------------------------- extract optional inputs
    if (nargin >= +5), fc = varargin{1}; end
    if (nargin >= +6), ec = varargin{2}; end
    
%------------------------------------ triangulate dual cells
   [tp,tv] = triadual2(cp,ce,ev);
%----------------------------- expand cell colour onto trias
    if (size(fc,1) == size(tp,1))
    %-------------------------------- alloc. colour per tria
        tc = zeros(size(tv,1),1);
    %-------------------------------- expand colour per tria
        for ci = 1 : size(tp,1)
            for ti = tp(ci,1):tp(ci,2)
                tc(ti) = fc(ci) ;
            end
        end
    else
    %--------------------------------------- uniform colours
        tc = fc;
    end 
%----------------------------------- only draw indexed edges
    in = zeros(size(ev,1),1);
    for ci = 1 : size(cp,1)
        in(ce(cp(ci,2):cp(ci,3))) = +1;
    end
%------------------------------------------- save hold state
    hs = ishold; hold on;
%------------------------------------- draw dual edges/faces
    if (size(tc,1) ~= size(tv,1))
    patch('faces',tv,'vertices',pv,'facecolor',fc,'edgecolor','none');
    else
    patch('faces',tv,'vertices',pv,'facevertexcdata',tc,...
        'facecolor','flat','edgecolor','none');
    end
    patch('faces',ev(in==+1,:),'vertices',pv,'facecolor','none','edgecolor',ec);
%------------------------------------------- push hold state
    if(~hs),hold off;end
    
end

