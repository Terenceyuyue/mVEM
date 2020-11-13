function h_all=mpt_plotPWQ(Pn,lyapunovQ,lyapunovL,lyapunovC,meshgridpoints,Options);
%MPT_PLOTPWQ Plots a PWQ function defined over polyhedral partition
%
% handle = mpt_plotPWQ(Pn,Q,L,C,meshgridpoints,Options);
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots a PWQ function defined over a given polyhedral partition.
%
% PWQ = x'Qx + Lx + C.
%
% This function could be used either to plot PWQ Lyapunov function
% obtained by mpt_getPWQLyapFct or to plot the piece-wise quadratic
% cost function of an explicit solution
%
% USAGE:
%   mpt_plotPWQ(Pn,Ai,Bi,Ci)           to plot PWQ function over Pn
%
% NOTE:
%   If Pn contains overlapping regions, this function will only plot the lowest
%   value associated to given points.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn                - Polyhedral partition given by the polytope array Pn
% Q,L,C             - Cells containing parameters of the PWQ function x'Qx+x'L+C
% meshgridpoints    - number of grid points in one axis,
%                     (default: 30)
% Options.drawnow   - Whether or not to use the drawnow command to refresh the
%                     current plot (default is true)
% Options.shade     - Level of transparency (0 = fully transparent, 1 = solid)
% Options.edgecolor - specifies the color of edges. Default: 'k'.
% Options.edgewidth - specifies the width of edges. Default: 0.5.
% Options.lpsolver  - Solver for LPs when (and if) computing bounding box of Pn,
%                     (default: mptOptions.lpsolver)
% Options.newfigure - If set to 1, opens a new figure window,
%                     (default: mptOptions.newfigure)
% Options.showPn    - If set to 1, plots on polyhedral sets Pn,
%                     (default: 1)
% Options.samecolors - If set to 1, pieces of PWA function will be plotted
%                      in the same colors as the partition below
%                      (default: 0)
% Options.min_x1    - Rectangular space where PWQ function is computed
%                     (default: bounding box on Pn)
% Options.max_x1    - Rectangular space where PWQ function is computed
%                     (default: bounding box on Pn)
% Options.min_x2    - Rectangular space where PWQ function is computed
%                     (default: bounding box on Pn)
% Options.max_x2    - Rectangular space where PWQ function is computed
%                     (default: bounding box on Pn)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% handle.PWQ        - handle of the PWQ function
% handle.Pn         - handle of the plotted partition
%
% see also MPT_PLOTPWA, MPT_PLOTJ
%

% Copyright is with the following author(s):
%
% (c) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (c) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch
% (c) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
% (c) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (c) 2003 Marco Luethi, Automatic Control Laboratory, ETH Zurich,
%          mluethi@ee.ethz.ch

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

error(nargchk(4,6,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if(nargin<5 | isempty(meshgridpoints))
    meshgridpoints=30;
end

if isstruct(meshgridpoints) & nargin < 6,
    Options = meshgridpoints;
    meshgridpoints = 30;
elseif nargin<6
    Options=[];
end

Options = mpt_defaultOptions(Options, ...
    'lpsolver', mptOptions.lpsolver, ...
    'newfigure', mptOptions.newfigure, ...
    'showPn', 1, ...
    'samecolors', 0, ...
    'verbose', mptOptions.verbose, ...
    'drawnow', 1 );

if dimension(Pn)>2,
    error('mpt_plotPWQ: only 2D partitions supported!');
end

oneDimCase = (dimension(Pn)==1);

pwq_or_q=1;

if(isempty(lyapunovL) | isempty(lyapunovC))
     pwq_or_q=0;
     if(~iscell(lyapunovQ))
         Qiscell=0;
     else
         Qiscell=1;
     end
     tQ=lyapunovQ;
       
     clear lyapunovQ
     for i=1:length(Pn)
         if(Qiscell)
             lyapunovQ=tQ;
         else
             lyapunovQ{i}=tQ;    
         end
         lyapunovL{i}=[];
         lyapunovC{i}=[];
     end
end

if Options.newfigure,
   figure;  % open new figure window
   Options.newfigure = 0;  % tell any subsequent function not to open a yet another figure
elseif ~ishold,
    % else reuse current plot or open new figure if no plot exists
    newplot;
end

if Options.verbose > 0,
    disp('Creating PWQ plot...');
end

if Options.samecolors,
    maxlen = length(Pn);
    auxcolors=hsv(maxlen);
    multiplier=7;
    if mod(size(auxcolors,1),multiplier)==0,
        multiplier=multiplier+1;
    end
    for i=1:maxlen,
        jj=mod(i*multiplier,size(auxcolors,1))+1; % prepare enough colors for all polytopes, cycle through a given color map
        colors(i,:)=auxcolors(jj,:);
    end
end

% Use predefined range for x, or find the extremal values of Pn

if oneDimCase,
    V = [];
    for ii=1:length(Pn),
        V = [V; extreme(Pn(ii))];
    end
    min_x1 = min(V);
    max_x1 = max(V);
    %   produce a meshgrid between the extremal values
    x1=linspace(min_x1,max_x1,meshgridpoints);
    x2=linspace(0,0,meshgridpoints);
    [X,Y]=meshgrid(x1,x2);
    Z=repmat(Inf,meshgridpoints,meshgridpoints); %if meshgridpoint is in not in any region, set function value to Inf
    %(no plot for this point)
    C=zeros(meshgridpoints,meshgridpoints);
else
    V = [];
    [B, L, U] = bounding_box(Pn, struct('noPolyOutput', 1));
    if isfield(Options,'min_x1')
        min_x1=Options.min_x1;
    else
        min_x1 = L(1);
    end
    if isfield(Options,'max_x1')
        max_x1 = Options.max_x1;
    else
        max_x1 = U(1);
    end
    if isfield(Options,'min_x2')
        min_x2 = Options.min_x2;
    else
        min_x2 = L(2);
    end
    if isfield(Options,'max_x2')
        max_x2 = Options.max_x2;
    else
        max_x2 = U(2);
    end
    %   produce a meshgrid between the extremal values
    x1=linspace(min_x1,max_x1,meshgridpoints);
    x2=linspace(min_x2,max_x2,meshgridpoints);
    [X,Y]=meshgrid(x1,x2);
    Z=repmat(Inf,meshgridpoints,meshgridpoints); %if meshgridpoint is in not in any region, set function value to Inf
    %(no plot for this point)
    C=zeros(meshgridpoints,meshgridpoints);
end

for i=1:meshgridpoints
    for j=1:meshgridpoints
        if oneDimCase,
            x = x1(j);
        else
            x = [x1(j); x2(i)];
        end
        [isin, inwhich] = isinside(Pn, x);
        if ~isin,
            % point is not inside of any region, skip to next grid point
            continue
        else
            mincost = Inf;            
            for k=inwhich(:)',
                % now compute cost in each region, then take the minimal one
                % (because we can have overlaps
                if((pwq_or_q==1) & (~isempty(lyapunovL{k})) & (~isempty(lyapunovC{k})))
                    % calculate Lyapunov-function value
                    mincost = min(mincost, x'*lyapunovQ{k}*x+(lyapunovL{k}(:))'*x+lyapunovC{k}); 
                else
                    % calculate quadratic-function value
                    mincost = min(mincost, x'*lyapunovQ{k}*x); 
                end
            end
            zaux = mincost;
            if zaux<Z(i,j)
                Z(i,j)=zaux;
                C(i,j)=zaux;   % C is used for color (surf-plot)
            end
        end
    end
end

% plot feasible set: can not use ifa_region_plot,
% because then it is not possible that region and corresponding part of PWQ function have same color
legLab = [];  % legend labels
handles = []; % legend handles
xlabel('x_1','Fontsize',16); % LaTeX math symbols are directly supported!
ylabel('x_2','Fontsize',16);
h=gcf;
h1 = get(h,'CurrentAxes');
set(h1,'Fontname','times');
set(h1,'Fontsize',14);
title(sprintf('The PWQ function over %d regions', length(Pn)),'FontSize',18);
grid; 
%end of plot invariant-regions

plot_holded = ishold;

if oneDimCase,
    if ~plot_holded,
        % only hold the plot if user didn't do it manually        
        hold on
    end
    plot(X(1,:), Z(1,:), 'LineWidth', 3);
    ylabel('f_{PWQ}','Fontsize',16);
else
    h = surf(X,Y,Z,C);
    if isfield(Options, 'shade'),
        set(h, 'FaceAlpha', Options.shade);
    else
        set(gcf, 'Renderer', 'painters');
    end
    if isfield(Options, 'edgecolor')
        set(h, 'EdgeColor', Options.edgecolor);
    end    
    if isfield(Options, 'edgewidth')
        set(h, 'LineWidth', Options.edgewidth);
    end
end
handle = [];
if Options.showPn
    if ~plot_holded,
        % only hold the plot if user didn't do it manually
        hold on
    end
    plotOpt = Options;
    plotOpt.drawnow = 0;
    if Options.samecolors,
        plotOpt.color = hsv(length(Pn));
        handle=plot(Pn, plotOpt);
    else
        handle=plot(Pn, plotOpt);
    end
end

if Options.samecolors,
    colormap(hsv);
end
if ~oneDimCase,
    view(-37.5,20) 
end
if ~plot_holded,
    % don't do 'hold off' if user holded the plot manually
    hold off;
end
grid on
if Options.drawnow,
    drawnow;
end

h_all.PWQ = h;
h_all.Pn  = handle;

if nargout==0,
    clear h_all
end
