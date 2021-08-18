function handle=mpt_plotPWA(PA,Fi,Gi,Options)
%MPT_PLOTPWA Plots a PWA function defined over a given polyhedral partition
%
% handle=mpt_plotPWA(Pn,L,C,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots a PWA function (L*x + C) defined over a given polyhedral partition.
%
% This function could be used either to plot value of control moves (Fi,Gi) or
% to plot the piece-wise linear cost over polyhedral partition (Ai,Bi)
%
% USAGE:
%   mpt_plotPWA(Pn,Fi,Gi)   to plot value of control moves
%   mpt_plotPWA(Pn,Bi,Ci)   to plot linear cost index
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn                      - Polyhedral partition of the state-space
% L,C                     - cell arrays containing a PWA function
% Options.shade           - Level of transparency (0 = fully transparent, 
%                           1 = solid). Default: 1.
% Options.drawnow         - Whether or not to use the drawnow command to refresh
%                           the current plot (default is true)
% Options.edgecolor       - specifies the color of edges. Default: 'k'.
% Options.edgewidth       - specifies the width of edges. Default: 0.5.
% Options.extreme_solver  - Which method to use for vertex enumeration 
%                           (see help extreme)
% Options.lpsolver        - LP solver to be used
% Options.abs_tol         - absolute tolerance
% Options.newfigure       - If set to 1, opens a new figure window
% Options.showPn          - If set to 1, plots on polyhedral sets Pn,
%                           (default: 1)
% Options.samecolors      - If set to 1, pieces of PWA function will be plotted
%                           in the same colors as the partition below
%                           (default: 0)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% handle  - handle of the plot object
%
% see also MPT_PLOTPWQ, MPT_PLOTU

% Copyright is with the following author(s):
%
% (c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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

error(nargchk(3,4,nargin));

if ~isa(PA,'polytope')
    error('mpt_plotPWA: First input argument MUST be a polytope');
end
if nargin<4,
    Options=[];
end

if ~iscell(Fi),
    error('mpt_plotPWA: Second argument must be a cell array!');
end
if ~iscell(Gi),
    error('mpt_plotPWA: Third argument must be a cell array!');
end

if length(PA)~=length(Fi) | length(PA)~=length(Gi)
    error('mpt_plotPWA: All input arguments must have same length!');
end

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

Options = mpt_defaultOptions(Options, ...
    'abs_tol', mptOptions.abs_tol, ...
    'axis', 'auto', ...
    'lpsolver', mptOptions.lpsolver, ...
    'extreme_solver', mptOptions.extreme_solver, ...
    'showPn', 0, ...
    'newfigure', mptOptions.newfigure, ...
    'samecolors', 0, ...
    'rowindex', 1, ...
    'clip_minmax', [], ...
    'drawnow', 1 );

% users can decide which line of Fi,Gi they want to consider
for i = 1:length(Fi),
    Fi{i} = Fi{i}(Options.rowindex, :);
    Gi{i} = Gi{i}(Options.rowindex, :);
end

index=0;

if Options.newfigure,
   figure;  % open new figure window
   Options.newfigure = 0;  % tell any subsequent function not to open a yet another figure
elseif ~ishold,
    % else reuse current plot or open new figure if no plot exists
    newplot;
end

handle=[];

maxlen=length(PA);

if Options.samecolors,
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

minu=Inf;
plotOpt = Options;
plotOpt.drawnow = 0;
if dimension(PA)==1,
    if Options.showPn,
        handle=plot(PA, plotOpt);
    end
    for ii=1:maxlen,
        [xc,rc]=chebyball(PA(ii));  % get chebyshev's radius
        if rc<=Options.abs_tol
            disp('mpt_plotPWA: Empty polytope detected!');
            continue
        end
        axis(Options.axis);
        [V,R,PA(ii)]=extreme(PA(ii),Options);    % compute extreme points and rays of polytope P
        V = sort(V);
        x1=V(1); x2=V(2);
        xv1 = x1*Fi{ii} + Gi{ii};
        xv2 = x2*Fi{ii} + Gi{ii};
        if ~isempty(Options.clip_minmax)
            xv1 = max(min(xv1, Options.clip_minmax(:, 2)), ...
                Options.clip_minmax(:, 1));
            xv2 = max(min(xv2, Options.clip_minmax(:, 2)), ...
                Options.clip_minmax(:, 1));
            h = line([x1;x2],[xv1;xv2],'LineWidth',4, 'LineStyle', '--', 'color', 'b');
        else
            h = line([x1;x2],[xv1;xv2],'LineWidth',5, 'LineStyle', '--', 'color', 'b');
%             h = line([x1;x2],[xv1;xv2],'LineWidth',3);
        end
        handle = [handle; h];
    end
    title(sprintf('PWA function over %d regions',length(PA)),'FontSize',18);
    xlabel('x','Fontsize',16); % LaTeX math symbols are directly supported!
    ylabel('f','Fontsize',16);
    grid;
    h=gcf;
    h1 = get(h,'CurrentAxes');
    set(h1,'Fontname','times');
    set(h1,'Fontsize',14);
    axis tight

    
elseif dimension(PA)==2,
    
    for ii=1:maxlen
        P=PA(ii);
        nb=nconstr(P);       % number of constraints
        dimP=dimension(P);   % dimension
        if dimP~=2,
            close
            error('mpt_plotPWA: Only 2D partitions can be plotted with this function!');
        end
        [xc,rc]=chebyball(P);  % get chebyshev's radius
        if rc<=Options.abs_tol
            disp('mpt_plotPWA: Empty polytope detected!');
            continue
        end
        axis(Options.axis);

        [V,R,PA(ii)]=extreme(P,Options);    % compute extreme points and rays of polytope P
        if size(R,1)>0
            error('mpt_plotPWA: Polytope is unbounded!'); % existence of rays means polytope is unbounded
        end

        % sort vertices in a cyclic way;
        x1=V(:,1);
        x2=V(:,2);

        ang=angle([(x1-xc(1))+(x2-xc(2))*sqrt(-1)]);
        [val,ind]=sort(ang);
        x1=x1(ind);
        x2=x2(ind);
        x3=[x1 x2]*Fi{ii}(1,:)'+Gi{ii}(1,:);   % the third dimension will be value of the control action, i.e. u=Fx+G
        if ~isempty(Options.clip_minmax)
            x3 = max(min(x3, Options.clip_minmax(:, 2)), ...
                Options.clip_minmax(:, 1));
        end
        if min(x3)<minu,
            minu=min(x3);
        end
        if Options.samecolors,
            h=patch(x1,x2,x3,colors(ii,:));
        else
            h=patch(x1,x2,x3,x3);
        end
        handle=[handle;h];
    end
    Jhandles = handle;
    if Options.showPn,
        % we are going to plot the polyhedral partition PA below
        for ii=1:maxlen,
            P=PA(ii);
            nb=nconstr(P);
            dimP=dimension(P);
            [xc,rc]=chebyball(P);
            if rc<=Options.abs_tol
                disp('mpt_plotPWA: Empty polytope detected!');
                continue
            end
            axis(Options.axis);

            [V,R,PA(ii)]=extreme(P,Options);
            if size(R,1)>0
                error('mpt_plotPWA: Polytope is unbounded!');
            end

            % sort vertices in a cyclic way;
            x1=V(:,1);
            x2=V(:,2);

            ang=angle([(x1-xc(1))+(x2-xc(2))*sqrt(-1)]);
            [val,ind]=sort(ang);
            x1=x1(ind);
            x2=x2(ind);
            x3=[x1 x2]*Fi{ii}(1,:)'+Gi{ii}(1,:);  % third dimension will be value of the control action, i.e. u=Fx+G
            if x3<minu,
                minu=x3;
            end
            if Options.samecolors,
                h=patch(x1,x2,minu*ones(size(x1)),colors(ii,:));
            else
                h=patch(x1,x2,minu*ones(size(x1)),'b');
            end
            handle=[handle;h];
        end
    end
    view(3);
    title(sprintf('PWA function over %d regions',length(PA)),'FontSize',18);
    xlabel('x_1','Fontsize',16); % LaTeX math symbols are directly supported!
    ylabel('x_2','Fontsize',16);
    zlabel('f','Fontsize',16);
    grid;
    h=gcf;
    if isfield(Options, 'shade')
        for ii = 1:length(Jhandles),
            set(Jhandles(ii), 'FaceAlpha', Options.shade);
        end
    else
        set(h, 'Renderer', 'painters');
    end
    if isfield(Options, 'edgecolor')
        for ii = 1:length(Jhandles),
            set(Jhandles(ii), 'EdgeColor', Options.edgecolor);
        end
    end
    if isfield(Options, 'edgewidth')
        for ii = 1:length(Jhandles),
            set(Jhandles(ii), 'LineWidth', Options.edgewidth);
        end
    end
    h1 = get(h,'CurrentAxes');
    set(h1,'Fontname','times');
    set(h1,'Fontsize',14);
    axis tight
end

if Options.drawnow,
    drawnow;
end
% If no outputs is asked, clear the variable.
if nargout == 0;
    clear('handle');
end
