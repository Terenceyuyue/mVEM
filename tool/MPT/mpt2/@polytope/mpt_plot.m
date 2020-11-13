function [handle,titlehandle]=mpt_plot(varargin)
%PLOT Plots polytopes in 2D or 3D
%
% handle = plot(varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% PLOT(P,OPTIONS) plots polytope P if 1<=dimension(P)<=3
%
% Example
%   plot(P)             - plots P in default (red) color
%   plot(P,'g',Q,'b')   - plots P in green and Q in blue color
%   plot(PA,P)          - plots array of polytopes PA and a single polytope P
%   plot(PA,Options)    - plots PA using options given in the Options structure
%
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Polytopes
% Options.extreme_solver  - which method to use for extreme points computation
%                           (see help extreme)
% Options.newfigure   - if 1, each plot command will open a new figure window
% Options.color='r'   - sets default color
% Options.colormap    - sets a different colormap (default: 'hsv')
% Options.gradcolor   - if 1, third coordinate of a 3D polytope is used to
%                       color the polytopes with a gradient coloring (default 0)
% Options.edgecolor   - color of edges (default is black 'k')
% Options.wire=1/0    - plots polytopes in wireframe
% Options.wirestyle   - a string giving style of the wireframe, e.g. '--' or '.-'
% Options.wirecolor   - color of the wireframe (black is default)
% Options.linewidth   - width of the border of polytope (1 is default)
% Options.shade=0-1   - level of transparency (0 - transparent faces, 1 - opaque faces)
% Options.zvalue      - when plotting 2D polytopes and this parameter is given,
%                       it will shift all plotted polytopes along the 3rd
%                       coordinate by this value
% Options.elevate     - if set to 1, every 1D or 2D polytope is elevated by an
%                       increasing value of the z-axis. E.g. plot(P1, P2, P3)
%                       will plot P1 with z=1, P2 with z=2 and P3 with z=3.
%                       default is 0. this option takes precedence from
%                       Options.zvalue.
% Options.xdimension  - whether to plot a section in 2D or 3D (allowed values: 2,3)
%                       default value is 3 if this field is not present
% Options.xsection    - through which states to cut the plot (values e.g. [4 5] for a
%                       cut along x4 and x5). If you have a 3D polytope and want
%                       to intersect it with a plane x3=0, just set
%                       Options.xsection=3 and this will do the trick
% Options.xvalues     - at which values to cut (e.g. [0 10] will cut the plot at x4=0, x5=10)
% Options.verbose     - level of verbosity
% Options.marker      - Which marker to use for 1D plots (no by default)
% Options.linestyle   - Style of the border of each polytope (could be etiehr
%                       '-', ':', ';' or 'none') (default is '-')
% Options.persp_sect  - If set to 1 and Options.xsection is given (i.e. cut is
%                       made, uses perspective plotting (off by default)
% Options.drawnow     - If true (default), the "drawnow" command is used to
%                       refresh current figure after all polytopes are drawn
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% handle      - handle of the figure
%
% see also POLYTOPE
%

% Copyright is with the following author(s):
%
% (c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (c) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (c) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
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

nii=nargin;
Q=varargin{1};
if ~isa(Q,'polytope')
    error('PLOT: input argument MUST be a polytope');
end

if nii>1
    if isstruct(varargin{nii})
        Options=varargin{nii};
        nii=nii-1;
    else
        Options=[];
    end
else
    Options=[];
end

global mptOptions;
if ~isstruct(mptOptions)
    mpt_error;
end

Options = mpt_defaultOptions(Options, ...
    'abs_tol', mptOptions.abs_tol, ...
    'verbose', mptOptions.verbose, ...
    'axis', 'auto', ...
    'lpsolver', mptOptions.lpsolver, ...
    'extreme_solver', mptOptions.extreme_solver, ...
    'newfigure', mptOptions.newfigure, ...
    'marker', '', ...
    'gradcolor', 0, ...
    'linestyle', '-', ...
    'elevate', 0, ...
    'drawnow', 1, ...
    'wire', 0, ...
    'wirestyle', '-', ...
    'linewidth', 1, ...
    'fontsize', 10, ...
    'edgecolor', 'k', ...
    'persp_sect', 0 );

lenQ=0; dimP=0;
for ii=1:nargin
    if isa(varargin{ii},'polytope')
        dimP=dimension(varargin{ii});
        lenQ=lenQ+length(varargin{ii}); % determine how many polytopes are contained in all input arguments
    end
end

if ~isfield(Options, 'color')       % default color
    % to differentiate constraints cycle through the available colors
    %colors='ymcrgbk'; % even white is not Ok, it clash with the background
    
    %Color maps.
    %    hsv        - Hue-saturation-value color map.
    %    hot        - Black-red-yellow-white color map.
    %    gray       - Linear gray-scale color map.
    %    bone       - Gray-scale with tinge of blue color map.
    %    copper     - Linear copper-tone color map.
    %    pink       - Pastel shades of pink color map.
    %    white      - All white color map.
    %    flag       - Alternating red, white, blue, and black color map.
    %    lines      - Color map with the line colors.
    %    colorcube  - Enhanced color-cube color map.
    %    vga        - Windows colormap for 16 colors.
    %    jet        - Variant of HSV.
    %    prism      - Prism color map.
    %    cool       - Shades of cyan and magenta color map.
    %    autumn     - Shades of red and yellow color map.
    %    spring     - Shades of magenta and yellow color map.
    %    winter     - Shades of blue and green color map.
    %    summer     - Shades of green and yellow color map.
    
    if isfield(Options, 'colormap')
        Options.color = Options.colormap(1:lenQ, :);
    else
        auxcolors = hsv(lenQ);
        colors = auxcolors;
        multiplier=7;
        if mod(size(auxcolors,1),multiplier)==0
            multiplier=multiplier+1;
        end
        for i=1:lenQ
            jj=mod(i*multiplier,size(auxcolors,1))+1; % prepare enough colors for all polytopes, cycle through a given color map
            colors(i,:)=auxcolors(jj,:);
        end
        colors = flipud(colors);
        Options.color = colors;
    end
else
    if size(Options.color,1)~=length(varargin{1})
        if size(Options.color,1)>length(varargin{1})
            Options.color=Options.color(1:lenQ,:);
        else
            Options.color=repmat(Options.color(1, :),lenQ,1);
        end
    end
end

if ~isfield(Options, 'shade')       % transparency of plots
    if dimP<=2
        Options.shade=1;
    else
        Options.shade=0.5;
    end
end
if ~isfield(Options,'wirecolor')
    Options.wirecolor='k';
    iswirecolor=0;
else
    iswirecolor=1;
end

Indices=[];
Values=[];
titlehandle=[];

if size(Options.color,1)<lenQ
    error('plot: Not enough colors!');
end
if isfield(Options,'xdimension')
    plottabledimension=min(3,Options.xdimension);
else
    plottabledimension=3;
end

if isfield(Options,'xsection')
    if length(Options.xsection)<dimP-plottabledimension
        error(['Options.xsection MUST be a row vector with ' num2str(dimP-plottabledimension) ' elements!']);
    end
    Indices=Options.xsection;
    plottabledimension=min(plottabledimension,dimP-length(Indices));
else
    Indices=plottabledimension+1:dimP;
end
if isfield(Options,'xvalues')
    Values=Options.xvalues;
else
    Values=zeros(size(Indices));
end
if length(Values)~=length(Indices)
    error('Options.xsection and Options.xvalues MUST be row vectors of same dimension!');
end

index=0;

if Options.newfigure
    figure;  % open new figure window
else
    newplot; % get current figure (or create new figure)
end

handle=[];

%emptypoly=polytope;
emptypoly = mptOptions.emptypoly;
titlehandle=[];
cindex=0;
emptypolymsgctr = 0;
sectionemptyctr = 0;
problemtolctr = 0;

if dimP>2 && Options.verbose>0
    disp('Plotting...');
end

unbctr = 0; % counter for unbounded polytopes
extctr = 0; % counter for polytopes for which vertex enumeration failed

while index<nii
    index=index+1;
    Q=varargin{index};
    if ~isa(Q,'polytope')
        error('PLOT: input argument MUST be a polytope');
    end
    usercolor=0;
    if index<nii
        usercolor=0;
        if ischar(varargin{index+1})
            color=varargin{index+1};
            usercolor=1;
            index=index+1;
        else
            color=Options.color;
        end
    else
        color=Options.color;
    end

    maxlen=length(Q);

    for ii=1:maxlen
        cindex=cindex+1;
        if maxlen>1
            P = Q.Array{ii};
            if usercolor==0
                color = Options.color(cindex,:);
            end
        else
            if usercolor==0
                color = Options.color(cindex,:);
            end
            P = Q;
        end
        if ~P.minrep
            P=reduce(P);
        end
        [nb,dimP]=size(P.H);
        
        axislabels=1:dimP;
        
        % do the section:
        cutmade=0;
        indices=[];
        if dimP>plottabledimension
            indices=Indices;
            values=Values;
            if length(values)>length(indices)
                values=values(indices-min(indices)+1);
            elseif length(values)<length(indices)
                values=[values, zeros(1,length(indices)-length(values))];
            end
            values(indices>dimP)=[];
            indices(indices>dimP)=[];
            if length(indices)<dimP-plottabledimension
                indices=[indices setdiff(plottabledimension+1:dimP,indices)];
            end
            if length(values)<dimP-plottabledimension
                values=[values zeros(1,dimP-plottabledimension)];
            end
            Hb=P.H;
            Kb=P.K;
            Hb(:,indices)=[];
            Kb=P.K-P.H(:,indices)*values';
            try
                T=evalc('P=polytope(Hb,Kb);');    % to catch any warnings
            catch
                P=emptypoly;
            end
            axislabels = setdiff(1:dimP,indices);
            [nb,dimP]=size(P.H);
            cutmade=1;
        end
        if P.RCheb<=Options.abs_tol
            if isempty(indices)                      % no section made
                emptypolymsgctr = 1;
            else                                       % section is empty
                sectionemptyctr = 1;
            end
            % warning will be displayed only once at the end of a file
            continue
        end

        tempH.A=P.H;
        tempH.B=P.K;
        
        if dimP==1
            [tempV.V,tempV.R]=extreme(P,Options);
            x1=tempV.V(:,1);
            if Options.elevate==1
                x2 = cindex*ones(size(x1));
            else
                x2 = zeros(size(x1));
            end
            if isempty(Options.marker)
                h=line(x1,x2,'Color',color(1,:),'LineWidth',5);
            else
                h=line(x1,x2,'Color',color(1,:),'LineWidth',5,'Marker',Options.marker);
            end
            handle=[handle;h];
            
        elseif dimP==2
            shade = Options.shade;
            [tempV.V,tempV.R]=extreme(P,Options);
            if size(tempV.R,1)>0
                unbctr = unbctr + 1;
                continue
            end
            
            if isempty(tempV.V)
                continue
            end
            
            % sort vertices in a cyclic way;
            x1=tempV.V(:,1);
            x2=tempV.V(:,2);
            
            ang=angle((x1-P.xCheb(1))+(x2-P.xCheb(2))*sqrt(-1));
            [val,ind]=sort(ang);
            x1=x1(ind);
            x2=x2(ind);
            
            if Options.wire
                h=line([x1; x1(1)],[x2; x2(1)]);
                if iswirecolor==0
                    wcolor=color(1,:);
                else
                    wcolor=Options.wirecolor;
                end
                set(h,'Color',wcolor);
                set(h,'LineWidth',Options.linewidth);
                set(h,'LineStyle',Options.wirestyle);
            else
                if cutmade
                    if length(values)==1
                        if Options.persp_sect
                            if indices==1
                                h=patch(values*ones(size(x1)),x1,x2,color(1,:),'EdgeColor',Options.edgecolor,'FaceAlpha',shade,'LineStyle',Options.linestyle);
                            elseif indices==2
                                h=patch(x1,values*ones(size(x1)),x2,color(1,:),'EdgeColor',Options.edgecolor,'FaceAlpha',shade,'LineStyle',Options.linestyle);
                            else
                                h=patch(x1,x2,values*ones(size(x1)),color(1,:),'EdgeColor',Options.edgecolor,'FaceAlpha',shade,'LineStyle',Options.linestyle);
                            end
                        else
                            h=patch(x1,x2,values*ones(size(x1)),color(1,:),'EdgeColor',Options.edgecolor,'FaceAlpha',shade,'LineStyle',Options.linestyle);
                        end
                        zlabel(sprintf('x_{%d}',indices(1)));
                    else
                        h=patch(x1,x2,color(1,:),'EdgeColor',Options.edgecolor,'FaceAlpha',shade,'LineStyle',Options.linestyle);
                    end
                else
                    if isfield(Options, 'zvalue')
                        % move the 2D patch along the 3rd coordinate
                        zz = Options.zvalue*ones(size(x1));
                        h=patch(x1,x2,zz,color(1,:),'EdgeColor',Options.edgecolor,'FaceAlpha',shade,'LineStyle',Options.linestyle);
                    elseif Options.elevate==1
                        % move the 2D patch along the 3rd coordinate
                        zz = cindex*ones(size(x1));
                        h=patch(x1,x2,zz,color(1,:),'EdgeColor',Options.edgecolor,'FaceAlpha',shade,'LineStyle',Options.linestyle);
                    else
                        h=patch(x1,x2,color(1,:),'EdgeColor',Options.edgecolor,'FaceAlpha',shade,'LineStyle',Options.linestyle);
                    end
                end
            end
            if cutmade
                Title='Cut through ';
                for qq=1:length(indices)
                    Title=[Title sprintf('x_{%d}=%.2f',indices(qq),values(qq)) ' '];
                end
                titlehandle=title(Title);
            end
            handle=[handle;h];
        elseif dimP==3
            if Options.extreme_solver==3
                try
                    [tempV,VA]=cddmex('adj_extreme',tempH); % get vertices + adjacency list
                catch
                    % CDD failed, try analytical solution
                    if Options.verbose > 1
                        fprintf('CDD failed to compute extreme points, switching to analytical computation...\n');
                    end
                    opt = Options;
                    opt.extreme_solver = 0;
                    [tempV.V,tempV.R,QQ,VA,adjF]=extreme(P,opt);
                end
            else
                [tempV.V,tempV.R,QQ,VA,adjF]=extreme(P,Options);
            end
            if size(tempV.R,1)>0
                unbctr = unbctr + 1;
                continue
            end
            nVert=size(tempV.V,1);
            f=NaN*ones(nb,nVert);
            h=zeros(nb,1);
            
            for border=1:nb
                % find all points on this hyperplane
                points=find(abs(tempV.V*P.H(border,:)'-P.K(border))<=2*Options.abs_tol);
                npoints=length(points);
                if npoints<3
                    extctr = extctr + 1;
                    if npoints==0
                        continue
                    end
                end       
                % sort points in a cyclic way
                f(border,1)=points(npoints);
                active=f(border,1);
                previous=0;
                for kk=2:npoints
                    aux=intersect(VA{active},points);
                    if length(aux)~=2
                        if isbounded(P)
                            problemtolctr = 1;
                        else
                            unbctr = unbctr + 1;
                            continue
                        end
                    else
                        if aux(1)~=previous
                            f(border,kk)=aux(1);
                        else
                            f(border,kk)=aux(2);
                        end
                        previous=active;
                        active=f(border,kk);
                    end
                end
                if Options.wire
                    if any(isnan(f(border,1:npoints)))
                        % line can't handle NaNs, but patch can, strange...
                        continue
                    end
                    h(border)=line([tempV.V(f(border,1:npoints),1); tempV.V(f(border,1),1)],...
                        [tempV.V(f(border,1:npoints),2); tempV.V(f(border,1),2)],...
                        [tempV.V(f(border,1:npoints),3); tempV.V(f(border,1),3)]);
                    if iswirecolor==0
                        wcolor=color(1,:);
                    else
                        wcolor=Options.wirecolor;
                    end
                    set(h(border),'Color',wcolor);
                else
                    if Options.gradcolor
                        h(border)=patch(tempV.V(:,1),tempV.V(:,2),tempV.V(:,3),tempV.V(:,3),...
                            'Vertices',tempV.V,'Faces',f(border,1:npoints),'FaceAlpha',...
                            Options.shade,'EdgeColor',Options.edgecolor,'LineStyle',Options.linestyle);
                    else
                        h(border)=patch('Vertices',tempV.V,'Faces',f(border,1:npoints),...
                            'FaceVertexCData',tempV.V(:,3),'FaceColor',color,'FaceAlpha',...
                            Options.shade,'EdgeColor',Options.edgecolor,'LineStyle',Options.linestyle);
                    end
                end
                view(3);
            end
            if cutmade
                Title='Cut through ';
                for qq=1:length(indices)
                    Title=[Title sprintf('x_{%d}=%.2f',indices(qq),values(qq)) ' '];
                end
                titlehandle=title(Title);
            end
            handle=[handle;h];
        end
        if dimP>1
            if length(h)>1
                for ih=1:length(h)
                    if h(ih)~=0
                        set(h(ih),'LineWidth',Options.linewidth);
                    end
                end
            else
                set(h,'LineWidth',Options.linewidth);
            end
        end
    end
end
axis(Options.axis);
if ~exist('axislabels','var')
    axislabels = 1:dimP;
end
if dimP==1
    xlabel(sprintf('x_{%d}',axislabels(1)),'FontSize',Options.fontsize);
elseif dimP==2
    xlabel(sprintf('x_{%d}',axislabels(1)),'FontSize',Options.fontsize);
    ylabel(sprintf('x_{%d}',axislabels(2)),'FontSize',Options.fontsize);
    grid on
elseif dimP==3
    xlabel(sprintf('x_{%d}',axislabels(1)),'FontSize',Options.fontsize);
    ylabel(sprintf('x_{%d}',axislabels(2)),'FontSize',Options.fontsize);
    zlabel(sprintf('x_{%d}',axislabels(3)),'FontSize',Options.fontsize);
    grid on
end

h=gcf;
h1 = get(h,'CurrentAxes');
set(h1,'Fontname','times');
set(h1,'Fontsize',min(Options.fontsize,14));


if problemtolctr > 0 && Options.verbose > 0
    disp('PLOT: Problem detected. Try to change value of abs_tol in mpt_init. 1e-7 should be a good value.');
end
if emptypolymsgctr > 0 && Options.verbose > 0
    disp('PLOT: Empty polytope(s) detected!');          
end
if sectionemptyctr>0 && Options.verbose > 0
    disp('PLOT: Section through polytope(s) at the given location is empty! Try to plot the projection (help projection) or specify Options.xvalues.');
end
if unbctr > 0
    if unbctr == 1
        fprintf('PLOT: the polytope is unbounded\n');
    else
        fprintf('PLOT: %d polytopes are unbounded\n',unbctr);
    end
end
if extctr > 0
    fprintf('PLOT: Problem detected. Most probably extreme point enumeration failed for %d polytopes...\n',extctr);
    disp('Try to change value for extreme_solver.');
end
if Options.drawnow
    drawnow
end
% If no outputs is asked, clear the variable.
if nargout == 0
    clear('handle');
    clear('titlehandle');
end

