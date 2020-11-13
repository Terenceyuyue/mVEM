function handle=mpt_plotPartition(ctrl,Options)
%MPT_PLOTPARTITION Plots a polyhedral partition obtained by mpt_control
%
% handle=mpt_plotPartition(ctrl)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots the polyhedral parition obtained as a solution of the 
% optimal control problem obtained by mpt_control
%
% USAGE:
%   mpt_plotPartition(ctrl)
%   mpt_plotPartition(ctrl, Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl             - Explicit controller (MPTCTRL object)
% Options.noInfCol - Supress green-shade coloring for time-optimal solutions
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% handle       - handle to the graph
%
% see also MPT_PLOTTRAJECTORY, MPT_PLOTPWA, MPT_PLOTPWQ, MPT_CONTROL
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(1,2,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<1,
    error('mpt_plotPartition: wrong number of input arguments!');
end

if nargin<2,
    Options = [];
end
Options = mpt_defaultOptions(Options, ...
    'noInfCol', 0, ...
    'newfigure', mptOptions.newfigure, ...
    'sameUcolors', 0 );


if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl)
        error('This function supports only explicit controllers!');
    end
    ctrlStruct = struct(ctrl);
else
    ctrlStruct = ctrl;
    if ~mpt_isValidCS(ctrlStruct)
        error('mpt_plotPartition: First argument has to be a valid controller structure! See mpt_control for details.');
    end
end


% if regions do not overlap, we can enforce break in isinside at least it finds the first region associated
% to a given state (since there cannot be more than one if there are no overlaps)

Options.fastbreak = ~ctrlStruct.overlaps;

titlehandle=[];

details = ctrlStruct.details;
if (iscell(details) & isfield(details{1},'setdist')) | ...
        (iscell(ctrlStruct.sysStruct.A) & ctrlStruct.probStruct.subopt_lev==1),
    % iterative PWA solution
    lC = length(ctrlStruct.Ci);
    details = zeros(lC,1);
    for ii=1:lC
        details(lC-ii+1) = ctrlStruct.Ci{ii};
    end
elseif isfield(details,'regionHorizon'),
    details = details.regionHorizon;
elseif isfield(details,'IterStore') & ctrlStruct.overlaps,
    details = details.IterStore;
else
    details = ctrlStruct.details;
end

if Options.newfigure,
    figure;
    Options.newfigure = 0;  % we need to tell polytope/plot() not to open a yet another window
end

PA = fliplr(ctrlStruct.Pn);
nx = dimension(PA);
dimP = nx;

nxt = ctrlStruct.details.x0format.required - ...
    ctrlStruct.details.x0format.reference - ...
    ctrlStruct.details.x0format.uprev;
nu = ctrlStruct.details.dims.nu;
nref = ctrlStruct.details.x0format.reference;
isdumode = isfield(ctrlStruct.sysStruct, 'dumode') | ...
    ctrlStruct.probStruct.tracking==1 | ...
    ctrlStruct.details.x0format.uprev>0;
istracking = ctrlStruct.probStruct.tracking > 0;

if istracking | isdumode,
    Options.xsection = min([3 nxt])+1:dimP;
    mergetitles = 1;
    % to avoid 'section is empty' warnings
    Options.verbose = 0;
elseif nx > 3,
    Options.xsection = 4:nx;
    % do not add 'section through x3=... x4=...' to the figure title
    mergetitles = 1; 
    % to avoid 'section is empty' warnings
    Options.verbose = 0;
else
    mergetitles = 1;
end

if isfield(Options,'color'),
    [handle,titlehandle]=plot(PA,Options);
    
elseif Options.sameUcolors,
    [nx, nu] = mpt_sysStructInfo(ctrlStruct.sysStruct);
    ctable = sub_uniqueOpt(ctrlStruct.Fi, ctrlStruct.Gi, nu);
    colors = hsv(length(ctable.Fi));
    C = colors(ctable.Reg, :);
    [handle, titlehandle] = plot(ctrlStruct.Pn, struct('color', C));
    

elseif isempty(details),  
    % plot for constrained Finite-time optimal control problem
    [handle,titlehandle]=plot(PA,Options);

elseif isfield(ctrlStruct,'simplified') & ctrlStruct.simplified == 1
    % plot for simplified solutions (mpt_simplify)
    [handle,titlehandle]=plot(PA,Options);

else % for Infinite-time and iterative solution
    if isa(details,'double') & ~Options.noInfCol,
        regionhorizon=fliplr(details);
        auxcolors=hsv;           % take the color map
        maxcol = length(PA)+1;
        for i=1:maxcol,      % prepare enough colors
            jj=mod(i*2,size(auxcolors,1))+1;
            colors(i,:)=auxcolors(jj,:);
        end 
        Options.newfigure=0;     % do not open a new figure  (because a call to hsv opens the windows automatically
        infcolors=[];
        for i=1:length(PA)
            color=colors(regionhorizon(i)+1,:);
            infcolors=[infcolors; color];
        end
        Options.color=infcolors;
        [handle,titlehandle]=plot(PA,Options);
    else
        [handle,titlehandle]=plot(PA,Options);
    end
end

% set figure properties:
if ~isempty(titlehandle) & mergetitles
    oldtitle = get(titlehandle,'String');
else
    oldtitle='';
end
uprev = ctrlStruct.details.x0format.uprev;
if istracking
    for i = nxt+1:nxt+uprev,
        oldtitle = strrep(oldtitle, ...
            sprintf('x_{%d}', i), sprintf('u_{%d}', i-nxt));
    end
    for i = nxt+uprev+1:nx,
        oldtitle = strrep(oldtitle, ...
            sprintf('x_{%d}', i), sprintf('ref_{%d}', i-nxt-uprev));
    end
    oldtitle = sprintf('(tracking)\n%s', oldtitle);

elseif isdumode
    for i = nxt+1:nxt+uprev,
        oldtitle = strrep(oldtitle, ...
            sprintf('x_{%d}', i), sprintf('u_{%d}', i-nxt));
    end
    oldtitle = sprintf('(deltaU mode)\n%s', oldtitle);

end

title(['Controller partition with ' num2str(length(PA)) ...
    ' regions. ' oldtitle], 'FontSize', 14);

ylab = ''; zlab = '';
if isfield(ctrlStruct.sysStruct, 'StateName'),
    xlab = ctrlStruct.sysStruct.StateName{1};
else
    xlab = 'x_1';
end
if isfield(ctrlStruct.sysStruct, 'StateName') & nxt >= 2,
    ylab = ctrlStruct.sysStruct.StateName{2};
elseif nxt >= 2,
    ylab = 'x_2';
end
if isfield(ctrlStruct.sysStruct, 'StateName') & nxt >= 3,
    zlab = ctrlStruct.sysStruct.StateName{3};
elseif nxt >= 3,
    zlab = 'x_3';
end
if ~isempty(xlab),
    xlabel(xlab, 'FontSize', 14);
end
if ~isempty(ylab),
    ylabel(ylab, 'FontSize', 14);
end
if ~isempty(zlab),
    zlabel(zlab, 'FontSize', 14);
end
h=gcf;
h1 = get(h,'CurrentAxes');
set(h1,'Fontname','times');
set(h1,'Fontsize',14);

grid on;
axis tight

if nargout == 0;
    clear('handle');
end
