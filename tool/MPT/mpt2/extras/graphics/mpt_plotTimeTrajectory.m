function [X,U,Y,cost,trajectory]=mpt_plotTimeTrajectory(ctrl,x0,horizon,Options)
%MPT_PLOTTIMETRAJECTORY Plots trajectories of states, inputs, outputs and disturbances
%
% [X,U,Y,cost,trajectory]=mpt_plotTimeTrajectory(ctrl,x0,horizon,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% For a given initial state x0, computes and graphicale plots time evolution of
% states, inputs and outputs, as well as disturbances.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl              - An Explicit or an on-line controller
% x0                - initial state
% horizon           - for how many steps should the state evolution be computed
%                     If horizon=Inf, computes evolution of the state to origin
% Options.reference - If tracking is requested, provide the reference point
%                     in this variable (e.g. Options.reference = [5;0])
% Options.randdist  - If set to 1, randomly generated additive disturbance
%                       vector will be added to the state equation
% Options.openloop  - If 1, the open-loop solution will be computed
% Options.minnorm   - If closed-loop trajectory is computed, we stop the
%                       evolution if norm of a state decreases below this value
% Options.verbose   - Level of verbosity
% Options.bigfonts  - 1 for big fonts, 0 for standard fonts
% Options.newfigure - If set to 1, opens a new figure window
% Options.lpsolver  - Solver for LPs (see help mpt_solveLP)
% Options.abs_tol   - absolute tolerance
% Options.legend    - if set to true (default), legend will be displayed on each
%                     subplot. set this option to 0 to disable the legends.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% X, U, Y    - matrices containing evolution of states, control inputs and
%              system outputs
% cost       - contains cost from the given initial state to the origin
% trajectory - indices specifying in which region of PA the given state lies
%
% see also MPT_COMPUTETRAJECTORY, MPT_PLOTTRAJECTORY
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch

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

error(nargchk(2,4,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3 | isempty(horizon),
    horizon=Inf;
end
if nargin<4,
    Options=[];
    if ~isa(horizon,'double'),
        error('mpt_plotTimeTrajectory: Third argument (horizon) must be an integer! See ''help mpt_plotTimeTrajectory'' for details.');
    end
end

if ~isfield(Options,'lpsolver') % lpsolver to be used
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'abs_tol') % absolute tolerance
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'openloop') % compute the open-loop solution? (values 1/0, where 0 represents the closed-loop solution)
    Options.openloop=0;
end
if ~isfield(Options,'maxCtr') % maximum number of iterations in the closed-loop control-law evaluation
    Options.maxCtr=500;
end
if ~isfield(Options,'minnorm') % consider state reached the origin if norm of the state is less then Option.minnorm and abort trajectory update
    %Options.minnorm=0.05;
    % setting of this variable performed in mpt_computeTrajectory
end
if ~isfield(Options,'verbose') % level of verbosity
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'newfigure')
    Options.newfigure=mptOptions.newfigure;
end
if ~isfield(Options,'bigfonts'),
    Options.bigfonts=0;
end
if ~isfield(Options,'randdist') % if random disturbances should be added to the state vector
    Options.randdist=1;    
end
if horizon<1,
    error('mpt_plotTimeTrajectory: horizon MUST be a positive integer!');
end
if ~isfield(Options, 'guierrors')
    Options.guierrors = 0;
end
if ~isfield(Options, 'legend')
    Options.legend = 1;
end

if ~isinf(horizon),
    Options.minnorm=-1;
end

if isstruct(ctrl)
    inptype = 3;
elseif isa(ctrl, 'mptctrl') & isexplicit(ctrl)
    inptype = 1;
elseif isa(ctrl, 'mptctrl') & ~isexplicit(ctrl)
    inptype = 2;
else
    error('mpt_plotTimeTrajectory: unknown input type!');
end
ctrlStruct = ctrl;

if inptype==3 & ~mpt_isValidCS(ctrlStruct)
    error('mpt_plotTimeTrajectory: First argument has to be a valid controller structure! See mpt_control for details.');
elseif inptype==1,
    ctrlStruct = struct(ctrlStruct);
end

opt.verbose=0;

if horizon<=0,
    error('mpt_plotTimeTrajectory: Options.horizon MUST be a positive integer!');
end

if nargin<2 | isempty(x0),
    error('mpt_plotTimeTrajectory: initial state x0 MUST be given!');
end
[a,b]=size(x0);
if inptype==3 | inptype==1,
    if (a~=size(ctrlStruct.Fi{1},2) | b~=1) & ~ctrlStruct.probStruct.tracking & ~isfield(ctrlStruct.sysStruct,'dumode'),
        if Options.guierrors,
            error(['Initial state must be a ' num2str(size(ctrlStruct.Fi{1},2)) ' x 1 vector!']);
        else
            error(['mpt_plotTimeTrajectory: x0 must be a ' num2str(size(ctrlStruct.Fi{1},2)) ' x 1 vector!']);
        end
    end
end

if ctrlStruct.probStruct.tracking,
    if ~isfield(Options,'reference')
        if Options.guierrors,
            error('Reference point must be provided for tracking controllers!');
        else
            error('For tracking, reference point must be provided in Options.reference !');
        end
    end
    if isfield(ctrlStruct.sysStruct, 'dims'),
        nu = ctrlStruct.sysStruct.dims.nu;
        nxt = ctrlStruct.sysStruct.dims.nx;
        nyt = ctrlStruct.sysStruct.dims.ny;
    else
        if iscell(ctrlStruct.sysStruct.B),
            nu = size(ctrlStruct.sysStruct.B{1},2);
            if ctrlStruct.probStruct.tracking==1,
                nxt = (size(ctrlStruct.sysStruct.A{1},2)-nu)/2;
            else
                nxt = size(ctrlStruct.sysStruct.A{1},2)/2;
            end
        else
            nu = size(ctrlStruct.sysStruct.B,2);
            if ctrlStruct.probStruct.tracking==1,
                nxt = (size(ctrlStruct.sysStruct.A,2)-nu)/2;
            else
                nxt = size(ctrlStruct.sysStruct.A,2)/2;
            end
        end
    end
end

if Options.bigfonts,
  linewidth = 3;
else
  linewidth = 2;
end

% compute state, input, output and disturbance evolution from initial state x0
if inptype==3,
    % input was a controller structure
    [X,U,Y,D,cost,trajectory,feasible,dyns]=mpt_computeTrajectory(ctrlStruct,x0,horizon,Options);
else
    % input was an MPTCTRL
    [X,U,Y,D,cost,trajectory,feasible,dyns]=mpt_computeTrajectory(ctrl,x0,horizon,Options);
end
    
if ~feasible
    if Options.guierrors,
        error('No Feasible control law found!');
    else
        disp('mpt_plotTimeTrajectory: No feasible control law found!');
    end
    if isempty(U),
        %clear X U cost trajectory
        X = [];
        U = [];
        Y = [];
        cost = [];
        trajectory = [];
        return
    end
end


if Options.newfigure,
    figure;
else
    newplot
    hold off
end

if Options.openloop,
    T=0:size(X,1)-1;
    Y = [Y; Y(end,:)];
    U = [U; U(end,:)];
    D = [D; D(end,:)];
    dyns = [dyns(:); dyns(end)];
    endp=0;
else
    endp=1;
    T=0:size(X,1)-2;
    
    T=0:size(X,1)-1;
    Y = [Y; Y(end,:)];
    U = [U; U(end,:)];
    if isempty(D),
        D = 0;
    else
        D = [D; D(end,:)];
    end
    dyns = [dyns(:); dyns(end)];
    endp=0;
end

%% first subplot - state trajectories
subplot(2,2,1);
nx=size(X,2);
colors=repmat('brcmgyk',1,5);
if nx>length(colors),
    % prepare enough colors
    colors = repmat(colors,1,ceil((nx-length(colors))/8)+1);
end
handleX=[];
legentXtext={};
if 0 & ctrlStruct.probStruct.tracking,
    % de-activated because X and U are transofrmed to contain only proper "true"
    % states and inputs in mpt_computeTrajectory
    ctr = 0;
    Xmat = [];
    for ii=1:nxt,
        % plot states
        ctr = ctr+1;
        Xmat = [Xmat X(1:end-endp,ii)];
        legendXtext{ctr}=sprintf('x_{%d}',ii);
        if isfield(ctrlStruct.sysStruct, 'StateName'),
            % use user defined text labels if defined
            if ~isempty(ctrlStruct.sysStruct.StateName),
                legendXtext{ctr} = ctrlStruct.sysStruct.StateName{ii};
            end
        end
    end
    Xrefmat = [];
    for ii=1:nxt,
        % plot references
        ctr = ctr+1;
        Xrefmat = [Xrefmat X(1:end-endp,nxt+nu+ii)];
        refname = sprintf('r_{%d}',ii);
        if isfield(ctrlStruct.sysStruct, 'StateName'),
            if ~isempty(ctrlStruct.sysStruct.StateName),
                % use user defined text labels if defined
                refname = sprintf('Reference %s', ctrlStruct.sysStruct.StateName{ii});
            end
        end
        legendXtext{ctr} = refname;
    end
    handleX = plot(T,Xmat,T,Xrefmat,'--','LineWidth',linewidth);
    for ii=1:nxt,
        set(handleX(ii),'Color',colors(ii));
    end
    for ii=nxt+1:2*nxt,
        set(handleX(ii),'Color',colors(ii-nxt));
    end
else
    Xmat = [];
    for ii=1:nx,
        % create legend description (i.e. x_1, x_2, ... )
        Xmat = [Xmat X(1:end-endp,ii)];
        legendXtext{ii}=sprintf('x_{%d}',ii);
        if isfield(ctrlStruct.sysStruct, 'StateName'),
            if ~isempty(ctrlStruct.sysStruct.StateName),
                % use user defined text labels if defined
                legendXtext{ii} = ctrlStruct.sysStruct.StateName{ii};
            end
        end
    end
    handleX = plot(T,Xmat,'LineWidth',linewidth);
    for ii=1:length(handleX),
        set(handleX(ii),'Color',colors(ii));
    end
end
grid on
if Options.bigfonts,
    title('Evolution of states','FontSize',18);
    xlabel('Sampling Instances','Fontsize',16);
    ylabel('States','Fontsize',16);
    h=gcf;
    h1 = get(h,'CurrentAxes');
    set(h1,'Fontname','times');
    set(h1,'Fontsize',14);
else
    title('Evolution of states');
    xlabel('Sampling Instances');
    ylabel('States');
end
if Options.legend,
    legend(handleX,legendXtext);
end


%% second subplot - output trajectories
subplot(2,2,2);
%hold on
ny=size(Y,2);
legentYtext={};
Ymat = [];

if ctrlStruct.probStruct.tracking,
    ctr = 0;
    for ii=1:nyt,
        % plot outputs
        ctr = ctr+1;
        Ymat = [Ymat Y(:,ii)];
        legendYtext{ctr}=sprintf('y_{%d}',ii);
        if isfield(ctrlStruct.sysStruct, 'OutputName'),
            if ~isempty(ctrlStruct.sysStruct.OutputName),
                % use user defined text labels if defined
                legendYtext{ctr} = ctrlStruct.sysStruct.OutputName{ii};
            end
        end
    end
else
    for ii=1:ny,
        % create legend description (i.e. y_1, y_2, ...)
        Ymat = [Ymat Y(:,ii)];
        legendYtext{ii}=sprintf('y_{%d}',ii);
        if isfield(ctrlStruct.sysStruct, 'OutputName'),
            if ~isempty(ctrlStruct.sysStruct.OutputName),
                % use user defined text labels if defined
                legendYtext{ii} = ctrlStruct.sysStruct.OutputName{ii};
            end
        end
    end
    if isfield(ctrlStruct.probStruct,'yref'),
        Ymatref = [];
        for ii=1:ny,
            Ymatref = [Ymatref ctrlStruct.probStruct.yref(ii)*ones(length(T),1)];
            refname = sprintf('r_{%d}',ii);
            if isfield(ctrlStruct.sysStruct, 'OutputName'),
                if ~isempty(ctrlStruct.sysStruct.OutputName),
                    % use user defined text labels if defined
                    refname = sprintf('Reference %s', ctrlStruct.sysStruct.OutputName{ii});
                end
            end
            legendYtext{end+1}=refname;
        end
    end
end

% plot lines
if isfield(ctrlStruct.probStruct,'yref'),
    % plot references if given
    handleY = plot(T,Ymat,T,Ymatref,'--','LineWidth',linewidth);
    for ii=1:ny,
        set(handleY(ii),'Color',colors(ii));
    end
    for ii=ny+1:2*ny,
        set(handleY(ii),'Color',colors(ii-ny));
    end
else
    % plot outputs
    handleY = plot(T,Ymat,'LineWidth',linewidth);
    for ii=1:length(handleY),
        set(handleY(ii),'Color',colors(ii));
    end
end
    
grid on
if Options.bigfonts,
    title('Evolution of outputs','FontSize',18);
    xlabel('Sampling Instances','FontSize',16);
    ylabel('Outputs','FontSize',16);
    h=gcf;
    h1 = get(h,'CurrentAxes');
    set(h1,'Fontname','times');
    set(h1,'FontSize',14);
else
    title('Evolution of outputs');
    xlabel('Sampling Instances');
    ylabel('Outputs');
end    
if Options.legend,
    legend(handleY,legendYtext);
end


% third subplot - control moves
subplot(2,2,3);
%hold on
nu=size(U,2);
handleU=[];
legentUtext={};
Umat = [];

for ii=1:nu,
    % create legend description (i.e. u_1, u_2, ...)
    [Tx,Uy]=stairs(T,U(:,ii));
    Umat = [Umat Uy];
    legendUtext{ii}=sprintf('u_{%d}',ii);
    if isfield(ctrlStruct.sysStruct, 'InputName'),
        if ~isempty(ctrlStruct.sysStruct.InputName),
            % use user defined text labels if defined
            legendUtext{ii} = ctrlStruct.sysStruct.InputName{ii};
        end
    end
end
handleU=plot(Tx,Umat,'LineWidth',linewidth);
for ii=1:length(handleU),
    set(handleU(ii),'Color',colors(ii));
end
grid on
if Options.bigfonts,
    title('Evolution of control moves','FontSize',18);
    xlabel('Sampling Instances','Fontsize',16);
    ylabel('Inputs','Fontsize',16);
    h=gcf;
    h1 = get(h,'CurrentAxes');
    set(h1,'Fontname','times');
    set(h1,'Fontsize',14);
else
    title('Evolution of control moves');
    xlabel('Sampling Instances');
    ylabel('Inputs');
end
if Options.legend,
    legend(handleU,legendUtext);
end


if all(all(D==0)),
    % all disturbances are zero, plot dynamics instead
    % fourth subplot - dynamics
    [nx,nu,ny,nPWA]=mpt_sysStructInfo(ctrlStruct.sysStruct);
    subplot(2,2,4);
    if length(dyns)<length(T),
        dyns = dyns(1)*ones(1,length(T));
    end
    D = dyns;
        
    %hold on
    nd=size(D,2);
    handleD=[];
    legentDtext={};
    [Tx,Uy]=stairs(T,D);
    handleD = plot(Tx,Uy,'Color',colors(1),'LineWidth',linewidth);
    %axis([0 max(T) 0 nPWA+1]);
    grid on
    if Options.bigfonts,
        title('Active dynamics','FontSize',18);
        xlabel('Sampling Instances','Fontsize',16);
        ylabel('Dynamics','Fontsize',16);
        h=gcf;
        h1 = get(h,'CurrentAxes');
        set(h1,'Fontname','times');
        set(h1,'Fontsize',14);
    else
        title('Active dynamics');
        xlabel('Sampling Instances');
        ylabel('Dynamics');
    end

else

    % fourth subplot - disturbances
    subplot(2,2,4);
    %hold on
    nd=size(D,2);
    handleD=[];
    legentDtext={};
    Dmat = [];
    if ctrlStruct.probStruct.tracking,
        ctr = 0;
        for ii=1:nxt,
            % plot states
            ctr = ctr+1;
            Dmat = [Dmat D(:,ii)];
            legendDtext{ctr}=sprintf('d_{%d}',ii);
            if isfield(ctrlStruct.sysStruct, 'StateName'),
                if ~isempty(ctrlStruct.sysStruct.StateName),
                    % use user defined text labels if defined
                    legendDtext{ctr} = sprintf('Disturbance on %s',ctrlStruct.sysStruct.StateName{ii});
                end
            end
        end
    else
        for ii=1:nd,
            % create legend description (i.e. d_1, d_2, ...)
            Dmat = [Dmat D(:,ii)];
            legendDtext{ii}=sprintf('d_{%d}',ii);
            if isfield(ctrlStruct.sysStruct, 'StateName'),
                if ~isempty(ctrlStruct.sysStruct.StateName),
                    % use user defined text labels if defined
                    legendDtext{ii} = sprintf('Disturbance on %s',ctrlStruct.sysStruct.StateName{ii});
                end
            end
        end
    end
    handleD = plot(T,Dmat,'LineWidth',linewidth);
    for ii=1:length(handleD),
        set(handleD(ii),'Color',colors(ii));
    end
    grid on
    if Options.bigfonts,
        title('Evolution of disturbances','FontSize',18);
        xlabel('Sampling Instances','Fontsize',16);
        ylabel('Disturbances','Fontsize',16);
        h=gcf;
        h1 = get(h,'CurrentAxes');
        set(h1,'Fontname','times');
        set(h1,'Fontsize',14);
    else
        title('Evolution of disturbances');
        xlabel('Sampling Instances');
        ylabel('Disturbances');
    end
    if Options.legend,
        legend(handleD,legendDtext);
    end
end


if nargout==0,
    % if no output arguments requested, return
    clear X U cost trajectory
end



hold off
