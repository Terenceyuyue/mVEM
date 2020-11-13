function [X,U,cost,trajectory]=mpt_plotTrajectory(ctrl,Options)
%MPT_PLOTTRAJECTORY Graphical interface for piloting trajectories for LTI and PWA systems subject to control
%
% [X,U,cost,trajectory]=mpt_plotTrajectory(ctrl)
% [X,U,cost,trajectory]=mpt_plotTrajectory(ctrl,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% This function allows the user to graphically enter an initial state in the
% state space. Then, for the given state, the state evolution trajectory is
% computed and displayed.
%
% To enter the point, click on an open plot of the state space.
%
% Note: Press right mouse button if you want to abort.
%
% If the system is subject to additive disturbance, random uncertainty
% is added to the state vector by default. To suppress this, use
% Options.randdist = 0;
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl              - Explicit controller (MPTCTRL object)
% Options.x0        - defines initial state (if provided, mouse interface
%                     will be disabled)
% Options.randdist  - If set to 1, randomly generated additive disturbance
%                       vector will be added to the state equation
% Options.lpsolver  - Solver for LPs (see help mpt_solveLP)
% Options.abs_tol   - absolute tolerance
% Options.openloop  - If 1, the open-loop solution will be computed
% Options.minnorm   - If closed-loop trajectory is computed, we stope evolution
%                       if norm of a state decreases below this value
% Options.verbose   - Level of verbosity
% Options.horizon   - If infinity, computes evolution of the state to origin.
%                       Set this value to a positive integer if you want
%                       to plot just the first Options.horizon points
% Options.newfigure - If set to 1, opens a new figure window
% Options.legend    - If set to 1, legend will be verboseed
% Options.showPn    - If set to 1 (default), plots the polyhedral partition over
%                     which the control law is defined.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% X, U       - matrices containing evolution of states and control inputs
% cost       - contains cost from the given initial state to the origin
% trajectory - indices specifying in which region of PA the given state lies
%
% see also MPT_COMPUTETRAJECTORY, MPT_PLOTTIMETRAJECTORY
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
    error('mpt_plotTrajectory: Wrong number of input arguments!');
end

if nargin<2,
    Options=[];
end

Options = mpt_defaultOptions(Options, ...
    'openloop', 0, ... % compute the open-loop solution? (values 1/0, where 0 represents the closed-loop solution)
    'maxCtr', 200, ... % maximum namer of iterations in the closed-loop control-law evaluation
    'verbose', mptOptions.verbose, ...
    'horizon', Inf, ... 
    'newfigure', mptOptions.newfigure, ...
    'randdist', 1, ... % if random disturbances should be added to the state vector
    'legend', 0, ...
    'x0', [], ...
    'showPn', 1 );

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl)
        error('This function supports only explicit controllers. Use mpt_plotTimeTrajectory() instead.');
    end
end

ctrlStruct = ctrl;
if isstruct(ctrlStruct)
    if ~mpt_isValidCS(ctrlStruct)
        error('mpt_plotTrajectory: First argument has to be a valid controller structure! See mpt_control for details.');
    end
elseif ~isa(ctrlStruct, 'mptctrl')
    error('MPT_PLOTTRAJECTORY: unknown type of first argument!');
end

    
Options.closedloop_color='k';
Options.openloop_color='r';

X=[]; U=[]; cost=0; trajectory=[];

sysStruct = ctrlStruct.sysStruct;
probStruct = ctrlStruct.probStruct;

opt.verbose=0;
if ~isfield(sysStruct,'verified'),
    sysStruct=mpt_verifySysStruct(sysStruct,opt);
end

if ~isfield(probStruct,'verified'),
    probStruct=mpt_verifyProbStruct(probStruct,opt);
end

if Options.horizon<=0,
    error('mpt_plotTrajectory: Options.horizon MUST be a positive integer!');
end

if ~isinf(Options.horizon),
    Options.minnorm=-1;
    Options.maxCtr=Options.horizon;
end

PA = ctrlStruct.Pn;
dimPA = dimension(PA(1));
if probStruct.tracking,
    if isfield(ctrlStruct.sysStruct, 'dims'),
        nu = ctrlStruct.sysStruct.dims.nu;
        nx = ctrlStruct.sysStruct.dims.nx;
        nyt = ctrlStruct.sysStruct.dims.ny;
    else
        if iscell(sysStruct.B),
            if(isfield(sysStruct,'Dy') & isfield(probStruct,'Qy'))
                nu = size(sysStruct.B{1},2);
                nx = size(sysStruct.A{1},2)-size(sysStruct.Dy{1},1);
            else 
                nu = size(sysStruct.B{1},2);
                if probStruct.tracking==1,
                    nx = (size(sysStruct.A{1},2)-nu)/2;
                else
                    nx = size(sysStruct.A{1},2)/2;
                end
            end
        else
            if(isfield(sysStruct,'Dy') & isfield(probStruct,'Qy'))
                nu = size(sysStruct.B,2);
                nx = size(sysStruct.A,2)-size(sysStruct.Dy,1);
            else
                nu = size(sysStruct.B,2);
                if probStruct.tracking==1,
                    nx = (size(sysStruct.A,2)-nu)/2;
                else
                    nx = size(sysStruct.A,2)/2;
                end
            end
        end
    end
    Options.xsection = nx+1:dimPA;
    Options.verbose = 0;
    dimPA = nx;
    Options.shade=1;
end
if isfield(ctrlStruct.sysStruct,'dumode')
    % if this flag is set, solution has been computed for extended state-space
    % to guarantee fullfilment of deltaU constraints in closed-loop
    if iscell(ctrlStruct.sysStruct.B),
        nu = size(ctrlStruct.sysStruct.B{1}, 2);
        nx = size(ctrlStruct.sysStruct.A{1}, 2) - nu;
    else
        nu = size(ctrlStruct.sysStruct.B, 2);
        nx = size(ctrlStruct.sysStruct.A, 2) - nu;
    end
    Options.xsection = nx+1:nx+nu;
    Options.verbose = 0;
    dimPA = nx;
    Options.shade=1;
end

if dimPA~=2,
    error('mpt_plotTrajectory: Sorry, only 2-dimensional trajectories can be plotted. Please use mpt_plotTimeTrajectory instead.');
end

if ~isempty(Options.x0),
    Options.x0 = Options.x0(:);
    if length(Options.x0)~=dimPA,
        error(['mpt_plotTrajectory: Options.x0 has to be a ' num2str(dimPA) 'x1 vector!']);
    end
end

if Options.showPn,
    h=mpt_plotPartition(ctrlStruct,Options);
end
grid on;
hold on
if probStruct.tracking,
    fprintf('\n');
    Options.reference = input('Enter reference state/output:   ');
    if(size(Options.reference,1)<size(Options.reference,2))
        Options.reference=Options.reference';
    end
    
    %     title('Click on the figure to specify the reference point (right-click to abort)');
    %     [x,y,button] = ginput(1);
    %     if button==3,        % if right-mouse-button pressed, abort
    %         if nargout==0,
    %             clear X U cost trajectory
    %         end
    %         return
    %     end
    %     Options.reference = [x;y];
    %    plot(x,y,'wo','LineWidth',3);
    
    if(nx>length(Options.reference))
        plot(Options.reference(1),0,'wo','LineWidth',3);
    else
        plot(Options.reference(1),Options.reference(2),'wo','LineWidth',3);
    end
end

if isempty(Options.x0)
    fprintf('\n');
    disp('Please click on the figure to specify the initial state (right-click to abort)');
    title('Please click on the figure to specify the initial state (right-click to abort)')
    hold on
end
while 1,
    if isempty(Options.x0)
        [x,y,button]   =   ginput(1);   %graphically enter one point
        if button==3,        % if right-mouse-button pressed, abort
            if nargout==0,
                clear X U cost trajectory
            end
            hold off
            return
        end
        x0 = [x;y];          % get mouse pointer coordinates
    else
        x0 = Options.x0;
    end
    handle_trajectory=[];
    maxCtr=100;
    index=1;

    if Options.openloop,
        legendText{index}='Open-Loop Trajectory';
        color=Options.openloop_color;
    else
        legendText{index}='Closed-Loop Trajectory';
        color=Options.closedloop_color;
    end
    index=index+1;

    % compute evolution of the system from the given initial state
    [X,U,Y,D,cost,trajectory,feasible]=mpt_computeTrajectory(ctrlStruct,x0,Options.horizon,Options);
    if feasible | size(X,1)>1,
        % if the problem is feasible
        if Options.verbose>=2,
            if Options.openloop,
                disp(['The trajectory cost of the openloop system is ' num2str(cost)]);
            else
                disp(['The trajectory cost of the closed-loop system  is ' num2str(cost)]);
            end
        end
        %plot results
        for i=1:(size(X,1)-1)
            h=gcf;
            handle_trajectory=plot(X(i:i+1,1),X(i:i+1,2),[color ':'],'LineWidth',2.5);  
            plot(X(i:i+1,1),X(i:i+1,2),[color '*'],'LineWidth',3);
        end   
        if ~Options.legend,
            title([legendText{1} ' for initial state [' num2str(x0(1)) ',' num2str(x0(2)) ']'],'FontSize',16);
        else
            title(['Trajectory for initial state [' num2str(x0(1)) ',' num2str(x0(2)) ']'],'FontSize',16);
        end
        if isfield(ctrlStruct.sysStruct,'StateName'),
            xlabel(ctrlStruct.sysStruct.StateName{1},'Fontsize',16); % LaTeX math symbols are directly supported!
            ylabel(ctrlStruct.sysStruct.StateName{2},'Fontsize',16);
        else
            xlabel('x_1','Fontsize',16); % LaTeX math symbols are directly supported!
            ylabel('x_2','Fontsize',16);
        end
        h=gcf;
        h1 = get(h,'CurrentAxes');
        set(h1,'Fontname','times');
        set(h1,'Fontsize',14);
        if Options.legend,
            legend(handle_trajectory,legendText);
        end
        if ~feasible,
            disp(['State [' num2str(x0(1)) ',' num2str(x0(2)) '] has no associated control law at step ' num2str(size(X,1)) ' !']);
        end
    else
        disp(['State [' num2str(x0(1)) ',' num2str(x0(2)) '] has no associated control law at step ' num2str(size(X,1)) ' !']);
    end
    %compute infinite horizon cost
    if ~isempty(Options.x0),
        break
    end
end 
if nargout==0,
    clear X U cost trajectory
end
