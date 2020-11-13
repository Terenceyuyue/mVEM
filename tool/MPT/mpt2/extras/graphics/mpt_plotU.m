function handle=mpt_plotU(ctrl,uind,Options)
%MPT_PLOTU For a given explicit controller, plots value of the control action
%
% handle=mpt_plotU(ctrl,Options)
% handle=mpt_plotU(ctrl,uind,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots value of the PWA control law with respect to a polyhedral partition
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl                    - Explicit controller (MPTCTRL object)
% uind                    - Which output to plot (for MIMO systems)
% Options.extreme_solver  - Which method to use for vertex enumeration (see help extreme)
% Options.lpsolver        - LP solver to be used
% Options.abs_tol         - absolute tolerance
% Options.newfigure       - If set to 1, opens a new figure window
% Options.verbose         - Level of verbosity
% Options.showPn          - If set to 1, plots on polyhedral sets Pn,
%                           (default: 0)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% handle  - handle of the plot
%
% see also MPT_PLOTPWA
%

% Copyright is with the following author(s):
%
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2004 Arne Linder, Faculty of Electrical, Information and Media Engineering, 
%          Wuppertal University, alinder@uni-wuppertal.de
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

error(nargchk(1,3,nargin));

if nargin<3,
    Options=[];
end
if nargin<=1,
    uind=[];
end
if nargin==2,
    if isstruct(uind),
        Options=uind;
        uind = [];
    elseif ~isa(uind,'double')
        error('Second argument must be either Options structure or index of control input!');
    end
end

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl),
        error('This function supports only explicit controllers!');
    end
    ctrlStruct = struct(ctrl);
else
    ctrlStruct = ctrl;
end

global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;  % absolute tolerance
end
if ~isfield(Options, 'axis')        % axis for plot
    Options.axis='auto';
end
if ~isfield(Options, 'lpsolver')        % axis for plot
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options, 'extreme_solver')        % axis for plot
    Options.extreme_solver=mptOptions.extreme_solver;
end
if ~isfield(Options, 'showPn')
    Options.showPn=0;
end
if ~isfield(Options,'newfigure')
    Options.newfigure=mptOptions.newfigure;
end
if ~isfield(Options,'verbose')
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'clip_minmax')
    Options.clip_minmax = [];
end


if ~mpt_isValidCS(ctrlStruct)
    error('mpt_plotU: First argument has to be a valid controller structure! See mpt_control for details.');
end

index=0;

handle = [];

if iscell(ctrlStruct.sysStruct.B),
    nu = size(ctrlStruct.sysStruct.B{1},2);
else
    nu = size(ctrlStruct.sysStruct.B,2);
end

if iscell(ctrlStruct.sysStruct.A),
    nx = size(ctrlStruct.sysStruct.A{1},2);
else
    nx = size(ctrlStruct.sysStruct.A,2);
end

PA = ctrlStruct.Pn;
Fi = ctrlStruct.Fi;
Gi = ctrlStruct.Gi;

if isempty(uind)
    % plot all inputs in separe sub-plots
    urange = 1:nu;
else
    if ~ismember(uind,1:nu)
        error('Index of control move outside of range!');
    end
    urange = uind;
end

maxlen=length(PA);
minu=Inf;
locOpt = Options;
locOpt.openloop = 0;

disp('Plotting value of control moves...');

if ctrlStruct.overlaps & Options.verbose>0,
    fprintf('\n');
    disp('-----------------------------------------------------------------------------');
    disp('mpt_plotU warning: regions do overlap, control action is not uniquely defined!');
    disp('Please use mpt_removeOverlaps first to get proper result.');
    disp('-----------------------------------------------------------------------------');    
    fprintf('\n');
end

if ~ctrlStruct.overlaps,
    % if regions do not overlap, we can break in isinside as soon as the first region has been idintified
    Options.fastbreak=1;
else
    Options.fastbreak=0;
end

Options.fastbreak=1;


if nx==1
    for num_u = urange,
        % handle range of inputs
        if isempty(uind) & nu > 1,
            % open new subplot only if we have more than one input
            subplot(nu, 1, num_u);
        end
        Fi_u = cell(1, maxlen);
        Gi_u = cell(1, maxlen);
        for ir = 1:maxlen,
            % extract feedback law associated to a given input
            Fi_u{ir} = Fi{ir}(num_u:num_u);
            Gi_u{ir} = Gi{ir}(num_u:num_u);
        end
        handle = mpt_plotPWA(PA,Fi_u,Gi_u,Options);
        
        title(sprintf('Value of the control action U^{*}_%d over %d regions', ...
            num_u,length(PA)),'FontSize',14);
        xlabel('x_1','Fontsize',14); % LaTeX math symbols are directly supported!
        ylabel(sprintf('U^{*}_%d(x)', num_u),'Fontsize',14);
        grid on; 
        h=gcf;
        h1 = get(h,'CurrentAxes');
        set(h1,'Fontname','times');
        set(h1,'Fontsize',14);
    end
    
else
    for num_u=urange            % Addition by Arne Linder for dealing with models with multiple inputs
        if isempty(uind) & nu > 1,
            % open new subplot only if we have more than one input            
            subplot(nu,1,num_u);
        end
        
        for ii=1:maxlen
            P=PA(ii);
            nb=nconstr(P);       % number of constraints
            dimP=dimension(P);   % dimension
            if dimP~=2,
                close
                error('mpt_plotU: Only 2D partitions can be plotted by this function!');
            end
            [xc,rc]=chebyball(P);  % get chebyshev's radius
            if rc<=Options.abs_tol
                disp('mpt_plotU: Empty polytope detected!');
                continue
            end
            axis(Options.axis);
            
            [V,R,PA(ii)]=extreme(P,Options);    % compute extreme points and rays of polytope P
            if size(R,1)>0
                error('mpt_plotU: Polytope is unbounded!'); % existence of rays means polytope is unbounded
            end
            
            % sort vertices in a cyclic way;
            x1=V(:,1);
            x2=V(:,2);
            
            ang=angle([(x1-xc(1))+(x2-xc(2))*sqrt(-1)]);
            [val,ind]=sort(ang);
            x1=x1(ind);
            x2=x2(ind);
            x3=[];
            for jj=1:length(x1),
                x = [x1(jj); x2(jj)];
                [feasible,inwhich] = isinside(ctrlStruct.Pn,x,Options);
                inwhich=ii;
                %[optU,feasible] = mpt_getInput(ctrlStruct,x,locOpt);
                if ~feasible,
                    error('mpt_plotU: problem detected! Increase value of abs_tol in the file mpt_init!');
                end
                optU = ctrlStruct.Fi{inwhich}(num_u,:)*x + ctrlStruct.Gi{inwhich}(num_u,:); % Modified by Arne Linder
                if ~isempty(Options.clip_minmax)
                    optU = max(min(optU, Options.clip_minmax(:, 2)), ...
                        Options.clip_minmax(:, 1));
                end
                x3 = [x3; optU]; % the third dimension will be value of the control action, i.e. u=Fx+G
            end
            % x3=[x1 x2]*Fi{ii}(1,:)'+Gi{ii}(1,:);   
            if min(x3)<minu,
                minu=min(x3);
            end
            
            h=patch(x1,x2,x3,x3);
            if isfield(Options, 'LineWidth'),
                set(h, 'LineWidth', Options.LineWidth);
            end
            handle=[handle;h];
        end
        if Options.showPn,
            % plot the full polyhedral partition PA below
            for ii=1:maxlen,
                P=PA(ii);
                nb=nconstr(P);
                dimP=dimension(P);
                [xc,rc]=chebyball(P);
                if rc<=Options.abs_tol
                    disp('mpt_plotU: Empty polytope detected!');
                    continue
                end
                axis(Options.axis);
                
                [V,R,PA(ii)]=extreme(P,Options);
                if size(R,1)>0
                    error('mpt_plotU: Polytope is unbounded!');
                end
                
                % sort vertices in a cyclic way;
                x1=V(:,1);
                x2=V(:,2);
                
                ang=angle([(x1-xc(1))+(x2-xc(2))*sqrt(-1)]);
                [val,ind]=sort(ang);
                x1=x1(ind);
                x2=x2(ind);
                x3=[];
                for jj=1:length(x1),
                    x = [x1(jj); x2(jj)];
                    [optU,feasible] = mpt_getInput(ctrlStruct,x,locOpt);
                    if ~feasible,
                        error('mpt_plotU: problem detected! Increase value of abs_tol in mpt_init!');
                    end
                    if ~isempty(Options.clip_minmax)
                        optU = max(min(optU, Options.clip_minmax(:, 2)), ...
                            Options.clip_minmax(:, 1));
                    end
                    x3 = [x3; optU]; % the third dimension will be value of the control action, i.e. u=Fx+G
                end
                %x3=[x1 x2]*Fi{ii}(1,:)'+Gi{ii}(1,:);  % third dimension
                %will be value of the control action, i.e. u=Fx+G
                if x3<minu,
                    minu=x3;
                end
                h=patch(x1,x2,minu*ones(size(x1)),'b');
                handle=[handle;h];
            end    
        end
        view(3);
        title(sprintf('Value of the control action U_%d over %d regions',num_u,length(PA)),'FontSize',14); % Modified by Arne Linder
        xlabel('x_1','Fontsize',14); % LaTeX math symbols are directly supported!
        ylabel('x_2','Fontsize',14);
        zlabel(sprintf('U^{*}_%d(x)',num_u),'Fontsize',14);
        grid;
        h=gcf;
        h1 = get(h,'CurrentAxes');
        set(h1,'Fontname','times');
        set(h1,'Fontsize',14);
    end
end

% If no outputs is asked, clear the variable.
if nargout == 0;
    clear('handle');
end
axis tight
