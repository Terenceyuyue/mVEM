function [islyapfun,XC,decayrate] = mpt_checkLyapFct(varargin)
%MPT_CHECKLYAPFCT Checks if a function is a Lyapunov function for a given ctrlStruct
%
% for PWQ Lyapunov functions:
%   islyapfun = mpt_checkLyapFct(ctrlStruct,LQ,LL,LC,Options)
%
% for PWA Lyapunov functions:
%   islyapfun = mpt_checkLyapFct(ctrlStruct,LL,LC,Options)
%
% for Polynomial Lyapunov functions:
%   islyapfun = mpt_checkLyapFct(ctrlStruct,V,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Grids the state-space into given number of points and for each state
% checks if decay rate of PWQ(PWA) Lyapunov function is negative.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl               - Explicit controller (MPTCTRL object)
% LQ,LL,LC           - Parameters of the PWQ(PWA) Lyapunov function
% Options.gridpoints - Number of grid points (if not provided, 30 is default)
% Options.N          - Number of simulation steps (5 by default)
% Options.verbose    - Level of verbosity
% Options.Pfinal     - polytope defining part of the state-space which should
%                      be considered for cost computation. (only reasonable
%                      if ctrlStruct.Pfinal is an empty polytope)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% islyapfun - 1 if the function satifies criterions put on Lyapunov-type
%             functions (positivity and negative decay rate); 0 otherwise
% XC        - matrix of data points (row-wise)
% decayrate - vector of decay rates corresponding to rows of XC
%
% see also MPT_GETPWQLYAPFCT. MPT_GETPWALYAPFCT
%

% Copyright is with the following author(s):
%
%(C) 2005 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%         grieder@control.ee.ethz.ch
%(C) 2004-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%              kvasnica@control.ee.ethz.ch

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

error(nargchk(2,6,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

Options = [];

if ~iscell(varargin{2}),
    error('mpt_checkLyapFct: Wrong parameters!');
end

ispoly= 0;
ispwa = 0;
ispwq = 0;
ctrl = varargin{1};

islyapfun = 1;

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl),
        error('This function supports only explicit controllers!');    
    end
    ctrlStruct = struct(ctrl);
else
    ctrlStruct = ctrl;
end
if nargin==2,
    % ctrlStruct, polyLyap
    ispoly = 1;
elseif nargin==3,
    if isstruct(varargin{3}),
        ispoly = 1;
        Options = varargin{3};
    elseif iscell(varargin{3}),
        % ctrlStruct, LL, LC
        ispwa = 1;
    else
        error('mpt_checkLyapFct: Wrong parameters!');
    end
elseif nargin==4,
    % ctrlStruct, LQ, LL, LC
    % or
    % ctrlStruct, LL, LC, Options
    if isstruct(varargin{4}),
        Options = varargin{4};
    elseif iscell(varargin{4}),
        ispwq = 1;
    else
        error('mpt_checkLyapFct: Wrong parameters!');
    end
elseif nargin==5,
    % ctrlStruct, LQ, LL, LC, Options
    if ~iscell(varargin{4}) | ~isstruct(varargin{5}),
        error('mpt_checkLyapFct: Wrong parameters!');
    end
    ispwq = 1;
    Options = varargin{5};
end

if(ispoly)
    polyLyap =varargin{2};
elseif ispwa,
    LL = varargin{2};
    LC = varargin{3};
elseif ispwq
    LQ = varargin{2};
    LL = varargin{3};
    LC = varargin{4};
end

if ~isfield(Options,'verbose')
    Options.verbose=mptOptions.verbose;
end
if ~isfield(Options,'N')
    Options.N = 5;
end
if ~isfield(Options,'gridpoints')
    Options.gridpoints = 20;
    if Options.verbose>0,
        fprintf('mpt_checkLyapFct: Number of grid points not given, assuming %d\n', Options.gridpoints);
    end

end
if ~isfield(Options, 'abs_tol')
    Options.abs_tol = mptOptions.abs_tol;
end

if ~mpt_isValidCS(ctrlStruct)
    error('mpt_checkLyapFct: First argument has to be a valid controller structure! See mpt_control for details.');
end

if isfield(Options,'Pfinal') %& ~isfulldim(ctrlStruct.Pfinal),
    if ~isa(Options.Pfinal,'polytope'),
        error('mpt_checkLyapFct: Options.Pfinal must be a polytope object!');
    end
    ctrlStruct.Pfinal = Options.Pfinal;
end

% first compute bounds on feasible state-space
if ~isfulldim(ctrlStruct.Pfinal),
    error('mpt_checkLyapFct: Please limit the state-space of interest by Options.Pfinal !');
end

if(ispoly)
    nx=dimension(ctrlStruct.Pfinal);
	%extract power of polynomial
	x = sdpvar(nx,1);
	found_deg=0;
	for deg=[1 2 4 6 8]
        xmon = monolist(x,deg);
        if(length(xmon)==length(polyLyap{1}))
            found_deg=1;
            break
        end
	end
	if(~found_deg)
        error('The Lyapunov function you wish to analyze is not in the correct format.')
	end
end

% get states which lie in the feasible set of a given controller
datapoints = mpt_feasibleStates(ctrl,Options.gridpoints,Options);


npoints = size(datapoints,1);
fprintf('Found %d datapoints in feasible set\n', npoints);
locOpt = Options;
locOpt.verbose = -1;
locOpt.fastbreak = 0;
locOpt.openloop = 0;

cost = zeros(npoints,1);

% if iterative solution was computed, regions are already order in such a way,
% that the first region found is in fact the one with least step distance to the
% origin, thus we can break in mpt_getInput as soon as at least one region is
% found (results in much faster run-time!
locOpt.fastbreak = (ctrlStruct.probStruct.subopt_lev>0);

if iscell(ctrlStruct.sysStruct.A),
    nx = size(ctrlStruct.sysStruct.A{1},2);
else
    nx = size(ctrlStruct.sysStruct.A,2);
end

decayrate = [];
XC = [];
drate_pos = 0;

locOpt = Options;
locOpt.verbose = 0;
isinOpt = Options;
isinOpt.fastbreak = 1;

allfeasible = 1;

for ii=1:npoints,
    if ii==1 | ii==npoints | mod(ii,round(npoints/4))==0,
        if Options.verbose > 0
            fprintf('%d/%d         \r',ii,npoints);
        end
    end

    x0 = datapoints(ii,:)';     % one state vector
    
    [X,U,Y,D,cost,traj,feas,dyns,reason,details] = mpt_computeTrajectory(ctrl, x0, Options.N, locOpt);

    if ~feas,
        allfeasible = 0;
        islyapfun = 0;
        disp('mpt_checkLyapFct: Closed-loop system is not invariant! Use ''mpt_invariantSet(ctrl)'' first.');
        return
    end
    if length(traj)<2 & feas==0,
        continue
    end
    if length(details{1})>1 & length(details{2})>1
        % initial state as well as the subsequent state lie on a boundary of two
        % regions, skip to next initial state
        continue
    end
    lyapval = [];
    region = [];
    for ix = 1:length(traj),
        inwhich = details{ix};
        lvmin = Inf;
        for iiw = 1:length(inwhich)
            region = inwhich(iiw);
            if ispoly,
                r_polyLyap = polyLyap{region};
                lv = lyapeval(X(ix,:)',1,{r_polyLyap},x,xmon);
            elseif ispwa,
                rLC = LC{region};
                rLL = LL{region};
                lv = lyapeval(X(ix,:)',1,{rLC},{rLL});
            elseif ispwq
                rLQ = LQ{region};
                rLC = LC{region};
                rLL = LL{region};
                lv = lyapeval(X(ix,:)',1,{rLC},{rLL},{rLQ});
            else
                error('mpt_checkLyapFct: Unknown type of lyapunov function');
            end
            if lv<lvmin
                lvmin = lv;
            end
        end
        lv = lvmin;
        ispositive = (lv >= Options.abs_tol*(X(ix,:)*X(ix,:)') - eps);
        if ~ispositive,
            fprintf('Function is not positive for state %s, aborting...\n',mat2str(X(ix,:), 5));
            islyapfun=0;
            XC = [];
            decayrate = [];
            return
        end
        lyapval = [lyapval; lv];
        if length(inwhich)>1,
            % do not compute decay rates further if the state is on a boundary
            % of two or more regions
            break
        end
    end
    drate = diff(lyapval);
    
    for ip = 1:length(drate),
        if drate(ip) > ( -Options.abs_tol * (X(ip,:)*X(ip,:)') ) + eps,
            fprintf('Decay rate=%f is not negative for state %s, aborting...\n',drate(ip),mat2str(X(ip,:), 5));
            XC = [];
            decayrate = drate(ip);
            islyapfun = 0;
            return
        end
    end

    if ~isempty(drate)
        decayrate = [decayrate; drate(1)];
        XC = [XC; x0'];
    end
end

if ~allfeasible
    fprintf('\nClosed-loop system is not invariant!!!\n');
end

if nargout==0,
    clear XC decayrate
end
return

%--------------------------------------------------------------
function val = lyapeval(x0,region,LC,LL,LQ),
% evaluates a PWA(PWQ) lyapunov function at point x0

if(~iscell(LL))
    x=LL;
    xmon=LQ;
    %evaluate Lyapunov function at state x1
    setsdpvar(x,x0);
    val=double(LC{region}*xmon);
elseif nargin==4,
    % PWA Lyapunov function
    val = LL{region}(:)'*x0 + LC{region};
elseif nargin==5,
    val = x0'*LQ{region}*x0 + LL{region}(:)'*x0 + LC{region};
end
