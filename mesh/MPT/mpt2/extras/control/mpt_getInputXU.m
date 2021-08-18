function [U,feasible,region,XU_region]=mpt_getInputXU(ctrlStruct,x0,Options)
%MPT_GETINPUTXU For a given state, extracts the an output from an XUset
%
% [U,feasible,Pn_region,XU_region]=mpt_getInputXU(ctrl,x0,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% For the given state x0, this function extracts an output based on different 
% strategies from an XU-stability set given by means of the controller structure
% ctrlStruct. It is assumed that the controller partition is non-overlapping in 
% the X space
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl                   - MPT controller
% ctrl.details.XU        - provided structure of the XUset
% x0                     - initial state
% Options                - optional arguments
%   .strategy            - (global) strategy ontop of the local strategy to pick
%                          the controller from the XUset
%                          strategy = 'avg'     - average value of u (default)
%                                     'max'     - maximize u
%                                     'min'     - minimize u
%                                     'mean'    - mean value of u
%                                     'feas'    - feasible u
%   .local_strategy      - local strategy to pick the controller from the XUset
%                          (governs how to solve an LP)
%                          local_strategy = 'min'     - minimize u
%                                           'max'     - maximize u
%                                           'dual'    - minimize & maximize U
%                                           'feas'    - feasible u
%                                           'rand'    - find random U
%                           default: is chosen according to the global strategy
%   .abs_tol             - absolute tolerance
%   .verbose             - Level of verbosity
%
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% U         - control input computed based on different strategies
% feasible  - returns 1 if  the there is at least one control law associated to
%             a given state x0, 0 otherwise
% Pn_region - index of the first region of ctrl.Pn which contains the optimal 
%             control input associated to the given state x0
% XU_region - indices of ctrl.details.XU which contain the given state x0
%
% see also MPT_COMPUTETRAJECTORY, MPT_PLOTTIMETRAJECTORY, MPT_XUSETPWALYAP
%

% Copyright is with the following author(s):
%
% (c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (c) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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

error(nargchk(2,3,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3,
    Options=[];
end

if ~isfield(Options, 'lpsolver'),
    lpsolver = mptOptions.lpsolver;
else
    lpsolver = Options.lpsolver;
end
if ~isfield(Options, 'strategy'),
    % strategy to choose input:
    %   'avg'     - average value of u
    %   'max'     - maximize u
    %   'min'     - minimize u
    %   'mean'    - mean value of u
    %   'feas'    - feasible u
    strategy = 'avg';
else
    strategy = Options.strategy;
end
if ~isfield(Options, 'local_strategy'),
    % governs how to solve an LP:
    %   'min'   - minimize U
    %   'max'   - maximize U
    %   'dual'  - minimize & maximize U
    %   'rand'  - find random U
    %   'feas'  - just solve feasibility problem
    
    if strcmp(strategy,'max')
        local_strategy = 'max';
    elseif strcmp(strategy,'min')
        local_strategy = 'min';
    elseif strcmp(strategy,'feas')
        local_strategy = 'feas';
    else
        local_strategy = 'dual';
    end
else
    local_strategy = Options.local_strategy;
end

if isfield(ctrlStruct.details,'XU')
    if ~isfield(ctrlStruct.details.XU, 'XUset') | ~isfield(ctrlStruct.details.XU, 'Idx_orig'),
        error('MPT_GETINPUTXU: provide XUset and/or Idx_orig in ctrl.details.XU')
    end
else
    error('MPT_GETINPUTXU: provide XUset and/or Idx_orig in ctrl.details.XU')
end
XUset = ctrlStruct.details.XU.XUset;
XU_Idx_orig = ctrlStruct.details.XU.Idx_orig;

if ~isa(XUset, 'polytope'),
    error('MPT_GETINPUTXU: first input must be a polytope object.');
end
if ~isa(x0, 'double')
    error('MPT_GETINPUTXU: second argument must be vector!');
end
if nargin==4 & ~isstruct(Options)
    error('MPT_GETINPUTXU: third argument must be an Options structure!');
end

if ~isfulldim(XUset),
    error('MPT_GETINPUTXU: first and second arguments must be fully dimensional polytopes.');
end

if length(XUset) ~= length(XU_Idx_orig),
    error('MPT_GETINPUTXU: first and second arguments must have same length.');
end

x0 = x0(:);
if strcmpi(strategy, 'feas'),
    feas_strategy = 1;
    if strcmpi(local_strategy, 'dual'),
        % do not solve 2 LP's if we just want a feasible U
        local_strategy = 'feas';
    end
else
    feas_strategy = 0;
end

U_all = [];
region = 0;
inwhich = [];
XU_region = 0;

% first find which regions contain given state
[feasible, inwhich] = isinside(ctrlStruct.Pn, x0);
isin = 0;
if feasible
    Pn_region = inwhich(1);
    inwhich = find(XU_Idx_orig(:)'==Pn_region);
    if isempty(inwhich)
        isin=0;
    else
        isin=1;
   end
end

if ~isin,
    U = [];
    feasible = 0;
    return
end


[nx, nu, ny, ndyn, nbool] = mpt_sysStructInfo(ctrlStruct.sysStruct);

% the state belongs to at least one region

XU_region = [];
for ir = 1:length(inwhich),
    XUregion = XUset(inwhich(ir));
    [Hxu, Kxu] = double(XUregion);
    
    if strcmpi(local_strategy, 'dual'),
        % minimize U
        [u_min, f_min] = sub_getU(x0, nu, nx, Hxu, Kxu, 'min', lpsolver);
        [u_max, f_max] = sub_getU(x0, nu, nx, Hxu, Kxu, 'max', lpsolver);
        f_r = f_min & f_max;
        u_r = [u_min u_max];
    else        
        % solve an LP to get feasible U
        [u_r, f_r] = sub_getU(x0, nu, nx, Hxu, Kxu, local_strategy, lpsolver);
        u_r = u_r(:);
        f_max = 0; f_min = 0;
    end
    
    if ~isempty(u_r)
        XU_region = [XU_region inwhich(ir)];
    end
    
    if f_r == 0,
        % Region XUregion does not contain state x0, skip it...
        continue
    end
    
    U_all = [U_all u_r];
    
    if feas_strategy,
        % we do not need to traverse through remaining regions in this
        % strategy is selected
        U = U_all(:);
        feasible = 1;
        region = Pn_region;
        return
    end
end

% exit if no feasible control law was found
if isempty(U_all),
    U = [];
    feasible = 0;
    return
end

U_all = U_all';
% from U_all take the minimum / maximum according to given strategy

% Options.strategy - strategy to choose input:
%   'max'     - maximize u
%   'min'     - minimize u
%   'mean'    - mean value of u
%   'average' - average value of u
%   'feas'    - feasible u

switch strategy
    case 'min',
        U = min(U_all);
        U = U(:);
        
    case 'max',
        U = max(U_all);
        U = U(:);
        
    case 'avg',
        U = sum(U_all) / size(U_all, 1);
        U = U(:);
        
    case 'mean',
        U = mean(U_all);
        U = U(:);
        
    case 'rand',
        Umax = max(U_all);
        Umin = min(U_all);
        U = Umin(:) + rand*(Umax(:) - Umin(:));

end

% check if the U(x0) value is really feasible
isin = isinside(XUset(inwhich), [x0; U]);
if ~isin,
    % U(x0) is not feasible, take the minimum
    fprintf('MPT_GETINPUTXU: "%s" value is not feasible, taking minimum value...\n', strategy);
    U = min(U_all);
    U = U(:);
end
        

feasible = 1;
region = Pn_region;


%----------------------------------------------------------------
function [u, feasible] = sub_getU(x0, nu, nx, Hxu, Kxu, strategy, lpsolver)

switch strategy,
    case 'min'
        % minimize u
        f = ones(1, nu);
    case 'max'
        % maximize u
        f = -ones(1, nu);
    case 'feas',
        % just find any feasible u
        f = zeros(1, nu);
    case 'rand',
        % find random u
        f = randn(1, nu);
end

% make a slice of XU polytope to x0, so we are left with a polytope in U space
Hn = Hxu(:, nx+1:end);
Kn = Kxu - Hxu(:, 1:nx) * x0;

% use "rescue" mode, since the LP should always be feasible
[xopt,fval,lambda,exitflag,how]=mpt_solveLPs(f,Hn,Kn,[],[],[],lpsolver);
if strcmpi(how, 'ok'),
    u = xopt(:);
    feasible = 1;
else
    u = [];
    feasible = 0;
end
    