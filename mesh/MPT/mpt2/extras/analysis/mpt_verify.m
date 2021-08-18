function [flag, N, Vreach] = mpt_verify(object, X0, Xf, N, U0, Options)
%MPT_VERIFY Verifies whether states enter a given set in a given number of steps
%
% [flag, N, Vreach] = mpt_verify(sysStruct, X0, Xf, N, U0)
% [flag, N, Vreach] = mpt_verify(controller, X0, Xf, N)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Checks if states of a dynamical system subject to:
%  1. inputs driven by an explicit control law
%  2. inputs which belong to a set of admissible inputs U0
% enter a given set Xf, assuming x0 \in X0 in N steps.
%
% USAGE:
%   flag = mpt_verify(sysStruct, X0, Xf, N)
%   flag = mpt_verify(sysStruct, X0, Xf, N, U0)
%   flag = mpt_verify(controller, X0, Xf, N)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% object       - either a sysStruct structure, or an explicit controller
% X0           - set of initial states (a polytope or a polyarray)
% Xf           - target set (a polytope or a polyarray)
% N            - number of steps over which to compute reachable sets
% U0           - set of admissible inputs as a polytope object. if omitted,
%                inputs will be constrained by sysStruct.umin/umax
% Options.usereachsets 
%              - if set to 0 (default is 1), verification question for systems
%                described by sysStruct structures will be performed by solving
%                an LP/MILP problem. Set this option to true if you want to base
%                the verification on reachable sets. All system constraints
%                as defined in sysStruct will be enforced.
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% flag        - 1 if Xf is reachable from X0, 0 otherwise
% N           - number of steps in which Xf is reachable from X0 ([] if Xf is
%               not reachable)
% Vreach      - V-representation of reachable sets (default), or the feasible
%               optimizer (if Options.usereachsets=0)
%
% see also MPT_REACHSETS, MPT_REACHXU
%

% Copyright is with the following author(s):
%
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

global mptOptions
if ~isstruct(mptOptions),
    mpt_error;
end

error(nargchk(4,6,nargin));

if isa(object, 'mptctrl'),
    if ~isexplicit(object),
        error('This function supports only explicit controllers!');
    end
end

if ~isa(X0, 'polytope'),
    error('Set of initial states must be a polytope!');
end
if ~isa(N, 'double'),
    error('Number of steps must be a positive integer!');
end
if N<1 | any(size(N)~=1),
    error('Number of steps must be a positive integer!');
end
if ~isa(Xf, 'polytope'),
    error('Set of final states must be a polytope!');
end
if ~isfulldim(X0),
    error('Set of initial states must be a fully dimensional polytope!');
end
if ~isfulldim(Xf),
    error('Set of final states must be a fully dimensional polytope!');
end
if nargin<=4,
    Options = [];
    U0 = [];
elseif nargin == 5,
    if isa(U0, 'struct'),
        Options = U0;
        U0 = [];
    else
        Options = [];
    end
end

% set default options
Options = mpt_defaultOptions(Options, ...
    'usereachsets', 1, ...
    'verbose', mptOptions.verbose );

% check if U0 is a fully dimensional polytope
if isa(U0, 'polytope'),
    if ~isfulldim(U0),
        error('Set of admissible inputs must be a fully dimensional polytope!');
    end
end

% run the verification
if mpt_issysstruct(object) & Options.usereachsets==0,
    % verification based on a feasibility LP/MILP
    if Options.verbose > 0,
        fprintf('Performing verification based on a feasibility LP/MILP...\n');
    end
    [flag, Vreach] = sub_feasibility_verification(object, X0, Xf, U0, N);
    
    if flag==0,
        % no feasible solution exists
        N = [];
    end
    
elseif ~mpt_issysstruct(object) & Options.usereachsets==0,
    error('Cannot use Options.usereachsets when input is a controller.');
    
else
    % otherwise use verification based on reachable sets
    if Options.verbose > 0,
        fprintf('Performing verification based on reachable sets...\n');
    end
    
    % tell mpt_reachSets to stop the evolution once the target is reached
    Options.Xf = Xf;
    
    % compute reachable sets
    if ~isa(U0, 'polytope'),
        [res, Vreach] = mpt_reachSets(object, X0, N, Options);
        
    else
        [res, Vreach] = mpt_reachSets(object, X0, U0, N, Options);
        
    end
    
    % reachable sets exist?
    if isa(res, 'double'),
        N = res;
        flag = 1;
    else
        N = [];
        flag = 0;
    end
    
end


%------------------------------------------------------------------------
function [feasible, optimizer] = sub_feasibility_verification(sysStruct, ...
    X0, Xf, U0, N);
% performs the verification by solving a feasibility problem if the input
% is a system structure

global mptOptions

% we need dimensions to set a proper problem structure
[nx, nu, ny, ndyn, nbool] = mpt_sysStructInfo(sysStruct);

% no target set now, will add it later
probStruct.Tconstraint = 0;

% weights can be arbitrary, we will solve feasibility problem
probStruct.Q = eye(nx);
probStruct.R = eye(nu);

% cost can be arbitrary, we won't use it
probStruct.norm = 1;

% prediction horizon
probStruct.N = N;
probStruct.subopt_lev = 0;

% construct constraints
[Con, Obj, Vars] = mpt_ownmpc(sysStruct, probStruct);

% add the constraint x(0) \in X0
Con = Con + set(ismember(Vars.x{1}, X0));

% add target set constraint x(N) \in Xf
% notice that YALMIP will automatically handle cases when Xf is a
% polytope array
Con = Con + set(ismember(Vars.x{end}, Xf));

% add polytopic constraint on all inputs if asked for
if isa(U0, 'polytope'),
    for i = 1:length(Vars.u),
        Con = Con + set(ismember(Vars.u{i}, U0));
    end
end

% finally solve the feasibility problem
Obj = [];         % no objective = just feasibility
yalmipOptions = mptOptions.sdpsettings;
sol = solvesdp(Con, Obj, yalmipOptions);

% verification successfull?
feasible = (sol.problem == 0);
if feasible,
    % yep, feasible optimizer exists, retrieve it
    optimizer = double(cat(1, Vars.u{:}));

else
    % nope, no feasible sequence of manipulated variables exists
    optimizer = [];

end
