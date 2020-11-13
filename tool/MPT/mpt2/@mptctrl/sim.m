function [X,U,Y,cost,feasible]=sim(ctrl, arg2, arg3, arg4, arg5)
%SIM Simulates a given controller
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Once initial conditions x0 are known, it is possible to simulate a given
% controller in closed-loop by running
%
%       sim(ctrl, x0)
%
% Here the simulation will start from a given state "x0" and will terminate once
% a given control objective is reached (usually when the origin is reached).
%
% Additionally, you can specify for how many steps you want to perform the
% simulation:
%
%       sim(ctrl, x0, N)
%
% Additional options can be provided as well:
%
%       sim(ctrl, x0, N, Options)
%
% Hint: set N=[] in the above command if you want to specify Options but do not
% want to specify number of simulation steps.
%
% You can also specify your own dynamical system to use for simulations. In such
% case control actions obtained by a given controller can be applied to a
% different system than that which was used for computing the controller:
%
%       sim(ctrl, system, x0, N, Options)
%
% Note that the "N" and "Options" arguments are optional. You can specify your
% own dynamics in two ways: 
% 1. by setting the "system" parameter to a system structure:
%
%       sim(ctrl, sysStruct, x0, N, Options)
%
% 2. by setting the "system" parameter to a handle of a function which will
%    provide updates of system states in a discrete-time fashion:
%        
%       sim(ctrl, @sim_function, x0, N, Options)
%
%    Take a look at 'help di_sim_fun' on how to write simulation functions
%    compatible with this function.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl    - An MPTCTRL object
% system  - System model to use for simulations. It can either be a sysStruct
%           structure or a handle of a function which takes "x_k" and "u_k" as
%           input arguments and produces "x_k+1" (state update) and "y_k"
%           (system outputs) vectors.
% x0      - Initial state
% N       - Number of simulation steps. If not specified, or defined as an empty
%           matrix ([]) or as Inf, simulation is performed until a given
%           regulation objective is reached (usually when the origin is reached)
% Options - Additional options. See 'help mpt_computeTrajectory' and 
%           'help mpt_getInput' for list of supported options.
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% X,U,Y    - Matrices which contain evolution of system states, inputs and
%            outputs during the simulation
% cost     - Cost associated to the simulated trajectory
% feasible - A true/false flag indicating whether the simulation was feasible
%
% see also MPT_COMPUTETRAJECTORY, MPT_GETINPUT, MPTCTRL/SIMPLOT
%

% Copyright is with the following author(s):
%
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

error(nargchk(2,5,nargin));

if ~isa(ctrl, 'mptctrl')
  error('SIM: First argument must be an MPTCTRL object.');
end

% sys_type = 0    - use sysStruct stored in ctrl
% sys_type = 1    - sysStruct provided
% sys_type = 2    - function handle to simulation function provided
sys_type = 0; 

switch nargin,
    case 2,
        % sim(ctrl, x0)
        
        if ~isa(arg2, 'double'),
            error('SIM: second input must be vector of initial conditions.');
        end
        x0 = arg2;
        
    case 3,
        % sim(ctrl, x0, N)
        % sim(ctrl, x0, Options)
        % sim(ctrl, system, x0)        
        
        % detect type of second input. it can either be the "x0" vector or the
        % "system" pointed (either sysStruct structure or a function handle)
        if isa(arg2, 'double'),
            sys_type = 0;
            x0 = arg2;
        elseif isa(arg2, 'struct'),
            sys_type = 1;
            sysStruct = arg2;
        elseif isa(arg2, 'function_handle'),
            sys_type = 2;
            sysHandle = arg2;
        else
            error('SIM: wrong type of second input argument.');
        end
        
        % if x0 was provided as 2nd input, 3rd input must be either "N" or
        % "Options":
        if sys_type==0,
            % corresponds to either one of these two:
            % sim(ctrl, x0, N)
            % sim(ctrl, x0, Options)

            if isa(arg3, 'double'),
                % 3d input is the horizon
                N = arg3;
            elseif isa(arg3, 'struct'),
                % 3d input is the Options structure
                Options = arg3;
            else
                error('SIM: wrong type of third input argument.');
            end
            
        else
            % corresponds to "sim(ctrl, system, x0)"
            % 3rd input must be x0
            if isa(arg3, 'double'),
                x0 = arg3;
            else
                error('SIM: wrong type of third input argument.');
            end
        end
        
    case 4,
        % sim(ctrl, x0, N, Options)
        % sim(ctrl, system, x0, N)
        
        % detect type of second input. it can either be the "x0" vector or the
        % "system" pointed (either sysStruct structure or a function handle)
        if isa(arg2, 'double'),
            sys_type = 0;
            x0 = arg2;
        elseif isa(arg2, 'struct'),
            sys_type = 1;
            sysStruct = arg2;
        elseif isa(arg2, 'function_handle'),
            sys_type = 2;
            sysHandle = arg2;
        else
            error('SIM: wrong type of second input argument.');
        end
        
        if sys_type==0,
            % sim(ctrl, x0, N, Options)
            if ~isa(arg3, 'double'),
                error('SIM: wrong type of third input argument.');
            end
            N = arg3;
            if any(size(N)>1),
                error('SIM: third input must be a scalar.');
            end
            if ~isempty(N),
                if N<1 | round(N)~=N,
                    error('SIM: number of simulation steps must be a positive integer.');
                end
            end
            
            if isempty(arg4) | isstruct(arg4),
                Options = arg4;
            else
                error('SIM: fourth input must be an Options structure.');
            end
            
        else
            % sim(ctrl, system, x0, N)
            % sim(ctrl, sysStruct, x0, Options)
            if isa(arg3, 'double')
                x0 = arg3;
            else
                error('SIM: third input must be vector of initial conditions.')
            end
            
            if isa(arg4, 'double'),
                N = arg4;
            elseif isa(arg4, 'struct'),
                N = [];
                Options = arg4;
            else
                error('SIM: wrong type of fourth input argument.');
            end
            
            if any(size(N)>1),
                error('SIM: fourth input must be a scalar.');
            end
            if ~isempty(N),
                if N<1 | round(N)~=N,
                    error('SIM: number of simulation steps must be a positive integer.');
                end
            end
        end
        
    case 5,
        % sim(ctrl, system, x0, N, Options)
        
        % 2nd input must be the "system" argument
        if isa(arg2, 'struct'),
            sys_type = 1;
            sysStruct = arg2;
        elseif isa(arg2, 'function_handle'),
            sys_type = 2;
            sysHandle = arg2;
        else
            error('SIM: wrong type of second input argument.');
        end
        
        % 3d input must be x0
        if isa(arg3, 'double'),
            x0 = arg3;
        else
            error('SIM: wrong type of third input argument.');
        end
        
        % 4th input must be N
        if ~isa(arg4, 'double'),
            error('SIM: wrong type of third input argument.');
        end
        N = arg4;
        if any(size(N)>1),
            error('SIM: fourth input must be a scalar.');
        end
        if ~isempty(N),
            if N<1 | round(N)~=N,
                error('SIM: number of simulation steps must be a positive integer.');
            end
        end
        
        % 5th input must be Options
        if isempty(arg5) | isstruct(arg5),
            Options = arg5;
        else
            error('SIM: fifth input must be an Options structure.');
        end
    
    otherwise,
        error('SIM: wrong number of input arguments.');

end

if exist('x0', 'var') ~= 1,
    error('SIM: initial state x0 must be provided.');
end
if exist('N', 'var') ~= 1,
    % let mpt_computeTrajectory abort the simulation when a given control
    % objective is reached
    N = [];
end
if exist('Options', 'var') ~= 1,
    Options = [];
end
if sys_type==1,
    % use auxiliary sysStruct
    Options.sysStruct = sysStruct;
    
elseif sys_type==2,
    % use auxiliary simulation function
    Options.sysHandle = sysHandle;
end

[X,U,Y,D,cost,trajectory,feasible] = sub_computeTrajectory(ctrl, x0, N, Options);
