function [X,U,Y,cost,feasible]=simplot(ctrl, arg2, arg3, arg4, arg5)
%SIMPLOT Plots a simulated trajectory
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Plots trajectory of a given closed-loop system starting from a given initial
% state x0. If the state is not provided and the controller is an explicit
% control law in exactly 2D, it is possible to select the initial state by
% mouse, i.e.:
%
%    simplot(ctrl)   or   simplot(ctrl, Options)
%
% If the initial state is provided, evolution of system states, inputs and
% outputs is plotted versus the time axis:
%
%   simplot(ctrl, x0)
%
% By default the simulation is stopped when a given control objective is reached
% (usually when the state reaches the origin). It is, however, possible to tell
% the function to only compute the evolution for a given number of steps:
%
%   simplot(ctrl, x0, N)
%
% Additional options can also be specified:
%
%   simplot(ctrl, x0, Options)   or   simplot(ctrl, x0, N, Options)
%
% Hint: set N=[] in the above command if you want to specify Options but do not
% want to specify number of simulation steps.
%
% You can also specify your own dynamical system to use for simulations. In such
% case control actions obtained by a given controller can be applied to a
% different system than that which was used for computing the controller:
%
%       simplot(ctrl, system, x0, N, Options)
%
% Note that the "N" and "Options" arguments are optional. You can specify your
% own dynamics in two ways: 
% 1. by setting the "system" parameter to a system structure:
%
%       simplot(ctrl, sysStruct, x0, N, Options)
%
% 2. by setting the "system" parameter to a handle of a function which will
%    provide updates of system states in a discrete-time fashion:
%        
%       simplot(ctrl, @sim_function, x0, N, Options)
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
% see also MPT_PLOTTRAJECTORY, MPT_PLOTTIMETRAJECTORY, MPTCTRL/SIM
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


% possible calls:
%
% simplot(ctrl)
% simplot(ctrl, Options)
% simplot(ctrl, x0)
% simplot(ctrl, x0, Options)
% simplot(ctrl, x0, N)
% simplot(ctrl, system, x0)
% simplot(ctrl, x0, N, Options)
% simplot(ctrl, system, x0, N)
% simplot(ctrl, system, x0, Options)
% simplot(ctrl, system, x0, N, Options)

error(nargchk(1,5,nargin));

if ~isa(ctrl, 'mptctrl')
  error('SIM: First argument must be an MPTCTRL object.');
end

switch nargin,
    case 1,
        % sim(ctrl)

        [X, U, cost] = mpt_plotTrajectory(ctrl);
        Y = []; feasible = [];
        if nargout < 1,
            clear X U Y cost feasible
        end
        return
        
    case 2,
        % simplot(ctrl, Options)
        % simplot(ctrl, x0)

        if isa(arg2, 'double'),
            % simplot(ctrl, x0)
            [X, U, Y, cost, feasible] = mpt_plotTimeTrajectory(ctrl, arg2);
        
        elseif isstruct(arg2),
            % simplot(ctrl, Options)
            [X, U, cost] = mpt_plotTrajectory(ctrl, arg2);
            Y = []; feasible = [];
            if nargout < 1,
                clear X U Y cost feasible
            end
            return
            
        else
            error('SIMPLOT: initial state must be given.');
        end
            
    case 3,
        % simplot(ctrl, x0, Options)
        % simplot(ctrl, x0, N)
        % simplot(ctrl, system, x0)
        
        if isa(arg2, 'double') & isa(arg3, 'double'),
            % simplot(ctrl, x0, N)
            if any(size(arg3)>1) | round(arg3)~=arg3,
                error('SIMPLOT: number of simulation steps must be a positive integer.');
            end
            [X, U, Y, cost, feasible] = mpt_plotTimeTrajectory(ctrl, arg2, arg3);
            
        elseif isa(arg2, 'double') & isa(arg3, 'struct'),
            % simplot(ctrl, x0, Options)
            [X, U, Y, cost, feasible] = mpt_plotTimeTrajectory(ctrl, arg2, [], arg3);
            
        elseif isa(arg2, 'function_handle') | isa(arg2, 'struct'),
            % simplot(ctrl, system, x0)
            if ~isa(arg3, 'double'),
                error('SIMPLOT: initial state must be given.');
            end
            if isa(arg2, 'function_handle'),
                Options.sysHandle = arg2;
            else
                Options.sysStruct = arg2;
            end
            [X, U, Y, cost, feasible] = mpt_plotTimeTrajectory(ctrl, arg3, [], Options);
            
        else
            error('SIMPLOT: wrong type of input arguments.');
        end

    case 4,
        % simplot(ctrl, x0, N, Options)
        % simplot(ctrl, system, x0, N)
        % simplot(ctrl, system, x0, Options)
        
        if isa(arg2, 'double') & isa(arg3, 'double') & isa(arg4, 'struct'),
            % simplot(ctrl, x0, N, Options)
            [X, U, Y, cost, feasible] = mpt_plotTimeTrajectory(ctrl, arg2, arg3, arg4);
            
        elseif isa(arg2, 'function_handle') | isa(arg2, 'struct'),
            % simplot(ctrl, system, x0, N)
            % simplot(ctrl, system, x0, Options)
            if ~isa(arg3, 'double'),
                error('SIMPLOT: third input must be vector of initial conditions.');
            end
            if isa(arg2, 'function_handle'),
                Options.sysHandle = arg2;
            else
                Options.sysStruct = arg2;
            end
            if isa(arg4, 'double'),
                % simplot(ctrl, system, x0, N)
                if any(size(arg4)>1) | round(arg4)~=arg4,
                    error('SIMPLOT: number of simulation steps must be a positive integer.');
                end
                [X, U, Y, cost, feasible] = mpt_plotTimeTrajectory(ctrl, arg3, arg4, Options);
                
            elseif isa(arg4, 'struct'),
                % simplot(ctrl, system, x0, Options)
                Options = arg4;
                if isa(arg2, 'function_handle'),
                    Options.sysHandle = arg2;
                else
                    Options.sysStruct = arg2;
                end
                [X, U, Y, cost, feasible] = mpt_plotTimeTrajectory(ctrl, arg3, [], Options);
                
            else
                error('SIMPLOT: wrong type of input arguments.');
            end
        else
            error('SIMPLOT: wrong type of input arguments.');
        end

    case 5,
        % simplot(ctrl, system, x0, N, Options)
        if ~isa(arg5, 'struct'),
            error('SIMPLOT: fifth input must be an options structure.');
        end
        if ~isa(arg4, 'double'),
            error('SIMPLOT: fourth input must be a positive integer.');
        end
        if any(size(arg4)>1) | round(arg4)~=arg4,
            error('SIMPLOT: fourth input must be a positive integer.');
        end
        if ~isa(arg3, 'double'),
            error('SIMPLOT: vector of initial conditions must be given.');
        end
        Options = arg5;
        if isa(arg2, 'function_handle'),
            Options.sysHandle = arg2;
        else
            Options.sysStruct = arg2;
        end
        [X, U, Y, cost, feasible] = mpt_plotTimeTrajectory(ctrl, arg3, arg4, Options);

    otherwise
        error('SIMPLOT: wrong number of input arguments.');
end

if nargout < 1,
    clear X U Y cost feasible
end
