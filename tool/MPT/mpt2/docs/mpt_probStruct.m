%MPT_PROBSTRUCT Description of the problem structure
%
% Problem structure (probStruct) is a structure which states an optimization
% problem to be solved by MPT. 
%
% ---------------------------------------------------------------------------
% ONE AND INFINITY NORM PROBLEMS
% ---------------------------------------------------------------------------
%
% The optimal control problem for linear performance index is given by:
%
%                         N-1
%                         ___
%      ||        ||       \   ||      ||      ||      ||
% min  ||P_N x(N)||    +  /   ||Q x(k)||    + ||R u(k)||     
%  U   ||        || p     --- ||      || p    ||      || p
%                         k=0
%
% subject to:
%    x(k+1) = fdyn( x(k), u(k), w(k) )
%               
%     x(N) \in Tset
%     
%
% where:
%   U   - vector of manipulated variables over which the optimization is
%         performed
%   N   - prediction horizon
%   p   - linear norm, can be 1 or Inf for 1- and Infinity-norm, respectively
%   Q   - weighting matrix on the states
%   R   - weighting matrix on the manipulated variables
%   P_N - weight imposed on the terminal state
%   Tset         - terminal set
%
%   the function fdyn( x(k), u(k), w(k) ) is the state-update function and
%   is different for LTI and for PWA systems (see help mpt_sysStruct for more
%   details on system dynamics).
%
%
% ---------------------------------------------------------------------------
% TWO NORM PROBLEMS
% ---------------------------------------------------------------------------
%
% In case of quadratic performance index, the optimal control problem takes
% the following form: 
%
%                        N-1
%                        ___
%          T             \        T              T
% min  x(N)  P_N x(N) +  /    x(k)  Q x(k) + u(k)  R u(k)
%  U                     ---
%                        k=0
%
% subject to the same constraints as mentioned above.
%
% If the problem is formulated for a fixed prediction horizon N, we refer to it 
% as to Constrained Finite Time Optimal Control Problem (CFTOC). On the other
% hand, if N is infinity, the Constrained Infinite Time Optimal Control Problem
% (CITOC) is formulated. Objective of the optimization is to choose the
% manipulated variables such that the performance index is minimized.
%
% ---------------------------------------------------------------------------
% LEVEL OF OPTIMALITY
% ---------------------------------------------------------------------------
%
% MPT provides different control strategies. The cost-optimal solution leads to
% a control law which minimizes a given performance index. This strategy is
% enforced by
%
%   probStruct.subopt_lev = 0
%
% The cost optimal solution for PWA systems is currently supported only for
% linear performance index (i.e. probStruct.norm = 1 or Inf)
%
%
% An another possibility is to use the time-optimal solution, i.e. the control
% law will push a given state to the origin (or to an invariant set around the
% origin) as fast as possible. This strategy usually leads to more simple
% control laws (i.e. less controller regions are generated). This approach is
% enforced by
%
%   probStruct.subopt_lev = 1
%
%
% The last option is to use a low-complexity control scheme. For LTI systems,
% this approach aims at constructing a one-step solution and subsequent PWQ
% Lyapunov function computation to testify stability properties. In case of PWA
% systems, the approach aims at reducing the number of switches in the
% closed-loop system. If you want to use this kind of solution, use:
%
%   probStruct.subopt_lev = 2
%
% ---------------------------------------------------------------------------
% NOTATION
% ---------------------------------------------------------------------------
%
% In order to specify which problem the user wants to solve, the following
% fields of the problem structure probStruct have to be provided: 
%
%   probStruct.N    - prediction horizon
%   probStruct.Q    - weights on states
%   probStruct.R    - weights on inputs
%   probStruct.norm - can be either 1 or Inf for linear performance index, or 2
%                     for quadratic cost objective
%
% By default, probStruct.Q and probStruct.R are assumed to be matrices. It is
% also possible to use time-varying weights in the problem formulation. In such
% case, define probStruct.Q and probStruct.R as cell arrays with probStruct.N
% elements.
% 
% In addition, several optional fields can be given:
% 
%   probStruct.Qy       - cost on outputs. Mandatory for output regulation and
%     output tracking problems.
%  
%   probStruct.Qd       - penalty on boolean variables (only for control of MLD
%                         systems) 
%
%   probStruct.Qz       - penalty on auxiliary variables (only for control of
%                         MLD systems) 
%
%   probStruct.y0bounds - a true/false (1/0) flag denoting
%     whether or not to impose constraints also on the initial system output
%     (default is 1 = enabled)
%
%   probStruct.tracking - {0/1/2} flag
%     0 - no tracking, resulting controller is a state regulator which drives
%         all system states (or outputs, if probStruct.Qy is given) towards
%         origin 
%     1 - tracking with deltaU formulation - controller will drive system states
%         (or outputs, if probStruct.Qy is given) to given reference. The
%         optimization is performed over the difference of manipulated variables
%         (u(k) - u(k-1)), which involves extension of the state vector by "nu"
%         additional states where "nu" is number of system inputs
%     2 - tracking without deltaU formulation - same as probStruct.tracking=1
%         with the exception that optimization is performed over u(k), i.e. no
%         deltaU formulation is used and no state vector extension is needed.
%         Note, however, that offset-free tracking cannot be guaranteed with
%         this setting.
%     (default is 0 = no tracking)
%
%   probStruct.yref     - instead of driving a state to zero, it is
%     possible to reformulate the control problem and rather force the output
%     to zero. to ensure this task, define probStruct.Qy which penalizes
%     difference of the actual output and the given reference. you will need to
%     define "probStruct.Qy" in order to use this option.
%   
%   probStruct.xref     - by default the toolbox designs a controller which
%     forces the state vector to convert to the origin. If you want to track
%     some a-priori given reference point, provide the reference state in this
%     variable. probStruct.tracking has to be 0 (zero) to use this option!
%
%   probStruct.uref     - Similarly a reference point for the manipulated
%     variable (i.e. the equilibrium $u$ for state probStruct.xref can be
%     specified here. If it is not given, it is assumed to be zero. 
%
%   probStruct.dref     - reference for boolean variables (only for control of
%     MLD systems). Can be used together with probStruct.Qd
%
%   probStruct.zref     - reference on auxiliary variables (only for control of
%     MLD systems). Can be used together with probStruct.Qz
%
%   probStruct.P_N      - weight on the terminal state. If not specified, it is
%     assumed to be solution of the Ricatti equation, or P_N = Q for linear cost
%
%   probStruct.Tset     - a polytope object describing the terminal set. If not
%     provided and probStruct.norm is 2, the LQR set around the origin will be
%     calculated automatically to guarantee stability properties.
%
%   probStruct.Tconstraint - an integer (0, 1, 2) denoting which stability
%     constraint to apply. 0 - no terminal constraint, 1 - use LQR terminal set
%     2 - use user-provided terminal set constraint). Note that if
%     probStruct.Tset is given, Tconstraint will be set to 2 automatically.
%
%   probStruct.xN       - terminal state constraint. if given, a hard constraint
%     on the last predicted state in the form "x(N)==probStruct.xN" will be
%     imposed. NOTE! only available for probStruct.subopt_lev = 0.
%
%   probStruct.Nc       - control horizon. Specifies number of free control
%     moves in the optimization problem. If "probStruct.Nc" < "probStruct.N",
%     all control moves between Nc and N will be constrained to be identical.
%     E.g. if Nc=2 and N=5, then u_0 and u_1 will be considered as free control
%     moves, and u_2==u_3==u_4 will be enforced.
%
%   probStruct.feedback - boolean variable, if set to 1, the problem is
%     augmented such that U = K x + c   where K is a state-feedback gain
%     (typically a LQR controller) and the optimization aims to identify the
%     proper offset "c".
%     (default is 0 = no)
%
%   probStruct.FBgain   - if the former option is activated, a specific
%     state-feedback gain matric K can be provided (otherwise a LQR controller
%     will be computed automatically)
%
% ---------------------------------------------------------------------------
% SOFT CONSTRAINTS
% ---------------------------------------------------------------------------
%
% Since MPT 2.6 it is possible to denote certain constraints as soft. This means
% that the respective constraint can be violated, but such a violation is
% penalized. To soften certain constraints, it is necessary to define the
% penalty on violation of such constraints:
%
%   probStruct.Sx   - if given as a "nx" x "nx" matrix, all state constraints
%                     will be treated as soft constraints, and violation will be
%                     penalized by the value of this field.
%   probStruct.Su   - if given as a "nu" x "nu" matrix, all input constraints
%                     will be treated as soft constraints, and violation will be
%                     penalized by the value of this field.
%   probStruct.Sy   - if given as a "ny" x "ny" matrix, all output constraints
%                     will be treated as soft constraints, and violation will be
%                     penalized by the value of this field.
% 
% In addition, one can also specify the maximum value by which a given
% constraint can be exceeded:
%
%   probStruct.sxmax - must be given as a "nx" x 1 vector, where each element
%                      defines the maximum admissible violation of each state
%                      constraints
%   probStruct.sumax - must be given as a "nu" x 1 vector, where each element
%                      defines the maximum admissible violation of each input
%                      constraints
%   probStruct.symax - must be given as a "ny" x 1 vector, where each element
%                      defines the maximum admissible violation of each output
%                      constraints
%
% The afforementioned fields also allow to tell that only a subset of state,
% input, or output constraint should be treated as soft constraints, while the
% rest of them remain hard. Say, for instance, that we have a system with 2
% states and we want to soften only the second state constraint. Then we would
% write:
%
%   probStruct.Sx = diag([1 1000])
%   probStruct.sxmax = [0; 10]
%
% Here probStruct.sxmax(1)=0, which tells MPT that the first constraint should
% be treated as a hard constraint, while we are allowed to exceed the second
% constraints at most by the value of 10, while every such violation will be
% penalized by the value of 1000.


% Copyright is with the following author(s):
%
%(C) 2003-2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%              kvasnica@control.ee.ethz.ch
%

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


error('This script is not executable. Type "edit mpt_probStruct" to see more information.');

%% Template for problem structure
% * Read help description of this file for more details
% * Remove all fields marked as optional if you don't need them.
% * Certain combinations are not allowed (e.g. tracking & time-optimal control) 


%% Prediction horizon
probStruct.N = [];

% Use an integer (1,2,...) for finite horizon problem
% Use Inf for infinite-time and time-optimal problems


%% Type of performance criterion
probStruct.norm = 1/2/Inf
% 1/Inf correspond to linear performance index (no need of QP solver)
% 2 indicates quadratic performance objective


%% Weighting matrices in the optimization objective
probStruct.Q = [];
probStruct.R = [];


%% Time-varying weighting matrices
probStruct.Q = { Q1, Q2, ..., QN};
probStruct.R = { R1, R2, ..., RN};
% Length of the above fields must be consistent with the prediction horizon


%% Type of control problem
probStruct.subopt_lev = 0;   % for cost-optimal solutions
probStruct.subopt_lev = 1;   % for time-optimal solutions
probStruct.subopt_lev = 2;   % for low complexity strategies


%% Output regulation
% Instead of driving a state to zero, it is possible to reformulate the control
% problem and rather force the output to zero. To ensure this task, define
% probStruct.Qy which penalizes difference of the actual output and the given
% reference
% (optional; assuming state regulation if not provided)
probStruct.Qy = [];


%% State tracking with arbitrary reference 
% (solved by augmentation of the state vector, may result in long computation.
% Consult the "MPT in 15 minutes"  section of MPT manual for more details)

% (optional; assuming tracking=0 if not given)
probStruct.tracking = 0;  % set to 1 to enable tracking


%% State tracking with a-priori known reference
% (optional)
probStruct.xref = [];
probStruct.uref = [];


%% Output tracking with a-priori known reference
% (optional)
% it is furthermore necessary to define probStruct.Qy for this type of tracking!
probStruct.yref = [];


%% Type of final state constraint
% (optional; assuming probStruct.Tconstraint=1 if argument not defined)

probStruct.Tconstraint = 0; % use final state constraint
probStruct.Tconstraint = 1; % final state forced to be in LQR set (stability guarantee)
probStruct.Tconstraint = 2; % final state forced to be in user defined target set


%% User defined target set
% (optional; must be given if Tconstraint=2)
probStruct.Tset = polytope(H,K);  % final state forced to enter set H x <= K at the end of prediction horizon


%% User defined penalty on final state
% (optional; must be given if Tconstraint=0)
probStruct.P_N = [];


