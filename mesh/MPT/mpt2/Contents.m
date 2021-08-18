% Multi-Parametric Toolbox
% Version 2.6.3 (R2008a) 14-Jan-2008
%
%
% Information
%   README.TXT             - Installation tips etc.
%   HISTORY.TXT            - Version history.
%   BUGS.TXT               - Known bugs and limitations
%
% Graphical User Interface
%   mpt_studio             - MPT Graphical User Interface
%   mpt_setup              - Main MPT configuration routine
%
% Control routines
%   mpt_control            - Computes an explicit or an on-line MPC controller
%   mpt_ownmpc             - "Design your own MPC problem" function
%
% Simulink library and code generation
%   mpt_sim                - MPT Simulink library
%   mpt_exportc            - Exports a given explicit controller to C-code
%
% Solvers
%   mpt_mplp               - Multi-Parametric LP solver
%   mpt_mpqp               - Multi-Parametric QP solver
%   mpt_mpmilp             - Multi-Parametric Mixed-Integer LP solver
%   mpt_mpmiqp             - Multi-Parametric Mixed-Integer QP solver
%   mpt_solveLP            - Interface to various LP solvers
%   mpt_solveQP            - Interface to various QP solvers
%   mpt_solveMILP          - Interface to various MILP solvers
%   mpt_solveMIQP          - Interface to various MIQP solvers
%
% Analysis and post-processing
%   mpt_invariantSet       - Computes invariant set for a given explicit controller
%   mpt_reachSets          - Computes sets of reachable states for a given system / controller
%   mpt_verify             - Verification of hybrid systems
%   mpt_simplify           - Simplifies a given explicit controller
%   mpt_searchTree         - Generate a binary search tree for fast region identification
%
% Setting up the optimization problem
%   mpt_constructMatrices  - Constructs matrices for the finite time constrained optimal control problem
%   mpt_blockingMatrices   - Constructs matrices for the CFTOC problem for move blocking strategies
%
% Lyapunov analysis routines
%   mpt_lyapunov           - Interface for Lyapunov-type function computation
%
% Application of control law
%   mpt_getInput           - For a given state, extracts the optimal output from controller structure
%   mpt_computeTrajectory  - Calculates time evolution of state trajectories subject to control
%   mpt_simSys             - Simulates evolution of a given system
%   mpt_plotTimeTrajectory - Plots open-loop or closed-loop trajectories
%   
% Visualization functions
%   mpt_plotPartition      - Plots a polyhedral partition obtained by mpt_control
%   mpt_plotTrajectory     - Graphical interface for piloting trajectories for LTI and PWA systems subject to control
%   mpt_plotPWA            - Plots a PWA function defined over a given polyhedral partition
%   mpt_plotPWQ            - Plots a PWQ function defined over polyhedral partition
%   mpt_plotJ              - Plots value function associated to a given controller
%   mpt_plotU              - For a given explicit controller, plots value of the control action
%   mpt_plotArrangement    - Plots hyperplane arrangement of a polytope in H-representation
%   mpt_plotSysStruct      - Plots the system partitions of a PWA system
%
% Toolbox interface methods
%   mpt_init               - Initializes the MPT toolbox
%   mpt_update             - Checks for updates of the Multi-Parametric Toolbox
%   mpt_version            - Returns version of MPT
%   mpt_options            - Read / modify solver settings for MPT
%
% Geometry functions
%   mpt_getInnerEllipsoid  - Computes the largest ellipsoid inscribed in a polytope
%   mpt_getOutterEllipsoid - Computes the smallest ellipsoid which covers the polytope P
%   mpt_plotellip          - Plots a given ellipsoid
%   mpt_delaunay           - Computes the delaunay triangulation of a polytope 
%   mpt_voronoi            - Computes the voronoi diagram via mpLP
%
% Others
%   hull                   - Converts vertices into an H-representation polytope
%   mousepoly              - allows to "draw" a polytope in 2D by mouse
%   unitbox                - Creates a hypercube of given dimension centered at origin
%   mpt_removeOverlaps     - Removes overlaps from (a set of) polyhedral partitions with associated linear cost
%   mpt_randLTISys         - Generates random LTI systems
%   mpt_randPWASys         - Generates random PWA systems
%
% Demos and Examples
%   runExample             - Demonstrates MPT control routines
%   runMoveBlocking        - Demonstrates MPT Move Blocking routines
%   mpt_demo1              - Explains basic manipulation with polytopes
%   mpt_demo2              - A tour through visualization capabilities of the toolbox
%   mpt_demo3              - Explains control routines for LTI systems
%   mpt_demo4              - Explains control routines for PWA systems
%   mpt_demo5              - Explains tracking functionality
%   mpt_demo6              - Illustrates implementation and visualization of the control law
%   reachdemo1, reachdemo2 - Reachability analysis demos
%   verifdemo1, verifdemo2 - Verification demos
%   ballandplatedemo       - "Ball and plate" demo
%   turbocardemo           - "Car with turbo" demo
%   twotanksdemo           - "Two tanks" demo
%
%   mpt_sysStruct          - Help file for the system structure sysStruct
%   mpt_probStruct         - Description of the problem structure
%
% Polytope library
%   see help mpt/polytope
%
% Authors: Michal Kvasnica, Pascal Grieder, Mato Baotic
% kvasnica@control.ee.ethz.ch, grieder@control.ee.ethz.ch, baotic@control.ee.ethz.ch
%
% Copyright (C) 2003-2006 Michal Kvasnica, Pascal Grieder, Mato Baotic
%
% For support, write to: mpt@control.ee.ethz.ch
