% Multi-Parametric Toolbox Auxiliary Functions
%
%   identifyRegion        - Plots the number of the polytopes into the current 
%   mousepoly             - Allows to specify polytope by mouse-clicks
%   mpt_addTset           - Adds terminal set constraint to Matrices
%   mpt_allcombs          - All combinations of discrete-valued inputs
%   mpt_checkLyapFct      - Checks if decay rate of PWQ(PWA) Lyapunov function is always negative
%   mpt_error             - Function called if MPT toolbox is not initialized
%   mpt_evalSystem        - Extracts data from sysStruct and probStruct structures
%   mpt_feasibleStates    - Returns equidistantly spaced data points in feasible set
%   mpt_getCommonLyapFct  - Computes common Lyapunov function for PWA system
%   mpt_getFullRankSubset - Removes rows from matrix Gt until it has full row rank
%   mpt_greedyMerging     - Greedy merging of polyhedra
%   mpt_isValidCS         - Checks if input argument is a valid controller structure
%   mpt_iscombequal       - Are two vectors combinatorially equal
%   mpt_lrs               - Matlab implementation of the LRS algorithm
%   mpt_lti2pwa           - Converts an LTI system to a PWA system
%   mpt_mergeCS           - Merges a cell array of ctrlStruct structures
%   mpt_mplp_ver1         - Explicitly solves the given linear program (LP)
%   mpt_mplp_ver2         - Explicitly solves the given linear program (LP)
%   mpt_mplp_ver3         - Explicitly solves the given linear program (LP)
%   mpt_mplp_ver4         - Explicitly solves the given linear program (LP)
%   mpt_mplp_ver5         - Explicitly solves the given linear program (LP)
%   mpt_optControlPWAold  - Solves the CFTOC problem for a given PWA system
%   mpt_patch2eps         - Exporting 3-dimensional transparent patches to eps for use in LaTeX
%   mpt_prepareDU         - Extends system and problem matrices to deal with deltaU constraints
%   mpt_prepareTracking   - Extends system and problem matrices to deal with tracking
%   mpt_qphess            - Computes the Hessian
%   mpt_solveLPi          - Interface to various LP solvers (version without errorchecks)
%   mpt_solverInfo        - returns information about a given solver
%   mpt_sysStructInfo     - Returns information about system structure
%   mpt_verifyProbStruct  - Verifies the probStruct structure
%   mpt_verifySysProb     - Verifies system and problem structures
%   mpt_verifySysStruct   - Verifies the sysStruct structure
%   plotc                 - Shortcut function which plots polyhedral partition of a given controller structure
%   unitbox               - Creates a unit box centered at origin
%
% Authors: Michal Kvasnica, Pascal Grieder, Mato Baotic
% kvasnica@control.ee.ethz.ch, grieder@control.ee.ethz.ch, baotic@control.ee.ethz.ch
%
% Copyright (C) 2003-2006 Michal Kvasnica, Pascal Grieder, Mato Baotic
%
% For support, write to: mpt@control.ee.ethz.ch
