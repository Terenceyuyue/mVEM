% Multi-Parametric Toolbox Control routines
%   mpt_control            - Main control routine. Computes explicit controller for a given problem
%
%   mpt_optControl         - CFTOC of LTI systems
%   mpt_optControlPWA      - CFTOC of PWA systems with linear cost
%   mpt_optQuadCtrl        - CFTOC of PWA systems with quadratic cost
%   mpt_optInfControl      - Solves the infinite-time constrained optimal control problem for LTI systems
%   mpt_optInfControlPWA   - Solves the infinite-time constrained optimal control problem for PWA systems
%   mpt_iterative          - Computes a time-optimal or low-complexity explicit controller for LTI systems
%   mpt_iterativePWA       - Computes a time-optimal or low-complexity explicit controller for PWA systems
%   mpt_oneStepCtrl        - Computes low complexity controller for LTI systems
%   mpt_simplexContr       - Computes a piecewise affine feedback law defined over simplices
%   mpt_optBoolCtrl        - Computes optimal controller for systems with discrete inputs
%   mpt_optMixedCtrl       - Computes optimal controller for systems discrete and continuous inputs
%   mpt_boolMinTime        - Computes minimum time controller for systems with discrete inputs
%   mpt_mixedMinTime       - Computes minimum time controller for systems with discrete and continuous inputs
%   mpt_solveMPC           - Solves the on line optimization MPC problem
%
% Authors: Michal Kvasnica, Pascal Grieder, Mato Baotic
% kvasnica@control.ee.ethz.ch, grieder@control.ee.ethz.ch, baotic@control.ee.ethz.ch
%
% Copyright (C) 2003-2006 Michal Kvasnica, Pascal Grieder, Mato Baotic
%
% For support, write to: mpt@control.ee.ethz.ch
