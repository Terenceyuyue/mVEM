% Multi-Parametric Toolbox Invariant set computation functions
%   mpt_reachSets          - Computes sets of reachable states for a given system / controller
%   mpt_maxCtrlSet         - Computes maximal (robust) control invariant set of
%                            maximal (robust) attractive set
%   mpt_infset             - Calculates the maximal positively invariant set for an LTI system
%   mpt_infsetPWA          - Computes (robust) positive invariant subset for PWA systems
%   mpt_computePWATset     - Computes a stabilizing control nvariant set (+ controllers) around the origin
%   mpt_getReachSubset     - Computes a subset of P which enters Pfin
%   mpt_getStabFeedback    - Computes a stabilizing feedback law for a PWA system
%
% Lyapunov analysis routines
%   mpt_getQuadLyapFct     - Computes common Lyapunov function for PWA system
%   mpt_getPWALyapFct      - Calculates Piecewise-Affine Lyapunov function
%   mpt_getPWQLyapFct      - Calculates Piecewise-Quadratic Lyapunov function
%   mpt_getPWQLyapFct      - Calculates Piecewise-Polynomial Lyapunov function
%   mpt_checkLyapFct       - Checks if decay rate of a PWA/PWQ Lyapunov function is always negative
%
% Authors: Michal Kvasnica, Pascal Grieder, Mato Baotic
% kvasnica@control.ee.ethz.ch, grieder@control.ee.ethz.ch, baotic@control.ee.ethz.ch
%
% Copyright (C) 2003-2006 Michal Kvasnica, Pascal Grieder, Mato Baotic
%
% For support, write to: mpt@control.ee.ethz.ch
