% Multi-Parametric Toolbox Simulink Library
%
% Simulink library
%   mpt_sim          - MPT Simulink library
%
% C-code functions
%   mpt_getInput.c       - C implementation of explicit controllers
%   mex_getInput         - MEX interface to mpt_getInput
%   mpt_exportc          - Exports a given explicit controller to C-code
%   mpt_getInput_sfunc.c - Pure C S-function for simulations of explicit controllers
%
% Internal functions:
%   mpt_getSfuncParam    - Prepares parameters for explicit controller S-function
%   mpt_simInput         - Evaluation of on-line controllers in Simulink blocks
%   mpt_simSys_sfunc     - S-function to simulate sysStruct system in Simulink
%
% Authors: Michal Kvasnica, Pascal Grieder, Mato Baotic
% kvasnica@control.ee.ethz.ch, grieder@control.ee.ethz.ch, baotic@control.ee.ethz.ch
%
% Copyright (C) 2003-2006 Michal Kvasnica, Pascal Grieder, Mato Baotic
%
% For support, write to: mpt@control.ee.ethz.ch
