% YALMIP
% Version 3 (R2009b) 11-June-2010
%
% Information
%
% Variables.
%   sdpvar       - Create SDPVAR variable 
%   intvar       - Create integer SDPVAR variable
%   binvar       - Create binary SDPVAR variable
%   blkvar       - Create container for block variable
%   double       - Convert SDPVAR object to double (i.e. extract solution)
%   is           - Check property of SDPVAR variable
%   setsdpvar    - Initialize the "double" value of an SDPVAR object
%   replace      - Replace free variables in an SDPVAR object with constant
%   sdisplay     - Tries to display an SDPVAR symbolically
%   monolist     - Creates multi-variate monomials up to specified degree
%   linearize    - Linearize SDPVAR object p(x) at the point double(x)
%   jacobian     - Calculate Jacobian for SDPVAR object
%   hessian      - Calculate Jacobian for SDPVAR object
%   coefficients - Extract coefficients and monomials from polynomials
%   unblkdiag    - Detect and extract blocks in SDP constraints
%
% Inequalites and constraints.
%   set          - Create a SET object (set of constraints). For on-line help: help sdpvar/set
%   sos          - Create sum-of-squares constraint
%   cone         - Second order cone constraint
%   rcone        - Rotated Lorentz cone constraint
%   integer      - Constrain variables to be integer
%   binary       - Constrain variables to be binary
%   unblkdiag    - Extracts diagonal blocks
%   dissect      - Dissect SDP constraint
%   checkset     - Checks constraint violations for current solution
%   dual         - Extract dual variable related to a constraint SET
%   sosd         - Extract SOS decomposition
%   double       - Convert SET object to double
%   linearize    - Linearize all constraints
%   is           - Checks property of constraint
%   plot         - Plot feasible set. For on-line help: help lmi/plot
%   hull         - construct convex hull of SET objects: help lmi/hull
%
% Optimization related.
%   sdpsettings  - Create an options structure
%   solvesdp     - Computes solution to optimization problem
%   solvesos     - Sum of squares decomposition of polynomial
%   solvemoment  - Lasserre's moment relaxation for polynomial programming
%   solvemp      - Solve multi-parametric program
%
% Auxillary
%   export       - Export YALMIP model to various solver formats
%   yalmip       - Various administrative stuff
%   savesdpafile - Saves a problem definition in the SDPA format
%   yalmipdemo   - Brief tutorial and examples.
%   yalmiptest   - Runs a number of test problems.
%

% Author Johan Löfberg
% $Id: Contents.m,v 1.53 2009-10-14 09:10:35 joloef Exp $   
