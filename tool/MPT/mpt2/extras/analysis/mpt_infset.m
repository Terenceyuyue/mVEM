function [Oinf,tstar,fd,isemptypoly] = mpt_infset(A,X,tmax,Pnoise,Options)
%MPT_INFSET Calculates the maximal positively invariant set for an LTI system
%
% [Oinf,tstar,fd,isemptypoly] = mpt_infset(A,X,tmax,Pnoise,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Calculates the maximal positively invariant set for an autonomous discrete-time
% LTI system
% Takes into account polytopic and additive system uncertainty, if defined
% in "sysStruct" / "Pnoise", respectively.
%
% ---------------------------------------------------------------------------
% INPUT
% --------------------------------------------------------------------------- 
% A                - The A matrix of the discrete-time LTI system x{k+1} = Ax{k}.
% X                - State constraints given as a polytope (x \in X)
% tmax             - Maximum number of iterations allowed.
% Options.lpsolver - Solver for LPs (see help mpt_solveLP for details)
% Options.verbose  - level of verbosity (see help mpt_init for details)
% Options.abs_tol  - absolute tolerance
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
%   OUTPUT:
% ---------------------------------------------------------------------------
% Oinf         - Maximal positively invariant set contained in X, i.e.
%                  Oinf = {x_k \in X: x_{k+1} = Ax_k \in X}.
% tstar        - Determinedness index.
% fd           - 1 if Oinf is finitely determined (tstar <= tmax).
%                0 if Oinf tstar > tmax.
% isemptypoly  - 1 if resulting polyhedron is empty; 0 otherwise
%
% ---------------------------------------------------------------------------
%   Literature:
% ---------------------------------------------------------------------------
% "Theory and Computation of Disturbance Invariant Sets for Discrete-Time Linear Systems"
% I. Kolmanovsky and E. G. Gilbert, Mathematical Problems in Egineering, vol (4), 1998, 
% pages 317-367
%
% AND
%
% "Linear systems with state and control constraints: the theory and
% applications of maximal output admissible sets", E. G. Gilbert and K. Tin Tan,
% IEEE Transcations on Automatic Control, 1991, vol 36, number 9, pages 1008--1020
%
%
% see also MPT_INFSETPWA

% Copyright is with the following author(s):
%
% (C) 2003-2007 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
%
%  This function is based on a script by Eric Kerrigan. Thanks Eric, for letting us
%  use you code:)

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

error(nargchk(3,5,nargin));

if nargin < 5
    Options = [];
end
if nargin < 4
    Pnoise = polytope;
end

Options = mpt_defaultOptions(Options, ...
    'maxIter', tmax, ...
    'verbose', -1 );

if ~iscell(A)
    A = { A };
end
ndyn = length(A);
if ndyn == 1 & length(X) > 1
    % handle special case with a single dynamics but the target is a
    % polytope array
    ndyn = length(X);
    Anew = cell(1, ndyn);
    [Anew{:}] = deal(A{1});
    A = Anew;
end

nx = size(A{1}, 1);
f = cell(1, ndyn);
[f{:}] = deal(zeros(nx, 1));

% mpt_infsetPWA() is more numerically robust, uses more heuristics and
% handles more general cases (such as "A" being a cell array or "X" being a
% polytope array)
Oinf = mpt_infsetPWA(X, A, f, Pnoise, Options);
fd = 1;
tstar = 1;
isemptypoly = ~isfulldim(Oinf);
