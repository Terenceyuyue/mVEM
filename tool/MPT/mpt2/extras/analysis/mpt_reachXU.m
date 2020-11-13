function [Xn, V] = mpt_reachXU(X0, U0, A, B, f, Options)
%MPT_REACHXU One-step reachability computation for an affine system
%
% Xn = mpt_reachXU(X0, U0, A, B);
% Xn = mpt_reachXU(X0, U0, A, B, f);
% [Xn, V] = mpt_reachXU(X0, U0, A, B);
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes set Xn such that
%
%  Xn = { xn | xn = A*x + B*u + f, x \in X0, u \in U0}
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% X0       - set of initial conditions on states
% U0       - set of initial conditions on inputs
% A, B, f  - state update matrices
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Xn       - polytope Xn (may be not full dimensional)
% V        - vertices of Xn
%
% see also POLYTOPE/RANGE, MPT_REACHSETS
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

global mptOptions

error(nargchk(4,6,nargin));

if nargin<6,
   Options = [];
end

if ~isfield(Options, 'verbose')
    Options.verbose = mptOptions.verbose;
end
if ~isfield(Options, 'userange'),
    Options.userange = 0;
end

if ~isa(X0, 'polytope') | ~isa(U0, 'polytope')
    error('First two input arguments must be polytope objects!');
end

if size(A, 2) ~= dimension(X0),
    error('A matrix must be of the same dimension as polytope X0!');
end

if size(B, 1) ~= size(A, 1),
    error('B matrix must have the same number of rows as A matrix!');
elseif size(B, 2) ~= dimension(U0),
    error('B matrix must be of the same dimension as polytope U0!');
end

if nargin < 5,
    f = zeros(size(A,1), 1);
end

if size(f, 2) ~= 1,
    error('Fifth input arguments must be a column vector!');
elseif size(f, 1) ~= size(A, 1),
    error('Affine "f" term must have the same number of rows as the A matrix!');
end

% create a polytope in X, U space
XU = X0 * U0;

if Options.userange,
    if size(A, 1) ~= size(A, 2),
        error('A matrix must be square!');
    end
    nu = size(B, 2);
    nx = size(A, 1);
    Acl = [A B; zeros(nu, nx) eye(nu)];
    fcl = [f; zeros(nu,1)];
    if det(Acl) < mptOptions.abs_tol,
        Options.userange = 0;
        disp('Mapping is lower-dimensional, switching to vertex based method...');
        [Xn, V] = mpt_reachXU(X0, U0, A, B, f, Options);
        return
    end
    R = range(XU, Acl, fcl);
    if isfulldim(R),
        Xn = projection(R, 1:nx);
    else
        Xn = R;
    end
    V = [];
    return
end

% compute extreme points of the polytope
Exu = extreme(XU);

if isempty(Exu)
    error('Polytope in extended space is not fully dimensional!');
end

Exn = [A B] * Exu';
Pxn = polytope(Exn');

[xcheb, rcheb] = chebyball(Pxn);

if ~isfulldim(Pxn) | rcheb == 1e9,
    if Options.verbose > -1,
        disp('Polytope is not fully dimensional!');
    end
    Pxn = polytope;
end

if nargin>=5,
    if isfulldim(Pxn),
        Pxn = Pxn + f;
    end
    Exn = Exn + repmat(f, 1, size(Exn, 2));
end

Xn = Pxn;
V = Exn';