function [xn,y,z,d,feasible]=mpt_mldsim(S, x0, u)
%MPT_MLDSIM Simulates an MLD system for one time step
%
% [xn,y,z,d,feasible]=mpt_mldsim(S, x0, u)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Simulates an MLD system for one time step.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% S        - structure contianing matrices of MLD system
% x0       - initial state
% u        - current input
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% xn       - next state
% y        - current output
% d        - current values of boolean variables
% z        - current values of auxiliary continuous variables
% feasible - 1 if a feasible solution exists, 0 or -1 otherwise
%
% see also MPT_SIMSYS

% Copyright is with the following author(s):
%
%(C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch

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

error(nargchk(3,3,nargin));

if ~isstruct(mptOptions)
    mpt_error;
end

x0 = x0(:);
u = u(:);

nx = S.nx;
nu = S.nu;
ny = S.ny;
nd = S.nd;
nz = S.nz;
ne = S.ne;

if length(x0) ~= nx,
    error('MPT_MLDSIM: Wrong dimension of x0.');
end
if length(u) ~= nu,
    error('MPT_MLDSIM: Wrong dimension of u.');
end

xn = [];
y = [];
d = [];
z = [];

if ne==0,
    % no inequalities, compute xn and y directly
    xn = S.A * x0 + S.B1 * u;
    y = S.C * x0 + S.D1 * u;
    if isfield(S, 'B5'),
        % add affine term if available
        xn = xn + S.B5;
    end
    if isfield(S, 'D5'),
        % add affine term if available
        y = y + S.D5;
    end
    feasible = 1;
    
elseif isempty(S.E2) & isempty(S.E3)
    % constraints, but no "z" and "d" variables
    B = S.E1*u + S.E4*x0 + S.E5;
    A = zeros(size(B));
    if any(B < A)
        disp('Warning: Constraints lead to infeasible or unbounded solution.');
        xn = zeros(nx, 1);
        y = zeros(ny, 1);
        d = zeros(nd, 1);
        z = zeros(nz, 1);
        feasible = 0;
    else
        xn = S.A*x0 + S.B1*u;
        y = S.C*x0 + S.D1*u;
        d = zeros(nd, 1);
        z = zeros(nz, 1);
        if isfield(S, 'B5')
            % add affine term if available
            xn = xn + S.B5;
        end
        if isfield(S, 'D5'),
            % add affine term if available
            y = y + S.D5;
        end
        feasible = 1;
    end
    
else
    % general case, solve feasiblity MILP to get "d" and "z"
    
    % no cost, just solve feasibility problem
    f = zeros(nd + nz, 1);
    
    % matrices containing inequality constraints
    A = [S.E2 S.E3];
    B = [S.E1*u + S.E4*x0 + S.E5];
    
    % type of variables:
    % first "nd" variables are boolean, the rest is continuous
    vartype = repmat('C', nd+nz, 1);
    vartype(1:nd) = 'B';
    
    % lower and upper bounds of "d" and "z" variables
    lb = [zeros(nd, 1); S.zl];
    ub = [ones(nd, 1); S.zu];
    
    % solve the feasibility MILP
    [dzmin,fmin,how,feasible]=mpt_solveMILP(f,A,B,[],[],lb,ub,vartype);

    if feasible==1,
        % problem is feasible
        d = dzmin(1:nd);
        z = dzmin(nd+1:end);
        
        % compute state update and output
        xn = S.A*x0 + S.B1*u + S.B2*d(:) + S.B3*z(:);
        y = S.C*x0 + S.D1*u + S.D2*d(:) + S.D3*z(:);
        if isfield(S, 'B5'),
            % add affine term if available
            xn = xn + S.B5;
        end
        if isfield(S, 'D5'),
            % add affine term if available
            y = y + S.D5;
        end
    else
        disp('Warning: Constraints lead to infeasible or unbounded solution.');
        xn = zeros(nx, 1);
        y = zeros(ny, 1);
        d = zeros(nd, 1);
        z = zeros(nz, 1);
    end
end
