function S = addMLDconstraints(S, sysStruct)
%ADDMLDCONSTRAINTS Adds state and output constraints to MLD matrices

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

% MLD model is of the form:
%
%  x(k+1) = A*x(k) + B_1*u(k) + B_2*d(t) + B_3*z(t)
%    y(k) = C*x(k) + D_1*u(k) + D_2*d(t) + D_3*z(t)
%  E_2*d(t) + E_3*z(t) <= E_1*u(k) + E_4*x(k) + E_5

if isfield(sysStruct, 'xmax') & S.nx > 0,
    % add state constraints
    xmax = sysStruct.xmax;
    xmin = sysStruct.xmin;
    if (all(isinf(xmax)) & all(isinf(xmin))),
        % do not add +/- Inf constraints
    else    
        nx = length(xmin);
        S.E4 = [S.E4; -eye(nx); eye(nx)];
        S.E5 = [S.E5; xmax; -xmin];
        S.E1 = [S.E1; zeros(2*nx, size(S.E1,2))];
        S.E2 = [S.E2; zeros(2*nx, size(S.E2,2))];
        S.E3 = [S.E3; zeros(2*nx, size(S.E3,2))];
        % we have added 2*nx inequality constraints, update global counter
        S.ne = S.ne + 2*nx;
    end
end

if isfield(sysStruct, 'ymax') & S.ny > 0,
    ymax = sysStruct.ymax;
    ymin = sysStruct.ymin;
    if (all(isinf(ymax)) & all(isinf(ymin))),
        % do not add +/- Inf constraints
    else
        ny = length(ymax);
        S.E1 = [S.E1; -S.D1; S.D1];
        S.E2 = [S.E2; S.D2; -S.D2];
        S.E3 = [S.E3; S.D3; -S.D3];
        S.E4 = [S.E4; -S.C; S.C];
        S.E5 = [S.E5; ymax; -ymin];
        % we have added 2*ny inequality constraints, update global counter
        S.ne = S.ne + 2*ny;
    end
end

