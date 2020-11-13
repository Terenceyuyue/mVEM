function data = mpt_mexData(ctrl)
%MPT_MEXDATA Prepares input data for mex_getInput
%
% mpt_exportc(ctrl, fname)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Prepares input data for mex_getInput
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl   - MPT explicit controller
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% data   - Structure containing data of a given explicit controller
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

if ~isa(ctrl, 'mptctrl')
    error('Input argument must be a MPT controller!');
end
if ~cancompile(ctrl),
    error('This controller cannot be exported.');
end
    
[H, K, F, G, nc, nr, nx, nu, ny, nxt, nref, Ts, dumode, tracking, abstol] = mpt_getSfuncParam(ctrl);

data.H = H;
data.K = K;
data.F = F;
data.G = G;
data.nc = nc;
data.nr = nr;
data.nx = nx;
data.nu = nu;
