function [H, K, F, G, nc, nr, nx, nu, ny, nxt, nref, Ts, dumode, tracking, abstol] = mpt_getSfuncParam(ctrl)
%MPT_GETSFUNCPARAM Prepares parameters for explicit controller S-function
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Creates parameters to be passed to mpt_getInput_sfunc_parm, a
% S-function which simulates a given explicit controller.
%
% Internal function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
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
    
error(nargchk(1,1,nargin)); 

global mptOptions
if ~isstruct(mptOptions)
    mpt_error;
end

if isa(ctrl, 'polytope'),
    % input is a polytope
    P = ctrl;
    nx = dimension(P);
    nr = length(P);
    abstol = mptOptions.abs_tol;
    [Hn,Kn] = double(P);
    if ~iscell(Hn),
        Hn = {Hn};
        Kn = {Kn};
    end
    H = []; K = []; nc = [];
    for ii = 1:length(Hn),
        h = Hn{ii}';
        k = Kn{ii};
        H = [H; h(:)];
        K = [K; k(:)];
        nc = [nc; size(Hn{ii}, 1)];
    end
    HH = H;
    KK = K;
    NCC = nc;
    
    % map computed parameters to output arguments
    H = nr;
    K = nx;
    F = abstol;
    G = HH;
    nc = KK;
    nr = NCC;
    return
end
    
if ~isa(ctrl, 'mptctrl')
    error('Input must be an MPTCTRL controller object!');
end
if ~cancompile(ctrl),
    error('This controller cannot be exported to C code.');
end

nr = length(ctrl);
Pn = ctrl.Pn;
[Hn,Kn] = double(Pn);
if ~iscell(Hn),
    Hn = {Hn};
    Kn = {Kn};
end
nctotal = 0;
for ii=1:length(Pn),
    nctotal = nctotal + nconstr(Pn(ii));
end
nr = length(Pn);
Fi = ctrl.Fi;
Gi = ctrl.Gi;

nref = ctrl.details.x0format.reference;
nxt = ctrl.details.x0format.required - ...
    ctrl.details.x0format.reference - ...
    ctrl.details.x0format.uprev;
nu = ctrl.details.dims.nu;
ny = ctrl.details.dims.ny;
dumode = isfield(ctrl.sysStruct, 'dumode') | ...
    ctrl.probStruct.tracking==1 | ...
    ctrl.details.x0format.uprev>0;
tracking = ctrl.probStruct.tracking;
nx = ctrl.details.x0format.required;

% all boolean variables have to be represented as doubles, otherwise C-code
% S-functions do not recognize them
dumode = double(dumode);
tracking = double(tracking);

for ii=1:length(Fi),
    Fi{ii} = Fi{ii}(1:nu,:);
    Gi{ii} = Gi{ii}(1:nu);
end

if isfield(ctrl.sysStruct, 'Ts'),
    Ts = ctrl.sysStruct.Ts;
else
    Ts = 1;
end

abstol = mptOptions.abs_tol;

H = []; K = []; F = []; G = []; nc = [];
for ii = 1:length(Hn),
    h = full(Hn{ii}');
    k = full(Kn{ii});
    f = full(Fi{ii}');
    g = full(Gi{ii});
    H = [H; h(:)];
    K = [K; k(:)];
    F = [F; f(:)];
    G = [G; g(:)];
    nc = [nc; size(Hn{ii}, 1)];
end
