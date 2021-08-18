function [sysStruct, probStruct] = mpt_yalmipDU(sysStruct, probStruct, verOpt)
%MPT_YALMIPDU Augmentes the system to cope with deltaU constraints
%
% [sysStruct, probStruct] = mpt_yalmipDU(sysStruct, probStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Extends system and problem matrices to deal with deltaU constraints
%
% internal function
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% sysStruct   - system definition
% probStruct  - problem definition
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sysStruct   - augmented system structure
% probStruct  - augmented problem structure
%

% Copyright is with the following author(s):
%
% (C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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
if ~isstruct(mptOptions),
    mpt_error;
end
if nargin < 3,
    verOpt = [];
end
if ~isfield(verOpt, 'verbose'),
    verOpt.verbose = mptOptions.verbose;
end

if ~isfield(sysStruct,'verified'),
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end

if ~isfield(probStruct,'verified'),
    probStruct=mpt_verifyProbStruct(probStruct,verOpt);
end

if isfield(sysStruct, 'data'),
    if isfield(sysStruct.data, 'onlymld'),
        if sysStruct.data.onlymld,
            % cannot compute an explicit controller if PWA model is not
            % available
            fprintf('\nPWA representation of the hybrid system must be available in order to use tracking.\n');
            fprintf('Call "sysStruct = mpt_sys(sysStruct.data.MLD)" to get the PWA representation.\n\n');
            error('Cannot deal with tracking if PWA representation is not available in sysStruct.');
        end
    end
end

%------------- deltaU formulation ------------------------------ 

[nx,nu,ny,ndyn,nbool] = mpt_sysStructInfo(sysStruct);
sysStruct.dims.nx = nx;
sysStruct.dims.nu = nu;
sysStruct.dims.ny = ny;

if probStruct.tracking > 0,
    error('For tracking, call mpt_yalmipTracking() instead.');
end

if all(isinf(sysStruct.dumax)) & any(isinf(sysStruct.dumin)) & ~isfield(probStruct, 'Rdu'),
    % nothing to do, no deltaU constraints nor penalty on them
    fprintf('mpt_yalmipDU: nothing to do here.\n');
    return
end

if nbool > 0,
    error('deltaU constraints/penalties cannot be used for systems with boolean inputs.');
end

% augment possible references
if isfield(probStruct, 'xref'),
    probStruct.xref = [probStruct.xref; zeros(nu, 1)];
end
if isfield(probStruct, 'yref'),
    probStruct.yref = [probStruct.yref; zeros(nu, 1)];
end

%+++++++++++++++++++++++++++++++++++++++++++++++++
% augment matrices to deal with deltaU formulation
%+++++++++++++++++++++++++++++++++++++++++++++++++
%
% Introduce new state vector z(k) = [x(k) u(k-1)].
% The new input is now delta u, i.e., du(k)=u(k)-u(k-1).
% Therefore the state update equation can now be written as:
%        [A B]         [B] 
%z(k+1)= [0 I] z(k) +  [I] du(k)
%
%[ y(k) ]   [C D]        [D]
%[u(k-1)] = [0 I] z(k) + [0] du(k)
ispwa = iscell(sysStruct.A);
if ~ispwa,
    sysStruct = mpt_lti2pwa(sysStruct);
end

for dyn=1:length(sysStruct.A),
    
    An{dyn} = [sysStruct.A{dyn} sysStruct.B{dyn}; ...
            zeros(nu,nx) eye(nu)];
    
    Bn{dyn} = [sysStruct.B{dyn}; eye(nu)];
    
    Cn{dyn} = [sysStruct.C{dyn} sysStruct.D{dyn}; ...
            zeros(nu, nx) eye(nu)];
    
    Dn{dyn} = [sysStruct.D{dyn}; zeros(nu)];
    
    if isfield(sysStruct, 'f'),
        fn{dyn} = [sysStruct.f{dyn}; zeros(nu, 1)];
    end
    if isfield(sysStruct, 'g'),
        gn{dyn} = [sysStruct.g{dyn}; zeros(nu, 1)];
    end

    guardX{dyn} = [sysStruct.guardX{dyn} sysStruct.guardU{dyn}];
    guardU{dyn} = sysStruct.guardU{dyn};
    guardX{dyn} = [guardX{dyn}; zeros(nu, nx) eye(nu); zeros(nu, nx) -eye(nu)];
    guardU{dyn} = [guardU{dyn}; zeros(2*nu, nu)];
    guardC{dyn} = [sysStruct.guardC{dyn}; sysStruct.umax; -sysStruct.umin];
end

%===================================================================
% update constraints

haveXbounds = isfield(sysStruct, 'xmax');
haveYbounds = isfield(sysStruct, 'ymax');

% use given bounds on references, or use state/output constraints
if haveXbounds,
    sysStruct.xmax = [sysStruct.xmax; sysStruct.umax];
    sysStruct.xmin = [sysStruct.xmin; sysStruct.umin];
end
if haveYbounds,
    sysStruct.ymax = [sysStruct.ymax; sysStruct.umax];
    sysStruct.ymin = [sysStruct.ymin; sysStruct.umin];
end

%===================================================================
% update slacks

sxmax = mpt_defaultField(probStruct, 'sxmax', zeros(nx, 1));
sumax = mpt_defaultField(probStruct, 'sumax', zeros(nu, 1));
symax = mpt_defaultField(probStruct, 'symax', zeros(ny, 1));
Sx = mpt_defaultField(probStruct, 'Sx', 1000*eye(nx));
Su = mpt_defaultField(probStruct, 'Su', 1000*eye(nu));
Sy = mpt_defaultField(probStruct, 'Sy', 1000*eye(ny));

sxmax = [sxmax; sumax];
Sx = [Sx zeros(nx, nu); zeros(nu, nx) Su];
symax = [symax; sumax];
Sy = [Sy zeros(ny, nu); zeros(nu, ny) Su];
if isfield(probStruct, 'sxmax') | isfield(probStruct, 'Sx'),
    probStruct.sxmax = sxmax;
    probStruct.Sx = Sx;
end
if isfield(probStruct, 'symax') | isfield(probStruct, 'Sy'),
    probStruct.symax = symax;
    probStruct.Sy = Sy;
end

%===================================================================
% update Pbnd
Bu = polytope([eye(nu); -eye(nu)], [sysStruct.umax; -sysStruct.umin]);
sysStruct.Pbnd = sysStruct.Pbnd * Bu;


%===================================================================
% update penalties

if isfield(probStruct, 'P_N')
    % P = [P 0; 0 0]
    [pnr, pnc] = size(probStruct.P_N);
    probStruct.P_N = [probStruct.P_N zeros(pnr, nu); zeros(nu, pnc) zeros(nu)];
end

if isfield(probStruct, 'Qy')
    % Qy = [Q 0; 0 R]
    [qyr, qyc] = size(probStruct.Qy);
    probStruct.Qy = [probStruct.Qy zeros(qyr, nu); zeros(nu, qyc) probStruct.R];
end

%       [ Q  0 ]
% Qn =  [ 0  R ]
[qr, qc] = size(probStruct.Q);
probStruct.Q = [probStruct.Q zeros(qr, nu); zeros(nu, qc) probStruct.R];

if isfield(probStruct, 'Tset'),
    if isfulldim(probStruct.Tset),
        Tset = polytope;
        for ii = 1:length(probStruct.Tset)
            Tset = [Tset probStruct.Tset(ii) * Bu];
        end
        probStruct.Tset = Tset;
    end
end

%write tracking data back into structure
if ~ispwa,
    An = An{1}; Bn = Bn{1}; Cn = Cn{1}; Dn = Dn{1};
    sysStruct = rmfield(sysStruct, 'f');
    sysStruct = rmfield(sysStruct, 'g');
    sysStruct = rmfield(sysStruct, 'guardX');
    sysStruct = rmfield(sysStruct, 'guardU');
    sysStruct = rmfield(sysStruct, 'guardC');
end
sysStruct.A = An;
sysStruct.B = Bn;
sysStruct.C = Cn;
sysStruct.D = Dn;
if ispwa,
    sysStruct.f = fn;
    sysStruct.g = gn;
    sysStruct.guardX = guardX;
    sysStruct.guardU = guardU;
    sysStruct.guardC = guardC;
end

if all(isinf(sysStruct.dumax)) & all(isinf(sysStruct.dumin)),
    % in tracking=1 we go for deltaU formulation, therefore we take the
    % constraints on deltaU from sysStruct.{dumin|dumax}. however it is very
    % frequent that these bounds are +/- Inf (either specified like that by
    % the user or set automatically by mpt_verifySysStruct(). therefore if
    % we detect such case, we compute the bounds on deltaU from constraints
    % on "u" itself
    umin = sysStruct.umin - sysStruct.umax;    
    umax = sysStruct.umax - sysStruct.umin;
else
    umin = sysStruct.dumin;
    umax = sysStruct.dumax;
end
sysStruct.umax = umax;
sysStruct.umin = umin;
sysStruct.dumax = Inf*ones(nu,1);
sysStruct.dumin = -Inf*ones(nu,1);
if isfield(probStruct, 'Rdu'),
    probStruct.R = probStruct.Rdu;
    probStruct = rmfield(probStruct, 'Rdu');
end

probStruct = rmfield(probStruct, 'verified');
sysStruct = rmfield(sysStruct, 'verified');
verOpt.verbose = verOpt.verbose - 1;
evalc('sysStruct = mpt_verifySysStruct(sysStruct, verOpt);');
evalc('probStruct = mpt_verifyProbStruct(probStruct, verOpt);');
sysStruct.dumode = 1;  % keep note that we augmented the system dynamics
