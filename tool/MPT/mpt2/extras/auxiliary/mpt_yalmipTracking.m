function [sysStruct, probStruct] = mpt_yalmipTracking(sysStruct, probStruct, verOpt)
%MPT_YALMIPTRACKING Extends system and problem matrices to deal with tracking
%
% [sysStructTr, probStructTr] = mpt_yalmipTracking(sysStruct, probStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Extends system and problem matrices to deal with tracking
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
% sysStructTr   - system definition with tracking
% probStructTr  - problem definition with tracking
%        .Rdu   - additional field in probStruct; weight on the delta u; 
%                 by default, this is identical to the weight on u;
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

if isfield(probStruct,'yref') | isfield(probStruct,'xref') | ...
        isfield(probStruct,'uref') | probStruct.tracking==0,
    % nothing to do here, we handle references directly in mpt_yalmipcftoc()
    return
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

%------------- free reference tracking ------------------------------ 

[nx,nu,ny,ndyn,nbool] = mpt_sysStructInfo(sysStruct);
sysStruct.dims.nx = nx;
sysStruct.dims.nu = nu;
sysStruct.dims.ny = ny;
use_deltaU = isfield(probStruct, 'Rdu') | ...
    ~(all(isinf(sysStruct.dumax)) & all(isinf(sysStruct.dumin)));

if nbool > 0 & (probStruct.tracking==1 | use_deltaU),
    fprintf('\ndeltaU constraints/penalties cannot be used for systems with boolean inputs.\n');
    error('Please use "probStruct.tracking=2".');

elseif probStruct.tracking==2 & verOpt.verbose > -1
    fprintf('================================================================================\n')
    fprintf('WARNING: offset-free tracking cannot be guaranteed with probStruct.tracking=2 !\n');
    if use_deltaU,
        fprintf('WARNING: deltaU constraint/penalties only guaranteed in on-line control !\n');
    end
    fprintf('================================================================================\n\n')
    
end

if ~isfield(probStruct,'Rdu') & probStruct.tracking==1,
    probStruct.Rdu = probStruct.R;
end

if (~isfield(probStruct,'Qy') | isempty(probStruct.Qy) | all(all(probStruct.Qy==0)))
    ycost=0;
else
    ycost=1;
end

if(ny~=nx & ~ycost)
    error('mpt_prepareTracking: Penalty on outputs probStruct.Qy must be defined if number of outputs is different from number of states.')
    % otherwise there are no constraints defined for the reference state
end


if ycost
    %reference state has dimension of output
    refdim = ny;
else
    refdim = nx;
end
sysStruct.dims.nref = refdim;


%++++++++++++++++++++++++++++++++++++
% augment matrices for tracking
%++++++++++++++++++++++++++++++++++++
%
% Introduce new state vector z(k) = [x(k) u(k-1) xref(k)], where the reference state xref is user defined
% The new input is now delta u, i.e., du(k)=u(k)-u(k-1).
% Therefore the state update equation can now be written as:
%        [A B 0]         [B] 
%z(k+1)= [0 I 0] z(k) +  [I] du(k)
%        [0 0 I]         [0]

ispwa = iscell(sysStruct.A);
if ~ispwa,
    sysStruct = mpt_lti2pwa(sysStruct);
end

for dyn=1:length(sysStruct.A),
    
    if probStruct.tracking==2,
        An{dyn} = [sysStruct.A{dyn} zeros(nx,ny); ...
                zeros(refdim,nx) eye(refdim)];
        
        Bn{dyn} = [sysStruct.B{dyn}; zeros(refdim,nu)];
        
        Cn{dyn} = [sysStruct.C{dyn} zeros(ny,refdim); ...
                zeros(refdim,nx) eye(refdim)];
        
        Dn{dyn} = [sysStruct.D{dyn}; zeros(refdim,nu)];
        
        if isfield(sysStruct, 'f'),
            fn{dyn} = [sysStruct.f{dyn}; zeros(refdim, 1)];
        end
        if isfield(sysStruct, 'g'),
            gn{dyn} = [sysStruct.g{dyn}; zeros(refdim, 1)];
        end
        
    else
        An{dyn} = [sysStruct.A{dyn} sysStruct.B{dyn} zeros(nx,refdim); ...
                zeros(nu,nx) eye(nu) zeros(nu,refdim); ...
                zeros(refdim,nx) zeros(refdim,nu) eye(refdim)];
        
        Bn{dyn} = [sysStruct.B{dyn}; eye(nu); zeros(refdim,nu)];
        
        Cn{dyn} = [sysStruct.C{dyn} sysStruct.D{dyn} zeros(ny,refdim); ...
                zeros(nu,nx) eye(nu) zeros(nu,refdim); ...
                zeros(refdim,nx+nu) eye(refdim)];
        
        Dn{dyn} = [sysStruct.D{dyn}; zeros(nu); zeros(refdim,nu)];
        
        if isfield(sysStruct, 'f'),
            fn{dyn} = [sysStruct.f{dyn}; zeros(nu+refdim, 1)];
        end
        if isfield(sysStruct, 'g'),
            gn{dyn} = [sysStruct.g{dyn}; zeros(nu+refdim, 1)];
        end
        
    end
    
    if probStruct.tracking==2,
        guardX{dyn} = [sysStruct.guardX{dyn} zeros(size(sysStruct.guardX{dyn},1), refdim)];
        guardU{dyn} = sysStruct.guardU{dyn};
        
    else
        guardX{dyn} = [sysStruct.guardX{dyn} sysStruct.guardU{dyn} ...
                zeros(size(sysStruct.guardX{dyn},1),refdim)];
        
        %sysStruct.guardU{dyn} = sysStruct.guardU{dyn}*0; %set to zero because the new input is delta U
        % sysStruct.guardU has to be included because otherwise the new control action will be disregarded. 
        guardU{dyn} = sysStruct.guardU{dyn};
    end
end


%===================================================================
% update constraints

haveXbounds = isfield(sysStruct, 'xmax');
haveYbounds = isfield(sysStruct, 'ymax');

% since reference is now part of the state vector, it is always a good idea to
% keep all elements bounded. therefore if we are tracking the system states
% (outputs), we use state (output) constraints as bounds on references
if haveXbounds,
    default_xrefmax = sysStruct.xmax;
    default_xrefmin = sysStruct.xmin;
else
    default_xrefmax = repmat(mptOptions.infbox, nx, 1);
    default_xrefmin = -repmat(mptOptions.infbox, nx, 1);
end
if haveYbounds,
    default_yrefmax = sysStruct.ymax;
    default_yrefmin = sysStruct.ymin;
else
    default_yrefmax = repmat(mptOptions.infbox, ny, 1);
    default_yrefmin = -repmat(mptOptions.infbox, ny, 1);
end
xrefmax = mpt_defaultField(sysStruct, 'xrefmax', default_xrefmax);
xrefmin = mpt_defaultField(sysStruct, 'xrefmin', default_xrefmin);
yrefmax = mpt_defaultField(sysStruct, 'yrefmax', default_yrefmax);
yrefmin = mpt_defaultField(sysStruct, 'yrefmin', default_yrefmin);

% replace any +/- Inf bounds with the value of mptOptions.infbox to keep at
% least some level of scaling
xrefmax(find(xrefmax==Inf)) = mptOptions.infbox;
xrefmin(find(xrefmin==-Inf)) = -mptOptions.infbox;
yrefmax(find(yrefmax==Inf)) = mptOptions.infbox;
yrefmin(find(yrefmin==-Inf)) = -mptOptions.infbox;

% use given bounds on references, or use state/output constraints
if haveXbounds,
    if probStruct.tracking==2
        if ycost
            % state vector has dimension nx+ny
            xmaxn = [sysStruct.xmax; yrefmax];
            xminn = [sysStruct.xmin; yrefmin];
        else
            % state vector has dimension nx+nx
            xmaxn = [sysStruct.xmax; xrefmax];
            xminn = [sysStruct.xmin; xrefmin];
        end
    else
        if ycost
            % state vector has dimension nx+nu+ny
            xmaxn = [sysStruct.xmax; sysStruct.umax; yrefmax];
            xminn = [sysStruct.xmin; sysStruct.umin; yrefmin];
        else
            % state vector has dimension 2*nx+nu
            xmaxn = [sysStruct.xmax; sysStruct.umax; xrefmax];
            xminn = [sysStruct.xmin; sysStruct.umin; xrefmin];
        end
    end
    sysStruct.xmax = xmaxn;
    sysStruct.xmin = xminn;
end

% use given bounds on references, or use state/output constraints
if haveYbounds,
    if probStruct.tracking==2,
        ymaxn = [sysStruct.ymax; yrefmax];
        yminn = [sysStruct.ymin; yrefmin];
    else
        ymaxn = [sysStruct.ymax; sysStruct.umax; yrefmax];
        yminn = [sysStruct.ymin; sysStruct.umin; yrefmin];
    end
end

sxmax = mpt_defaultField(probStruct, 'sxmax', zeros(nx, 1));
sumax = mpt_defaultField(probStruct, 'sumax', zeros(nu, 1));
symax = mpt_defaultField(probStruct, 'symax', zeros(ny, 1));
Sx = mpt_defaultField(probStruct, 'Sx', 1000*eye(nx));
Su = mpt_defaultField(probStruct, 'Su', 1000*eye(nu));
Sy = mpt_defaultField(probStruct, 'Sy', 1000*eye(ny));
% now augment bounds on slack variables
if isfield(probStruct, 'sxmax') | isfield(probStruct, 'Sx'),
    if probStruct.tracking==2,
        if ycost,
            n_sxmax = [sxmax; zeros(ny, 1)];
            n_Sx = [Sx zeros(nx, ny); zeros(ny, nx) Sy];
        else
            n_sxmax = [sxmax; zeros(nx, 1)];
            n_Sx = [Sx zeros(nx, nx); zeros(nx, nx) Sx];
        end
        
    else
        if ycost,
            n_sxmax = [sxmax; sumax; zeros(ny, 1)];
            n_Sx = [Sx zeros(nx, nu) zeros(nx, ny); ...
                    zeros(nu, nx) Su zeros(nu, ny); ...
                    zeros(ny, nx) zeros(ny, nu) Sy];

        else
            n_sxmax = [sxmax; sumax; zeros(nx, 1)];
            n_Sx = [Sx zeros(nx, nu) zeros(nx, nx); ...
                    zeros(nu, nx) Su zeros(nu, nx); ...
                    zeros(nx, nx) zeros(nx, nu) Sx];
            
        end
    end
    probStruct.sxmax = n_sxmax;
    probStruct.Sx = n_Sx;
end
if isfield(probStruct, 'symax') | isfield(probStruct, 'Sy'),
    if probStruct.tracking==2,
        if ycost,
            n_symax = [symax; zeros(ny, 1)];
            n_Sy = [Sy zeros(ny, ny); zeros(ny, ny) Sy];
            
        else
            n_symax = [symax; zeros(nx, 1)];
            n_Sy = [Sy zeros(ny, nx); zeros(nx, ny) Sx];
            
        end
    else
        if ycost,
            n_symax = [symax; sumax; zeros(ny, 1)];
            n_Sy = [Sy zeros(ny, nu) zeros(ny, ny); ...
                    zeros(nu, ny) Su zeros(nu, ny); ...
                    zeros(ny, ny) zeros(ny, nu) Sy];
            
        else
            n_symax = [symax; sumax; zeros(nx, 1)];
            n_Sy = [Sy zeros(ny, nu) zeros(ny, nx); ...
                    zeros(nu, ny) Su zeros(nu, nx); ...
                    zeros(nx, ny) zeros(nx, nu) Sx];
            
        end
    end
    probStruct.symax = n_symax;
    probStruct.Sy = n_Sy;
end


%===================================================================
% update Pbnd
if ycost,
    Br = polytope([eye(refdim); -eye(refdim)], [yrefmax; -yrefmin]);
else
    Br = polytope([eye(refdim); -eye(refdim)], [xrefmax; -xrefmin]);
end
Bu = polytope([eye(nu); -eye(nu)], [sysStruct.umax; -sysStruct.umin]);
% create polytopes in higher dimension
if probStruct.tracking==2,
    sysStruct.Pbnd = sysStruct.Pbnd * Br;
else
    sysStruct.Pbnd = sysStruct.Pbnd * Bu * Br;
end


%===================================================================
% update penalties

if isfield(probStruct, 'P_N')
    % P = [P -P; -P P]   for 2 norm
    % P = [P -P]         otherwise
    probStruct.P_N = sub_augmentpenalty(probStruct.P_N, probStruct, nu);
end

%the cost is updated accordingly to punish (x-xref)' Q (x-xref) + delta_u' Rdu delta_u
%       [ Q  0  -Q]
% Qn =  [ 0  0   0]
%       [-Q  0   Q]
if ycost,
    % if probStruct.Qy is defined, we minimize ||Qy*(y-yref)|| , i.e. there will
    % be no penalty w.r.t. states. therefore we can set probStruct.Q to anything
    if probStruct.tracking == 1,
        probStruct.Q = eye(nx+nu+refdim);
    else
        probStruct.Q = eye(nx+refdim);
    end
else
    probStruct.Q = sub_augmentpenalty(probStruct.Q, probStruct, nu);
end

if isfield(probStruct, 'Qy') & ycost,
    probStruct.Qy = sub_augmentpenalty(probStruct.Qy, probStruct, nu);
end

if isfield(probStruct, 'Tset'),
    if isfulldim(probStruct.Tset),
        Tset = polytope;
        for ii = 1:length(probStruct.Tset)
            if probStruct.tracking==2,
                Tset = [Tset probStruct.Tset(ii) * Br];
            else
                Tset = [Tset probStruct.Tset(ii) * Bu * Br];
            end
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
sysStruct.ymax = ymaxn;
sysStruct.ymin = yminn;
if ispwa,
    sysStruct.f = fn;
    sysStruct.g = gn;
    sysStruct.guardX = guardX;
    sysStruct.guardU = guardU;
end

if probStruct.tracking==1,
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
    probStruct.R = probStruct.Rdu;
    probStruct = rmfield(probStruct, 'Rdu');
end


sysStruct = rmfield(sysStruct,'verified');
probStruct = rmfield(probStruct, 'verified');
verOpt.verbose = verOpt.verbose - 1;
evalc('sysStruct = mpt_verifySysStruct(sysStruct, verOpt);');
evalc('probStruct = mpt_verifyProbStruct(probStruct, verOpt);');
sysStructTr = sysStruct;
probStruct.tracking_augmented = 1;  % keep a note that the system has already been augmented
probStructTr = probStruct;



%----------------------------------------------------------
function P = sub_augmentpenalty(P, probStruct, nu)

if ~iscell(P),
    P = {P};
end

Porig = P;
for ii = 1:length(Porig),
    P = Porig{ii};
    [nr,nc] = size(P);
    if probStruct.tracking==2,
        if probStruct.norm==2,
            P = [P -P; -P P];
        else
            P = [P -P];
        end
        
    else
        % tracking=1
        if probStruct.norm==2,
            P = [P zeros(nr, nu) -P; ...
                    zeros(nu, nc) zeros(nu) zeros(nu, nc); ...
                    -P zeros(nr, nu) P];
            
        else
            P = [P zeros(nr, nu) -P];
            
        end
        
    end
    Porig{ii} = P;
end
