function [sysStruct, probStruct] = mpt_prepareTracking(sysStruct, probStruct)
%MPT_PREPARETRACKING Extends system and problem matrices to deal with tracking
%
% [sysStructTr, probStructTr] = mpt_prepareTracking(sysStruct, probStruct)
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
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

if ~isfield(sysStruct,'verified'),
    verOpt.verbose=1;
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end

if ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    probStruct=mpt_verifyProbStruct(probStruct,verOpt);
end


if isfield(probStruct,'yref'),
    % augmentation moved to mpt_constructMatrices
    return
end

if isfield(probStruct,'xref') | isfield(probStruct,'uref')
    % fixed state tracking
    % substitution is made, new state vector xn = x - xref
    % new input un = u - uref
    %
    % then it follows:
    % xn(k+1) = A xn(k) + B un(k) + (A xref + B uref + f - xref)

    [nx,nu] = mpt_sysStructInfo(sysStruct);
    
    if ~isfield(probStruct,'xref'),
        xref = zeros(nx,1);
    else
        xref = probStruct.xref(:);
    end
    if ~isfield(probStruct,'uref'),
        uref = zeros(nu,1);
    else
        uref = probStruct.uref(:);
    end
    if size(xref,1)~=nx,
        error(['xref has to be a column vector of length ' nx '!']);
    end
    if size(uref,1)~=nu,
        error(['uref has to be a column vector of length ' nx '!']);
    end
    
    if isfulldim(probStruct.Tset),
        % shift the terminal set. handle properly cases when probStruct.Tset is
        % a polytope array
        Tset = polytope;
        for ii = 1:length(probStruct.Tset),
            [Hset,Kset] = double(probStruct.Tset(ii));
            Tset = [Tset polytope(Hset, Kset - Hset*xref)];
        end
        probStruct.Tset = Tset;
    end
    
    if isfulldim(sysStruct.Pbnd),
        % shift Pbnd
        [Hbnd,Kbnd] = double(sysStruct.Pbnd);
        sysStruct.Pbnd = polytope(Hbnd, Kbnd - Hbnd*xref);
    end
    
    if iscell(sysStruct.A),
        for dyn=1:length(sysStruct.A),
            sysStruct.f{dyn} = sysStruct.A{dyn}*xref + sysStruct.B{dyn}*uref - xref + sysStruct.f{dyn};
            sysStruct.guardC{dyn} = sysStruct.guardC{dyn} - sysStruct.guardX{dyn}*xref - sysStruct.guardU{dyn}*uref;
        end
    end
    sysStruct.umax = sysStruct.umax - uref;
    sysStruct.umin = sysStruct.umin - uref;
    
    % output constraints are substituted in mpt_constructMatrices because they
    % depend on the output equation which could be different in different
    % segements for PWA systems 

    if any(~isinf(sysStruct.dumax)) | any(~isinf(sysStruct.dumin)),
        % include uref into xref if we have deltaU formulation
        probStruct.xref = [xref; uref];
        probStruct.uref = zeros(size(uref));
    else
        probStruct.xref = xref;
        probStruct.uref = uref;
    end
    probStruct.tracking = 0;
    probStruct.xref_augmented = 1;
    sysStructTr = sysStruct;
    probStructTr = probStruct;
    return
end
    

%------------- free state tracking ------------------------------ 

if probStruct.tracking==2
    fprintf('===============================================================================\n')
    fprintf('WARNING: offset-free tracking cannot be guaranteed with probStruct.tracking=2 !\n');
    fprintf('===============================================================================\n\n')
end

if ~isfield(probStruct,'Rdu') & probStruct.tracking==1,
    probStruct.Rdu=probStruct.R;
end

if iscell(sysStruct.A),
    if ~isfield(probStruct,'P_N'),
        if isfield(probStruct, 'Qy'),
            probStruct.P_N = probStruct.Qy;
            disp('mpt_prepareTracking: no final point penalty specified, using output penalty Qy...');
        else
            probStruct.P_N = probStruct.Q;
            disp('mpt_prepareTracking: no final point penalty specified, using state penalty Q...');
        end
    end
end

if (~isfield(probStruct,'Qy') | isempty(probStruct.Qy) | all(all(probStruct.Qy==0)))
    ycost=0;
else
    ycost=1;
end

[A,B,C,D,Q,R,ymin,ymax,umin,umax,dumin,dumax,bndA,bndb]=mpt_evalSystem(sysStruct,probStruct);

[nx,nu,ny,ndyn,nbool] = mpt_sysStructInfo(sysStruct);

sysStruct.dims.nx = nx;
sysStruct.dims.nu = nu;
sysStruct.dims.ny = ny;

if(ny~=nx & ~ycost)
    error('mpt_prepareTracking: Penalty on outputs probStruct.Qy must be defined if number of outputs is different from number of states.')
    %since otherwise there are no constraints defined for the reference state
end
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
if ispwa,
    for dyn=1:length(sysStruct.A),
        if ycost
            %reference state has dimension of output
            if probStruct.tracking==2,
                An{dyn} = [sysStruct.A{dyn} zeros(nx,ny); zeros(ny,nx) eye(ny)];
                Bn{dyn} = [sysStruct.B{dyn}; zeros(ny,nu)];
                Cn{dyn} = [sysStruct.C{dyn} zeros(ny,ny); zeros(ny,nx) eye(ny)];
                Dn{dyn} = [sysStruct.D{dyn}; zeros(ny,nu)];

                if isfield(sysStruct, 'Cy') & 0,
                    Cny{dyn} = [sysStruct.Cy{dyn} -eye(size(sysStruct.Cy{dyn},1))];
                    Dny{dyn} = sysStruct.Dy{dyn};
                else
                    Cny{dyn} = [sysStruct.C{dyn} -eye(ny)];
                    Dny{dyn} = sysStruct.D{dyn};
                end

            else
                An{dyn} = [sysStruct.A{dyn} sysStruct.B{dyn} zeros(nx,ny); zeros(nu,nx) eye(nu) zeros(nu,ny); zeros(ny,nx) zeros(ny,nu) eye(ny)];
                Bn{dyn} = [sysStruct.B{dyn}; eye(nu); zeros(ny,nu)];
                Cn{dyn} = [sysStruct.C{dyn} sysStruct.D{dyn} zeros(ny,ny); zeros(nu,nx) eye(nu) zeros(nu,ny); zeros(ny,nx+nu) eye(ny)];
                Dn{dyn} = [sysStruct.D{dyn}; zeros(nu); zeros(ny,nu)];
                
                Cny{dyn} = [sysStruct.C{dyn} zeros(ny,nu) -eye(ny); zeros(nu,nx) eye(nu) zeros(nu,ny)];
                Dny{dyn} = [sysStruct.D{dyn}; zeros(nu,nu)];
            end
        else
            %reference state has dimension of "normal" state
            if probStruct.tracking==2,
                An{dyn} = [sysStruct.A{dyn} zeros(nx,nx); zeros(nx) eye(nx)];
                Bn{dyn} = [sysStruct.B{dyn}; zeros(nx,nu)];
                Cn{dyn} = [sysStruct.C{dyn} zeros(ny,nx); zeros(nx) eye(nx)];
                Dn{dyn} = [sysStruct.D{dyn}; zeros(ny,nu)];
            else
                An{dyn} = [sysStruct.A{dyn} sysStruct.B{dyn} zeros(nx,nx); zeros(nu,nx) eye(nu) zeros(nu,nx); zeros(nx,nx+nu) eye(nx)];
                Bn{dyn} = [sysStruct.B{dyn}; eye(nu); zeros(nx,nu)];
                Cn{dyn} = [sysStruct.C{dyn} sysStruct.D{dyn} zeros(ny,nx); zeros(nu,nx) eye(nu) zeros(nu,nx); zeros(nx,nx+nu) eye(nx)];
                Dn{dyn} = [sysStruct.D{dyn}; zeros(nu); zeros(ny,nu)];
            end
        end % ycost
        
        
        if isfield(sysStruct,'f'),
            if probStruct.tracking==2,
                if ycost,
                    sysStruct.f{dyn} = [sysStruct.f{dyn}; zeros(ny,1)];
                else
                    sysStruct.f{dyn} = [sysStruct.f{dyn}; zeros(nx,1)];
                end
                sysStruct.g{dyn} = [sysStruct.g{dyn}; zeros(ny,1)];
            else
                if ycost,
                    sysStruct.f{dyn} = [sysStruct.f{dyn}; zeros(nu+ny,1)];
                else
                    sysStruct.f{dyn} = [sysStruct.f{dyn}; zeros(nu+nx,1)];
                end
                sysStruct.g{dyn} = [sysStruct.g{dyn}; zeros(nu+ny,1)];
            end
        else
            if probStruct.tracking==2,
                if ycost,
                    sysStruct.f{dyn} = zeros(nx+ny,1);
                else
                    sysStruct.f{dyn} = zeros(2*nx,1);
                end
            else
                if ycost,
                    sysStruct.f{dyn} = zeros(nx+ny+nu,1);
                else
                    sysStruct.f{dyn} = zeros(2*nx+nu,1);
                end
            end
        end
        f{dyn} = sysStruct.f{dyn};
        g{dyn} = sysStruct.g{dyn};
        
        if probStruct.tracking==2,
            [gXc, gXx] = size(sysStruct.guardX{dyn});
            guardX{dyn} = [sysStruct.guardX{dyn} zeros(gXc, length(f{dyn})-gXx)];
            guardU{dyn} = sysStruct.guardU{dyn};
        else
            if ycost,
                % the new state is [x u yref]
                guardX{dyn} = [sysStruct.guardX{dyn} sysStruct.guardU{dyn} zeros(size(sysStruct.guardX{dyn},1),ny)];
            else
                % the new state is [x u xref]
                guardX{dyn} = [sysStruct.guardX{dyn} sysStruct.guardU{dyn} zeros(size(sysStruct.guardX{dyn},1),nx)];
            end
            %sysStruct.guardU{dyn} = sysStruct.guardU{dyn}*0; %set to zero because the new input is delta U
            % sysStruct.guardU has to be included because otherwise the new control action will be disregarded. 
            guardU{dyn} = sysStruct.guardU{dyn};
        end
    end
else
    if ycost
        %reference state has dimension of output
        
        if probStruct.tracking==2,
            An = [sysStruct.A zeros(nx,ny); zeros(ny,nx) eye(ny)];
            Bn = [sysStruct.B; zeros(ny,nu)];
            Cn = [sysStruct.C zeros(ny,ny); zeros(ny,nx) eye(ny)];
            Dn = [sysStruct.D; zeros(ny,nu)];
            
            if isfield(sysStruct, 'Cy') & 0,
                Cny = [sysStruct.Cy -eye(size(sysStruct.Cy,1))];
                Dny = sysStruct.Dy;
            else
                Cny = [sysStruct.C -eye(ny)];
                Dny = sysStruct.D;
            end

        else
            An = [sysStruct.A sysStruct.B zeros(nx,ny); zeros(nu,nx) eye(nu) zeros(nu,ny); zeros(ny,nx) zeros(ny,nu) eye(ny)];
            Bn = [sysStruct.B; eye(nu); zeros(ny,nu)];
            Cn = [sysStruct.C sysStruct.D zeros(ny,ny); zeros(nu,nx) eye(nu) zeros(nu,ny); zeros(ny,nx+nu) eye(ny)];
            Dn = [sysStruct.D; zeros(nu); zeros(ny,nu)];
            
            Cny = [sysStruct.C zeros(ny,nu) -eye(ny); zeros(nu,nx) eye(nu) zeros(nu,ny)];
            Dny = [sysStruct.D; zeros(nu,nu)];
        end
    else
        %reference state has dimension of "normal" state
        if probStruct.tracking==2,
            An = [sysStruct.A zeros(nx); zeros(nx) eye(nx)];
            Bn = [sysStruct.B; zeros(nx,nu)];
            Cn = [sysStruct.C zeros(ny,nx); zeros(nx) eye(nx)];
            Dn = [sysStruct.D; zeros(nx,nu)];
        else
            An = [sysStruct.A sysStruct.B zeros(nx,nx); zeros(nu,nx) eye(nu) zeros(nu,nx); zeros(nx,nx+nu) eye(nx)];
            Bn = [sysStruct.B; eye(nu); zeros(nx,nu)];
            Cn = [sysStruct.C sysStruct.D zeros(ny,nx); zeros(nu,nx) eye(nu) zeros(nu,nx); zeros(nx,nx+nu) eye(nx)];
            Dn = [sysStruct.D; zeros(nu); zeros(nx,nu)];
        end
    end
    
    
    if isfield(sysStruct,'f'),
        if probStruct.tracking==2,
            if ycost,
                sysStruct.f = [sysStruct.f; zeros(ny,1)];
            else
                sysStruct.f = [sysStruct.f; zeros(nx,1)];
            end
        else
            if ycost,
                sysStruct.f = [sysStruct.f; zeros(nu+ny,1)];
            else
                sysStruct.f = [sysStruct.f; zeros(nu+nx,1)];
            end
        end
    else
        if probStruct.tracking==2,
            if ycost,
                sysStruct.f = zeros(nx+ny,1);
            else
                sysStruct.f = zeros(2*nx,1);
            end
        else
            if ycost,
                sysStruct.f = zeros(nx+ny+nu,1);
            else
                sysStruct.f = zeros(2*nx+nu,1);
            end
        end
    end
    f = sysStruct.f;
end

if isfield(probStruct, 'P_N')
    P = probStruct.P_N;
    if probStruct.norm==1,
        P = P/2;
    end
    if probStruct.tracking==2,
        if probStruct.norm==2,
            if ycost,
                P = [P P; P P];
            else
                P = [P -P; -P P];
            end
        else
            if ycost,
                P = [P P];
            else
                P = [P -P];
            end
        end
    else
        if ycost,
            if size(P,1)~=ny,
                error('In tracking, probStruct.P_N must be a matrix of the dimension of the output y!');
            end
            P = [P zeros(ny,nu) -P; zeros(nu, ny) zeros(nu) zeros(nu, ny); -P zeros(ny,nu) P];
        else
            if size(P,1)~=nx,
                error('In tracking, probStruct.P_N must be a matrix of the dimension of the state x!');
            end
            P = [P zeros(nx,nu) -P;zeros(nu,nx) zeros(nu) zeros(nu,nx); -P zeros(nx,nu) P]; %terminal cost
        end
    end
    probStruct.P_N = P;
end

% use given bounds on references, or use state/output constraints
if isfield(sysStruct, 'yrefmax'),
    yrefmax = sysStruct.yrefmax;
else
    yrefmax = sysStruct.ymax;
    yrefmax(find(yrefmax==Inf)) = mptOptions.infbox;
end

% use given bounds on references, or use state/output constraints
if isfield(sysStruct, 'yrefmin'),
    yrefmin = sysStruct.yrefmin;
else
    yrefmin = sysStruct.ymin;
    yrefmin(find(yrefmin==-Inf)) = -mptOptions.infbox;
end

if isfield(sysStruct, 'xmax'),
    % use given bounds on references, or use state/output constraints
    if isfield(sysStruct, 'xrefmax'),
        xrefmax = sysStruct.xrefmax;
    else
        xrefmax = sysStruct.xmax;
        xrefmax(find(xrefmax==Inf)) = mptOptions.infbox;
    end

    % use given bounds on references, or use state/output constraints
    if isfield(sysStruct, 'xrefmin'),
        xrefmin = sysStruct.xrefmin;
    else
        xrefmin = sysStruct.xmin;
        xrefmin(find(xrefmin==-Inf)) = -mptOptions.infbox;
    end
end

if probStruct.tracking==2,
    ymaxn = [sysStruct.ymax; yrefmax];
    yminn = [sysStruct.ymin; yrefmin];
else
    ymaxn = [sysStruct.ymax; sysStruct.umax; yrefmax];
    yminn = [sysStruct.ymin; sysStruct.umin; yrefmin];
end

if isfield(sysStruct, 'xmax')
    % augment state constraints if present
    xmax = sysStruct.xmax;
    xmin = sysStruct.xmin;
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
            
%the cost is updated accordingly to punish (x-xref)' Q (x-xref) + delta_u' Rdu delta_u
%       [Q   0  -Q]
% Qn =  [0   0   0]
%       [-Q  0   Q]

if probStruct.tracking==2
    if ycost
        if(length(probStruct.Qy)==ny) | 1
            QQy = probStruct.Qy;
            if probStruct.norm==2,
                %probStruct.Qy = [QQy zeros(ny); zeros(ny) zeros(ny)];
                probStruct.Qy = [QQy zeros(size(QQy,1)); zeros(size(QQy,1)) QQy];
                Qn = probStruct.Qy;
            else
                probStruct.Qy = [QQy QQy];
            end
            probStruct.Qy = QQy;
        end
        Qn = probStruct.Q;
    else
        if probStruct.norm==2,
            QQ = probStruct.Q;
            Qn = [QQ -QQ; -QQ QQ];
        else
            QQ = probStruct.Q;
            Qn = [QQ -QQ];
        end
    end
else
    if ycost
        if(length(probStruct.Qy)==ny)
            probStruct.Qy = [probStruct.Qy zeros(ny,nu) zeros(ny); zeros(nu,ny) zeros(nu) zeros(nu,ny); zeros(ny) zeros(ny,nu) probStruct.Qy];
        end
        Qn = probStruct.Q;
    else
        %%Qn = [probStruct.Q zeros(nx,nu) -probStruct.Q; zeros(nu,nx) probStruct.R zeros(nu,nx); -probStruct.Q zeros(nx,nu) probStruct.Q];
        if probStruct.norm==1,
            QQ = probStruct.Q/2;
        else
            QQ = probStruct.Q;
        end
        Qn = [QQ zeros(nx,nu) -QQ; zeros(nu,nx) zeros(nu) zeros(nu,nx); -QQ zeros(nx,nu) QQ];
    end
end

%Update the outer bound Pbnd
if probStruct.tracking==2,
    [Hbnd, Kbnd] = double(sysStruct.Pbnd);
    if ycost,
        Haux = [eye(ny); -eye(ny)];
        Kaux = [sysStruct.ymax; -sysStruct.ymin];
        sysStruct.Pbnd = polytope( [Hbnd zeros(size(Hbnd,1),ny); zeros(2*ny, size(Hbnd,2)) Haux], [Kbnd; Kaux]);
    else
        sysStruct.Pbnd = polytope( [Hbnd zeros(size(Hbnd)); zeros(size(Hbnd)) Hbnd], [Kbnd; Kbnd]);
    end
else
    [Hbnd,Kbnd] = double(sysStruct.Pbnd);
    [nxc,nxx] = size(Hbnd);
    nuc = 2*size(umax,1);
    I2u = [eye(nu); -eye(nu)];
    I2y = [eye(ny); -eye(ny)];
    Hn = [Hbnd zeros(nxc,nu) zeros(nxc,ny); zeros(nuc,nx) I2u zeros(nuc,ny); zeros(2*ny,nx) zeros(2*ny,nu) I2y];
    Kn = [Kbnd; umax; -umin; ymax; -ymin];
    sysStruct.Pbnd = polytope(Hn,Kn);
end

if ispwa & isfield(probStruct,'Tset'),
    if isfulldim(probStruct.Tset),
        Tset = polytope;
        for ii = 1:length(probStruct.Tset),
            [Hs,Ks]=double(probStruct.Tset(ii));
            [nxc,nxx] = size(Hs);
            if probStruct.tracking==2,
                Hn = [Hs zeros(size(Hs)); zeros(size(Hs)) -Hs];
                Kn = [Ks; Ks];
            else
                Hn = [Hs zeros(nxc,nu) zeros(nxc,nxx);...
                        zeros(nuc,nxx) I2u zeros(nuc,nxx);
                    zeros(nxc,nxx) zeros(nxc,nu) -Hs];
                Kn = [Ks; umax; -umin; Ks];
            end
            Tset = [Tset polytope(Hn,Kn)];
        end
        probStruct.Tset = Tset;
    end
end

%write tracking data back into structure
sysStruct.A = An;
sysStruct.B = Bn;
sysStruct.C = Cn;
sysStruct.D = Dn;
if ycost
    sysStruct.Cy = Cny;
    sysStruct.Dy = Dny;
end
sysStruct.ymax = ymaxn;
sysStruct.ymin = yminn;

if probStruct.tracking==1,
    sysStruct.umax = sysStruct.dumax;
    sysStruct.umin = sysStruct.dumin;
    sysStruct.dumax = Inf*ones(nu,1);
    sysStruct.dumin = -Inf*ones(nu,1);
    probStruct.R = probStruct.Rdu;
end

probStruct.Q = Qn;
sysStruct = rmfield(sysStruct,'verified');

if ispwa,
    sysStruct.f = f;
    sysStruct.g = g;
    sysStruct.guardX = guardX;
    sysStruct.guardU = guardU;
end
    
probStruct.tracking_augmented = probStruct.tracking;
tmpopt.verbose = 0;
evalc('sysStruct = mpt_verifySysStruct(sysStruct,tmpopt);');
evalc('probStruct = mpt_verifyProbStruct(probStruct,tmpopt);');
sysStructTr = sysStruct;
probStructTr = probStruct;
