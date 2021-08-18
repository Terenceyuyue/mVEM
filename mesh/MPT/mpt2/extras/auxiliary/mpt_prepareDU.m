function [sysStruct, probStruct] = mpt_prepareDU(sysStruct, probStruct)
%MPT_PREPAREDU Extends system and problem matrices to deal with deltaU constraints
%
% [sysStruct, probStruct] = mpt_prepareDU(sysStruct, probStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Extends system and problem matrices to deal with deltaU constraints in
% closed-loop
%
% Introduce new state vector z(k) = [x(k) u(k-1)]
% The new input is now delta u, i.e., du(k)=u(k)-u(k-1).
% Therefore the state update equation can now be written as:
%
%         [A B]         [B] 
% z(k+1)= [0 I] z(k) +  [I] du(k)
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
% sysStruct   - system definition with tracking
% probStruct  - problem definition with tracking
%      .Rdu   - additional field in probStruct; weight on the delta u; 
%               by default, this is identical to the weight on u;
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

if ~isfield(sysStruct,'verified'),
    verOpt.verbose=1;
    sysStruct=mpt_verifySysStruct(sysStruct,verOpt);
end

if ~isfield(probStruct,'verified'),
    verOpt.verbose=1;
    probStruct=mpt_verifyProbStruct(probStruct,verOpt);
end

%++++++++++++++++++++++++++++++++++++++++
% augment matrices for deltaU constraints
%++++++++++++++++++++++++++++++++++++++++
%
% Introduce new state vector z(k) = [x(k) u(k-1)]
% The new input is now delta u, i.e., du(k)=u(k)-u(k-1).
% Therefore the state update equation can now be written as:
%        [A B]         [B] 
%z(k+1)= [0 I] z(k) +  [I] du(k)

[nx,nu,ny,ndyn] = mpt_sysStructInfo(sysStruct);

if (~isfield(probStruct,'Qy') | isempty(probStruct.Qy) | all(all(probStruct.Qy==0)))
    ycost=0;
else
    ycost=1;
end

ispwa = iscell(sysStruct.A);
if ispwa,
    for dyn=1:length(sysStruct.A),
        An{dyn} = [sysStruct.A{dyn} sysStruct.B{dyn}; zeros(nu,nx) eye(nu)];
        Bn{dyn} = [sysStruct.B{dyn}; eye(nu)];
        Cn{dyn} = [sysStruct.C{dyn} sysStruct.D{dyn}; zeros(nu,nx) eye(nu)];
        Dn{dyn} = [sysStruct.D{dyn}; zeros(nu)];
        
        if isfield(sysStruct,'f'),
            sysStruct.f{dyn} = [sysStruct.f{dyn}; zeros(nu,1)];
            sysStruct.g{dyn} = [sysStruct.g{dyn}; zeros(nu,1)];
        else
            sysStruct.f{dyn} = zeros(nx+nu,1);
        end
        f{dyn} = sysStruct.f{dyn};
        g{dyn} = sysStruct.g{dyn};
        guardX{dyn} = [sysStruct.guardX{dyn} sysStruct.guardU{dyn}]; %the new state is [x u xref], therefore update the guardlines
        %sysStruct.guardU{dyn} = sysStruct.guardU{dyn}*0; %set to zero because the new input is delta U
        % sysStruct.guardU has to be included because otherwise the new control action will be disregarded. 
        guardU{dyn} = sysStruct.guardU{dyn};
    end

%     if isfield(probStruct, 'P_N'),
%         nn = size(probStruct.P_N,1);
%         probStruct.P_N = [probStruct.P_N zeros(nx,nu);zeros(nu,nx) probStruct.R]; %terminal cost
%     end
else
    An = [sysStruct.A sysStruct.B; zeros(nu,nx) eye(nu)];
    Bn = [sysStruct.B; eye(nu)];
    Cn = [sysStruct.C sysStruct.D; zeros(nu,nx) eye(nu)];
    Dn = [sysStruct.D; zeros(nu)];
    Cy_augm = 0;
    if isfield(probStruct, 'yref') & ~isfield(sysStruct, 'Cy'),
        sysStruct.Cy = Cn;
        sysStruct.Dy = Dn;
        Cy_augm = 1;
    end
%     if isfield(probStruct, 'P_N')
%         nn = size(probStruct.P_N,1);
%         probStruct.P_N = [probStruct.P_N zeros(nn,nu); zeros(nu,nn) probStruct.R]; %terminal cost
%     end
    if isfield(sysStruct, 'Cy') & ~Cy_augm,
        if ~isfield(sysStruct, 'Dy')
            error('sysStruct.Dy must be defined if sysStruct.Cy is given!');
        end
        sysStruct.Cy = [sysStruct.Cy sysStruct.Dy];
    end
    if isfield(sysStruct,'f'),
        sysStruct.f = [sysStruct.f; zeros(nu,1)];
    else
        sysStruct.f = zeros(nx+nu,1);
    end
    f = sysStruct.f;
end

if isfield(probStruct, 'P_N')
    P = probStruct.P_N;
    [nr1, nr2] = size(probStruct.R);
    if probStruct.norm==1,
        P = P/2;
    end
    if ycost,
        if size(P,1)~=ny,
            error('In tracking, probStruct.P_N must be a matrix of the dimension of the output y!');
        end
        P = [P zeros(ny,nu); zeros(nr1, ny) probStruct.R];
    else
        if size(P,1)~=nx,
            error('In tracking, probStruct.P_N must be a matrix of the dimension of the state x!');
        end
        P = [P zeros(nx,nu); zeros(nr1,nx) probStruct.R]; %terminal cost
    end
    probStruct.P_N = P;
end

ymaxn = [sysStruct.ymax; sysStruct.umax];
yminn = [sysStruct.ymin; sysStruct.umin];
if isfield(sysStruct, 'xmax')
    % augment state constraints if present
    sysStruct.xmax = [sysStruct.xmax; sysStruct.umax];
    sysStruct.xmin = [sysStruct.xmin; sysStruct.umin];
end

%the cost is updated accordingly to punish (x-xref)' Q (x-xref) + delta_u' Rdu delta_u
%       [Q   0]
% Qn =  [0   R]

if ycost
    if(length(probStruct.Qy)==ny) & ~isfield(probStruct, 'yref')
        probStruct.Qy = [probStruct.Qy zeros(ny,nu); zeros(nu,ny) probStruct.R];
    elseif isfield(probStruct, 'yref')
        if probStruct.norm==2,
            mult = min(probStruct.Qy(probStruct.Qy~=0)/1e4);
        else
            mult = 0;
        end
        probStruct.Qy = [probStruct.Qy zeros(ny,nu); zeros(nu,ny) mult*eye(size(probStruct.R,2))];
        probStruct.yref = [probStruct.yref; zeros(nu,1)];
    end
end
Qn = [probStruct.Q zeros(nx,nu); zeros(nu,nx) probStruct.R];

if isfield(probStruct, 'Rdu'),
    probStruct.R = probStruct.Rdu;
end

%Update the outer bound Pbnd
[Hbnd,Kbnd] = double(sysStruct.Pbnd);
[nxc,nxx] = size(Hbnd);
nuc = 2*size(sysStruct.umax,1);
I2u = [eye(nu); -eye(nu)];
I2y = [eye(ny); -eye(ny)];
Hn = [Hbnd zeros(nxc,nu); zeros(nuc,nx) I2u];
Kn = [Kbnd; sysStruct.umax; -sysStruct.umin];
sysStruct.Pbnd = polytope(Hn,Kn);

if isfield(probStruct,'Tset'),
    if isfulldim(probStruct.Tset),
        [Hs,Ks]=double(probStruct.Tset);
        [nxc,nxx] = size(Hs);
        Hn = [Hs zeros(nxc,nu); zeros(nuc,nxx) I2u];
        Kn = [Ks; sysStruct.umax; -sysStruct.umin];
        probStruct.Tset=polytope(Hn,Kn);
    end
end

%write tracking data back into structure
sysStruct.A = An;
sysStruct.B = Bn;
sysStruct.C = Cn;
sysStruct.D = Dn;
sysStruct.ymax = ymaxn;
sysStruct.ymin = yminn;

if isfield(sysStruct,'dymin')
    sysStruct.dymin = [sysStruct.dymin; sysStruct.umin];
end
if isfield(sysStruct,'dymax')
    sysStruct.dymax = [sysStruct.dymax; sysStruct.umax];
end

sysStruct.umax = sysStruct.dumax;
sysStruct.umin = sysStruct.dumin;
sysStruct.dumax = Inf*ones(nu,1);
sysStruct.dumin = -Inf*ones(nu,1);
probStruct.Q = Qn;
sysStruct = rmfield(sysStruct,'verified');

if ispwa,
    sysStruct.f = f;
    sysStruct.g = g;
    sysStruct.guardX = guardX;
    sysStruct.guardU = guardU;
end
    
if isfield(sysStruct, 'InputName')
    if ~iscell(sysStruct.InputName),
        sysStruct.InputName = {sysStruct.InputName};
    end
end
if isfield(sysStruct,'StateName')
    sysStruct.StateName = {sysStruct.StateName{:}, sysStruct.InputName{:}};
end
if isfield(sysStruct,'OutputName')
    sysStruct.OutputName = {sysStruct.OutputName{:}, sysStruct.InputName{:}};
end
sysStruct.dumode = 1;

tmpopt.verbose = 0;
evalc('sysStruct = mpt_verifySysStruct(sysStruct,tmpopt);');
evalc('probStruct = mpt_verifyProbStruct(probStruct,tmpopt);');
sysStructTr = sysStruct;
probStructTr = probStruct;
