function [sys,x0,str,ts] = mpt_simSys_sfunc(t,x,u,flag,sysStruct,X0,Ts,nx,nu,ny,constr)
%MPT_SIMSYS_SFUNC S-function to simulate sysStruct system in Simulink
%
% [sys,x0,str,ts] = mpt_simSys_sfunc(t,x,u,flag,sysStruct,X0,Ts)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Allows to simulate a system described in sysStruct in Simulink, providing
% external input U
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% t           - time
% x           - state
% u           - input
% flag        - S-function flag
% sysStruct   - system structure
% X0          - initial condition
% Ts          - sampling time
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sys         - state update
% x0          - initial condition
% ts          - sampling time
%
% see also MPT_SIMSYS

% Copyright is with the following author(s):
%
%(C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%              kvasnica@control.ee.ethz.ch

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

if ~isstruct(sysStruct)
    error('sysStruct must be a structure!');
end

% check if we have direct feedthrough of inputs
directu = 0;

% we need direct feedthrough in case of:
%  - non-zero D matrix
%  - guards on inputs
%  - MLD systems

if isfield(sysStruct, 'data'),
    if isfield(sysStruct.data, 'MLD'),
        if isfield(sysStruct.data.MLD, 'D1'),
            D1 = sysStruct.data.MLD.D1;
            D2 = sysStruct.data.MLD.D2;
            D3 = sysStruct.data.MLD.D3;
            directu = any(any(D1~=0)) | any(any(D2~=0)) | any(any(D3~=0));
        end
    end
else
    if iscell(sysStruct.D),
        for ii = 1:length(sysStruct.D);
            if any(any(sysStruct.D{ii}~=0)) | any(any(sysStruct.guardU{ii}~=0)),
                directu = 1;
                break
            end
        end
    else
        directu = any(any(sysStruct.D~=0));
    end
end

switch flag,
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts] = mdlInitializeSizes(sysStruct,X0,Ts,nx,nu,ny,directu);
        
        %%%%%%%%%%
        % Update %
        %%%%%%%%%%
    case 2,
        sys = mdlUpdate(t,x,u,sysStruct,nx,nu,ny,constr); 
        
        %%%%%%%%%%
        % Output %
        %%%%%%%%%%
    case 3,
        sys = mdlOutputs(t,x,u,sysStruct,nx,nu,ny,directu,constr);
        
        %%%%%%%%%%%%%
        % Terminate %
        %%%%%%%%%%%%%
    case 9,
        sys = []; % do nothing
        
        %%%%%%%%%%%%%%%%%%%%
        % Unexpected flags %
        %%%%%%%%%%%%%%%%%%%%
    otherwise
        error(['unhandled flag = ',num2str(flag)]);
end

%end dsfunc

%
%=======================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=======================================================================
%
function [sys,x0,str,ts] = mdlInitializeSizes(sysStruct,X0,Ts,nx,nu,ny,directu)

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = nx+nu;
sizes.NumOutputs     = nx+ny;
sizes.NumInputs      = nu;
sizes.DirFeedthrough = directu;
sizes.NumSampleTimes = 1;

sys = simsizes(sizes);

if isempty(X0),
    x0 = zeros(sizes.NumDiscStates,1);
else
    x0 = [X0(:); zeros(nu, 1)];
end
str = [];
ts  = [Ts 0]; 

% end mdlInitializeSizes

%
%=======================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=======================================================================
%
function sys = mdlUpdate(t,x,u,sysStruct,nx,nu,ny,constr)

u = u(:);
if constr & isfield(sysStruct, 'umax'),
    uviol = find(u > sysStruct.umax);
    u(uviol) = sysStruct.umax(uviol);
    uviol = find(u < sysStruct.umin);
    u(uviol) = sysStruct.umin(uviol);
end

[X,U,Y] = mpt_simSys(sysStruct,x(1:nx),u(:));
X = X(:);
if constr & isfield(sysStruct, 'xmax')
    % enforce constraints
    xviol = find(X > sysStruct.xmax);
    X(xviol) = sysStruct.xmax(xviol);
    xviol = find(X < sysStruct.xmin);
    X(xviol) = sysStruct.xmin(xviol);
end
  
sys = [X(:); u(:)];
% disp(mat2str([t X' Y' U']))

%end mdlUpdate

%
%=======================================================================
% mdlOutputs
% Return the output vector for the S-function
%=======================================================================
%
function sys = mdlOutputs(t,x,u,sysStruct,nx,nu,ny,directu,constr)

% if t == 0
%     C = sysStruct.C;
%     if iscell(C)
%         C = C{1};
%     end
%     sys = [x(1:nx); C*x(1:nx)];
%     return
% end
t
if directu,
    % use direct feedthrough
    u = u(:);
    if constr & isfield(sysStruct, 'umax'),
        uviol = find(u > sysStruct.umax);
        u(uviol) = sysStruct.umax(uviol);
        uviol = find(u < sysStruct.umin);
        u(uviol) = sysStruct.umin(uviol);
    end
else
    % otherwise use the same U which was used for state update
    u = x(nx+1:end);
end

[X,U,Y] = mpt_simSys(sysStruct,x(1:nx),u);
Y = Y(:);
if constr & isfield(sysStruct, 'ymax')
    % enforce constraints
    yviol = find(Y > sysStruct.ymax);
    Y(yviol) = sysStruct.ymax(yviol);
    yviol = find(Y < sysStruct.ymin);
    Y(yviol) = sysStruct.ymin(yviol);
end

sys = [x(1:nx); Y];

%end mdlOutputs
