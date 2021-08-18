function sys = char(ctrl)
%CHAR Provides textual description of a given controller
%
% sys = char(ctrl)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Provides textual description of a given controller
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl - MPT controller object
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% sys  - String containing the information
%
% see also MPTCTRL
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

sys = '';

if isempty(ctrl.sysStruct)
    sys = 'Empty controller object';
    return
end
dim = dimension(ctrl.Pn);
nR = length(ctrl.Pn);

optlev = ctrl.probStruct.subopt_lev;
N = ctrl.probStruct.N;
sysStruct = ctrl.sysStruct;
probStruct = ctrl.probStruct;

if optlev==0 & isinf(N)
    ctrltype = 'infinite-horizon optimal control problem';
elseif optlev==0 & ~isinf(N)
    ctrltype = sprintf('horizon %d optimal control problem', N);
elseif optlev==1
    ctrltype = 'minimum-time problem';
else
    ctrltype = sprintf('low-complexity solution (Horizon %d)', N);
end
   
if mpt_isnoise(sysStruct.noise) | isfield(sysStruct, 'Aunc')
    ctrltype = [ctrltype ' robust'];
end
cost_str = sprintf('%s (%s-norm)', ctrltype, num2str(probStruct.norm));


%========================================================================
% determine if the controller is stabilizable
if isexplicit(ctrl)
    stability = isstabilizable(ctrl);
    if stability
        stability_str = 'the closed-loop system is stable';
    else
        stability_str = 'no stability/convergence guarantees given';
    end
end

%========================================================================
% determine if the controller is invariant
if isexplicit(ctrl),
    invariance = isinvariant(ctrl);
    if invariance==1
        invariance_str = 'the closed-loop system is invariant';
    elseif invariance==-1
        invariance_str = 'no invariance statement for on-line solutions';
    else
        invariance_str = 'the closed-loop system may be not invariant';
    end
end

%========================================================================
% determine regulation objective
regulation_str = '';
if probStruct.Tconstraint==2 & isfulldim(probStruct.Tset)
    regulation_str = 'towards a given target set';
elseif probStruct.tracking>0 & isfield(probStruct, 'Qy')
    regulation_str = 'towards free output reference';
elseif probStruct.tracking>0 & ~isfield(probStruct, 'Qy')
    regulation_str = 'towards free state reference';
elseif isfield(probStruct, 'Qy'),
    regulation_str = 'towards zero output';
    if isfield(probStruct, 'yref')
        if any(probStruct.yref ~= 0)
            regulation_str = 'towards given output reference';
        end
    end
elseif isfield(probStruct, 'xref')
    regulation_str = 'towards given state reference';
else
    regulation_str = 'towards origin';
end
if probStruct.tracking==1
    regulation_str = [regulation_str ' (offset-free)'];
elseif probStruct.tracking==2
    regulation_str = [regulation_str ' (with possible offset)'];
end

sysStruct = ctrl.details.origSysStruct;

%========================================================================
% determine type of controlled system
[nx,nu,ny,ndyn] = mpt_sysStructInfo(ctrl.details.origSysStruct);
if iscell(sysStruct.A),
    if isfield(sysStruct, 'nonlinhandle'),
        % piecewise nonlinear system
        systype = sprintf('Piecewise nonlinear system with %d dynamics', ndyn);
    else
        % PWA system
        systype = sprintf('PWA system with %d dynamics', ndyn);
    end
elseif isfield(sysStruct, 'nonlinhandle'),
    systype = 'Nonlinear system';
else
    systype = 'LTI system';
end

% make plurals, i.e. if nx>1 -> 'states', otherwise 'state'
plurx = ['state' repmat('s', 1, (nx>1))];
pluru = ['input' repmat('s', 1, (nu>1))];
plury = ['output' repmat('s', 1, (ny>1))];
sysinfo_str = sprintf('%s, %d %s, %d %s, %d %s', systype, ...
    nx, plurx, nu, pluru, ny, plury);


%========================================================================
% determine dimension of input vector

if probStruct.tracking>0,
    sysStruct = ctrl.sysStruct;

    if isfield(sysStruct, 'dims'),
        nx = sysStruct.dims.nx;
        nu = sysStruct.dims.nu;
        ny = sysStruct.dims.ny;
    else
        [nx, nu, ny] = mpt_sysStructInfo(sysStruct);
    end
    
    % make plurals, i.e. if nx>1 -> 'states', otherwise 'state'
    plurx = ['state' repmat('s', 1, (nx>1))];
    pluru = ['input' repmat('s', 1, (nu>1))];
    plury = ['output' repmat('s', 1, (ny>1))];

    if probStruct.tracking==1,
        if isfield(probStruct, 'Qy');
            statedim = sprintf('%d system %s + %d past %s + %d reference %s', ...
                nx, plurx, nu, pluru, ny, plury);
        else
            statedim = sprintf('%d system %s + %d past %s + %d reference %s', ...
                nx, plurx, nu, pluru, nx, plurx);
        end
    else
        if isfield(probStruct, 'Qy');
            statedim = sprintf('%d system %s + %d reference %s', ...
                nx, plurx, ny, plury);
        else
            statedim = sprintf('%d system %s + %d reference %s', ...
                nx, plurx, nx, plurx);
        end
    end
else
    [nx,nu,ny] = mpt_sysStructInfo(sysStruct);
    % make plurals, i.e. if nx>1 -> 'states', otherwise 'state'
    plurx = ['state' repmat('s', 1, (nx>1))];
    pluru = ['input' repmat('s', 1, (nu>1))];
    pluru = ['output' repmat('s', 1, (ny>1))];

    if any(~isinf(sysStruct.dumin)) | any(~isinf(sysStruct.dumax))
        % we have deltaU formulation
        statedim = sprintf('%d system %s + %d past %s', nx, plurx, nu, pluru);
    else
        statedim = sprintf('%d system %s', nx, plurx);
    end
end

sys = '';
if isexplicit(ctrl)
    sys = strvcat(sys, sprintf('Explicit MPC controller defined over %d regions in %dD', nR, dim));
else
    sys = strvcat(sys, sprintf('On-line MPC controller'));
end
sys = strvcat(sys, ' ');

sys = strvcat(sys, sprintf('   Controlled system: %s', sysinfo_str));
sys = strvcat(sys, sprintf('   Control objective: %s', cost_str));
sys = strvcat(sys, sprintf('          Regulation: %s', regulation_str));

if isexplicit(ctrl)
    sys = strvcat(sys, sprintf('          Invariance: %s', invariance_str));
    sys = strvcat(sys, sprintf('           Stability: %s', stability_str));
end
    
return
