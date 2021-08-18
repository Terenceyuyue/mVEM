function ctrl = mpt_mpsol2ctrl(mpsol, arg1, arg2)
% MPSOL2CTRL Convert a solution of solvemp to MPT's controller object
%
% ctrl = mpt_mpsol2ctrl(mpsol, nu)
% ctrl = mpt_mpsol2ctrl(mpsol, nu, ny)
% ctrl = mpt_mpsol2ctrl(mpsol, sysStruct, probStruct)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Converts a solution obtained by solvemp() to a valid MPT's controller object.
% If the function is called without the "sysStruct" and "probStruct"
% arguments, a dummy LTI dynamics will be created, therefore the
% obtained controller should NOT be used in connection with the sim() and
% simplot() functions. You can, however, use it to quickly obtain a control
% input associated to a given state x0 by calling
%
%   u = ctrl(x0)
%
% You can, of course, also use the controller in Simulink. If "sysStruct"
% and "probStruct" arguments are provided, the controller can be used in
% any MPT function.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% mpsol      - solution obtained by solvemp()
% nu         - number of system inputs
% ny         - number of system outputs. if not given, assumes ny=1
% sysStruct  - system structure
% probStruct - problem structure
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% ctrl     - an explicit controller object
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

error(nargchk(2,3,nargin));

if ~iscell(mpsol),
    mpsol = {mpsol};
end

if isempty(mpsol{1}),
    error('Empty solution provided.');
end

% detect number of states
nx = size(mpsol{1}.Fi{1}, 2);

if nargin==3 & isstruct(arg1) & isstruct(arg2),
    [sysStruct, probStruct] = mpt_verifySysProb(arg1, arg2);
    
elseif nargin==3 & isa(arg1, 'double') & isa(arg2, 'double'),
    nu = arg1;
    ny = arg2;
    % create a dummy sysStruct and probStruct
    [sysStruct, probStruct] = sub_dummy_sys_prob(nx, nu, ny);
    
elseif nargin==2 & isa(arg1, 'double'),
    nu = arg1;
    ny = 1;
    % create a dummy sysStruct and probStruct
    [sysStruct, probStruct] = sub_dummy_sys_prob(nx, nu, ny);
    
else
    error('Wrong type or number of input arguments.');
    
end


% check whether we have overlaps, remove them if necessary
quadcost = sub_have_quad_cost(mpsol);

if length(mpsol) > 1 & quadcost == 0,
    % cost is linear, we can remove overlaps
    mpsol = mpt_removeOverlaps(mpsol);
    mpsol.overlaps = 0;
    
elseif length(mpsol) > 1 & quadcost == 1,
    % we can't remove overlaps, merge partitions
    mpsol = mpt_mergeCS(mpsol);
    mpsol.overlaps = 1;
    
else
    mpsol = mpsol{1};
    
end

if quadcost == 0,
    % set the quadratic term to zero matrices of appropriate size
    [mpsol.Ai{:}] = deal(zeros(nx));
end

% create a valid controller
mpsol.sysStruct = sysStruct;
mpsol.probStruct = probStruct;
mpsol.details.origSysStruct = sysStruct;
mpsol.details.origProbStruct = probStruct;
mpsol.details.runTime = 0;
ctrl = mptctrl(mpsol);


%----------------------------------------------------------------------
function [sst, pst] = sub_dummy_sys_prob(nx, nu,  ny)
% creates a dummy sysStruct and probStruct

sysStruct.A = eye(nx);
sysStruct.B = eye(nx, nu);
sysStruct.C = eye(ny, nx);
sysStruct.D = zeros(ny, nu);
sysStruct.umax = repmat(Inf, nu, 1);
sysStruct.umin = repmat(-Inf, nu, 1);
sysStruct.ymax = repmat(Inf, ny, 1);
sysStruct.ymin = repmat(-Inf, ny, 1);

probStruct.N = 1;
probStruct.Q = eye(nx);
probStruct.R = eye(nu);
probStruct.norm = 1;
probStruct.subopt_lev = 0;

[sst,pst] = mpt_verifySysProb(sysStruct, probStruct);


%----------------------------------------------------------------------
function quadcost = sub_have_quad_cost(mpsol)
% checks whether the partition has a PWQ cost

quadcost = 0;
for i = 1:length(mpsol),
    for j = 1:length(mpsol{i}.Ai),
        Ai = mpsol{i}.Ai{j};
        if ~isempty(Ai) | any(any(Ai~=0)),
            quadcost = 1;
            return
        end
    end
end
