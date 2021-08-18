function invCtrl = mpt_invariantSet(ctrl, Options)
%MPT_INVARIANTSET Computes (robust) positive invariant subset of an explicit controller
%
%   invCtrl = mpt_invariantSet(ctrl)
%   invCtrl = mpt_invariantSet(ctrl, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes the maximal (robust) positive invariant set of a given explicit
% controller.
%
% NOTE: This function is NOT available for on-line controllers.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl              - explicit controller (an MPTCTRL object)
% Options.verbose   - Level of verbosity {0|1|2}
% Options.nohull    - If set to 1, do not compute convex unions (may
%                     significantly prolong run-time)
% Options.maxIter   - maximum number of iterations. Set is not invariant if
%                     iteration is aborted prior to convergence (default is 200)
% Options.useTmap   - If set to true (default), transition map will be
%                     computed to rule out certain transitions
% Options.mergefinal - If set to true (default is false), tries to simplify
%                      final result by merging regions
% Options.sphratio  - Gives factor which governs maximum number of separating
%                     hyperplanes computed in transition maps. Number of
%                     separating  hyperplnaes computed at each step is given by
%                     length(Pn)*length(targetPn) / Options.ratio
%                     Default value is 20.
%                     Set this option to 0 if you don't want to impose any limit
%                     on number of separating hyperplanes.
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% invCtrl   - invariant subset of a given explicit controller
%
% see also MPT_INFSETPWA
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

global mptOptions
if ~isstruct(mptOptions)
    mpt_error;
end

if ~isa(ctrl, 'mptctrl')
    error('MPT_INVARIANTSET: First input argument must be a valid MPTCTRL object!');
end

if nargin<2,
    Options = [];
elseif ~isstruct(Options)
    error('MPT_INVARIANTSET: Second input argument must be an Options structure!');
end

if ~isfield(Options, 'verbose')
    Options.verbose = mptOptions.verbose;
end
if ~isfield(Options, 'abs_tol')
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options, 'useTmap')
    Options.useTmap = 1;
end
if ~isfield(Options, 'mergefinal'),
    Options.mergefinal = 0;
end

if ~isexplicit(ctrl)
    error('MPT_INVARIANTSET: This function is not available for on-line controllers!');
end

if isinvariant(ctrl)
    % exit if controller is already invariant
    if Options.verbose>0
        disp('Controller is already invariant...');
    end
end

sysStruct = ctrl.sysStruct;
probStruct = ctrl.probStruct;
Wnoise = sysStruct.noise;

nR = length(ctrl);
Pn = ctrl.Pn;
Fi = ctrl.Fi;
Gi = ctrl.Gi;
Acell=cell(1, nR);
Fcell=cell(1, nR);
if isfield(probStruct, 'FBgain'),
    % handle pre-stabilization with feedback
    FBgain = probStruct.FBgain;
else
    FBgain = zeros(size(Fi{1}));
end

if (iscell(sysStruct.A))
    if Options.verbose>0,
        disp('Linking dynamics to regions...');
    end
    nu = size(sysStruct.B{1}, 2);
    FBgain = FBgain(1:nu, :);    
    for ii=1:nR
        [x,R] = chebyball(Pn(ii));            % compute center of the chebyshev's ball
        Fi{ii} = Fi{ii}(1:nu, :);
        Gi{ii} = Gi{ii}(1:nu);
        for jj=1:length(sysStruct.A)          % go through all dynamics description
            if max(sysStruct.guardX{jj}*x+sysStruct.guardU{jj}*((Fi{ii}+FBgain)*x+Gi{ii})-sysStruct.guardC{jj})<Options.abs_tol,    % check which dynamics is active in the region
                Acell{ii}=sysStruct.A{jj}+sysStruct.B{jj}*(Fi{ii} + FBgain);
                Fcell{ii}=sysStruct.B{jj}*Gi{ii}+sysStruct.f{jj};
            end
        end
    end
    if Options.verbose > 0,
        disp('Linking finished, now computing invariant set...');
    end
else
    hasaffine = isfield(sysStruct, 'f');
    nu = ctrl.details.dims.nu;
    for ii=1:length(Pn)
        Acell{ii} = sysStruct.A + sysStruct.B*(Fi{ii}(1:nu,:) + FBgain(1:nu,:));
        Fcell{ii} = sysStruct.B*Gi{ii}(1:nu);
        if hasaffine,
            Fcell{ii} = Fcell{ii} + sysStruct.f;
        end
    end
end

% call mpt_infsetPWA to compute invariant subset
[Pninv,dynamics]=mpt_infsetPWA(Pn, Acell, Fcell, Wnoise, Options);

if ~isfulldim(Pninv) | isempty(dynamics)
    % no invariant set exists
    if Options.verbose > 0,
        fprintf('Invariant set is empty\n');
    end
    invCtrl = mptctrl;
else
    ctrl = struct(ctrl);
    invCtrl = ctrl;
    invCtrl.Pn = Pninv;
    invCtrl.Pfinal = Pninv;
    invCtrl.Fi = {ctrl.Fi{dynamics}};
    invCtrl.Gi = {ctrl.Gi{dynamics}};
    invCtrl.Ai = {ctrl.Ai{dynamics}};
    invCtrl.Bi = {ctrl.Bi{dynamics}};
    invCtrl.Ci = {ctrl.Ci{dynamics}};
    invCtrl.dynamics = ctrl.dynamics(dynamics);
    invCtrl.details.isinvariant = 1;
    invCtrl.simplified = 0;
    invCtrl = mptctrl(invCtrl);
end
