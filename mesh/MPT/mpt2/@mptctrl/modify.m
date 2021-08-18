function [ctrl,Idx_orig] = modify(ctrl, action, indices, Options)
%MODIFY Modifies an MPTCTRL object
%
%   [ctrl,Idx_orig] = modify(ctrl, action, indices)
%   [ctrl,Idx_orig] = modify(ctrl, action, Options)
%   [ctrl,Idx_orig] = modify(ctrl, action, indices, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Modifies a given MPTCTRL object by either removing and adding regions.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl        - explicit controller (an MPTCTRL object)
% action      - what to do:
%                 'remove'     - remove selected regions
%                 'pick'       - pick a subset of regions
%                 'removeflat' - remove flat regions
% indices     - indices of regions to remove/pick
% Options     - additional options structure
%   .r_small  - maximal chebychev ball size for detection of a flat region 
%               (default: mptOptions.rel_tol)
%   .bbox     - minimal bounding box size for detection of a flat region
%               (default: mptOptions.rel_tol)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% ctrl      - updated MPTCTRL object with given regions removed
% Idx_orig  - corresponding region number in the original partition ctrl
%

% Copyright is with the following author(s):
%
% (c) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (c) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(2,4,nargin));

if ~isa(ctrl, 'mptctrl')
    error('First input must be an MPTCTRL object.');
end
if ~isa(action, 'char')
    error('Second input must be a string');
end
if nargin < 3
    indices = [];
    Options = [];
elseif nargin < 4
    Options = [];
end

if isstruct(indices)
    Options = indices;
    indices = [];
end

if ~isa(indices, 'double')
    error('Third input must be a double.');
end
if ~isexplicit(ctrl)
    error('First input must be an explicit controller.');
end


switch lower(action)
    case 'remove',
        % be sure to catch cases like indices=[1 2 3 1]
        indices = unique(indices);

        % check if indices are not out of allowed range
        error(sub_checkindices(indices, length(ctrl.Pn)));
        
        % remove given regions
        Idx_orig = setdiff(1:length(ctrl.Pn),indices);
        ctrl = sub_keep(ctrl, setdiff(1:length(ctrl.Pn), indices));
               
    case 'pick',
        % be sure to catch cases like indices=[1 2 3 1]
%         indices = unique(indices);

        % check if indices are not out of allowed range
        error(sub_checkindices(indices, length(ctrl.Pn)));
        
        % pick selected regions
        ctrl = sub_keep(ctrl, indices);
        Idx_orig = indices;

    case 'removeflat',
        % remove flat regions
        [ctrl,Idx_orig] = sub_removeflat(ctrl, Options);
        
    otherwise
        error(sprintf('Unknown operation ''%s''.', action));
end

%---------------------------------------------------------------
function err = sub_checkindices(indices, lengthPn)

if any(indices > lengthPn)
    err = 'Index exceeds total number of regions';
elseif any(indices <= 0)
    err = 'Indices must be only positive.';
else
    err = '';
end



%---------------------------------------------------------------
function ctrl = sub_keep(ctrl, keep)
% keeps regions given by indices 'keep'

nkeep = length(keep);

if nkeep > 0,
    Fi = cell(1, nkeep);
    Gi = cell(1, nkeep);
    Ai = cell(1, nkeep);
    Bi = cell(1, nkeep);
    Ci = cell(1, nkeep);
    
    % just keep elements at indices 'keep'
    [Fi{:}] = deal(ctrl.Fi{keep});
    [Gi{:}] = deal(ctrl.Gi{keep});
    [Ai{:}] = deal(ctrl.Ai{keep});
    [Bi{:}] = deal(ctrl.Bi{keep});
    [Ci{:}] = deal(ctrl.Ci{keep});
    
    ctrl.Pn = ctrl.Pn(keep);
    % note that Pfinal cannot be a single polytope any longer, because we don't
    % know which polytopes we are removing. we _could_ compute a set difference
    % between Pfinal and the polytopes to be removed, but that would simply be
    % expensive
    ctrl.Pfinal = ctrl.Pn;
    ctrl.Fi = Fi;
    ctrl.Gi = Gi;
    ctrl.Ai = Ai;
    ctrl.Bi = Bi;
    ctrl.Ci = Ci;
    ctrl.dynamics = ctrl.dynamics(keep);
    ctrl.details.modified = 'removed';
    
    % now handle fields which are specific to certain control strategies. these
    % fields are mostly important for plotting purposes.
    if isfield(ctrl.details, 'regionHorizon'),
        % infinite-time solution for LTI systems
        nrh = length(ctrl.details.regionHorizon);
        if  nrh >= nkeep & all(keep <= nrh),
            ctrl.details.regionHorizon = ctrl.details.regionHorizon(keep);
        end
    end
    if isfield(ctrl.details, 'IterStore'),
        % minimum-time solution for LTI systems
        nis = length(ctrl.details.IterStore);
        if nis >= nkeep & all(keep <= nis)
            ctrl.details.IterStore = ctrl.details.IterStore(keep);
        end
    end
    if isfield(ctrl.details, 'activeConstraints'),
        % mpqp solutions
        nac = length(ctrl.details.activeConstraints);
        if nac >= nkeep & all(keep <= nac),
            ac = cell(1, nkeep);
            [ac{:}] = deal(ctrl.details.activeConstraints{keep});
            ctrl.details.activeConstraints = ac;
        end
    end
    if isfield(ctrl.details, 'Horizon'),
        % cftoc for pwa systems - open loop solution
        if iscell(ctrl.details.Horizon),
            fprintf('Warning: Open-loop solution is stored but can''t be modified. This can lead to problems.\n');
        end
    end
    
else
    % no regions to keep
    ctrl = mptctrl;
end


%---------------------------------------------------------------
function [ctrl,Idx_orig] = sub_removeflat(ctrl, Options)

global mptOptions
if ~isstruct(mptOptions),
    mpt_error;
end

if ~isfield(Options, 'r_small'),
    Options.r_small = mptOptions.rel_tol;
end
if ~isfield(Options, 'bbox'),
    Options.bbox = mptOptions.rel_tol;
end

Pn = ctrl.Pn;
[xc, Rc] = chebyball(Pn);
ind      = find(Rc<=Options.r_small);

Idx_flat = [];
for ii=ind(:)'
    if Rc(ii)<Options.r_small
        [R,l,u] = bounding_box(Pn(ii),struct('noPolyOutput',1));
        d = abs(u-l);
        if any(d>Options.bbox)
            Idx_flat = [Idx_flat ii];
        end
    end
end

% remove flat regions from the controller
Idx_orig = setdiff(1:length(ctrl.Pn), Idx_flat);
ctrl = sub_keep(ctrl, Idx_orig);
