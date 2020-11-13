function [simpleCtrl,details] = mpt_simplify(ctrl, how, Options)
%MPT_SIMPLIFY simplifies a given explicit controller by merging regions with identical control law
%
%   simplectrl = mpt_simplify(ctrl)
%   simplectrl = mpt_simplify(ctrl, how)
%   simplectrl = mpt_simplify(ctrl, Options)
%   simplectrl = mpt_simplify(ctrl, how, Options)
%   [simplectrl, details] = mpt_simplify(ctrl, how, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Simplifies a given explicit controller by merging regions which have the
% same control law. By doing so, all stability and feasibility properties
% of the controller are maintained, but complexity is greatly reduced.
%
% NOTE: Information about the value function will be lost when doing region
%       simplification!
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl             - Explicit controller (an MPTCTRL object)
% how              - which method to use for merging. allowed values are:
%                      'greedy'  - greedy merging based on heuristics 
%                                  (default setting) 
%                      'optimal' - optimal merging based on boolean
%                                  minimization
% Options.trials   - for greedy merging, defines number of trials to
%                    improve the solution (default is 1, corresponds to 1 run)
% Options.verbose  - level of verbosity {0|1|2}
% Options.nu       - number of system inputs (when simplifying partitions which
%                    are not necessarily controller objects)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% simpleCtrl       - simplified explicit controller
% details.before   - number of polytopes before merging
% details.after    - number of polytopes after merging
% details.runTime  - run time of the algorithm
% details.alg      - string, either 'greedy' or 'optimal'
%
% see also MERGE
%

% Copyright is with the following author(s):
%
% (C) 2004-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch

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

error(nargchk(1,3,nargin));

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl)
        disp('MPT_SIMPLIFY: No simplification can be performed for on-line controllers.');
        simpleCtrl = ctrl;
        details = [];
        return
    end
end

if nargin==1,
    how = 'greedy';
    Options = [];
elseif nargin==2,
    if isstruct(how)
        Options = how;
        how = 'greedy';
    elseif ischar(how)
        Options = [];
    else
        error('MPT_SIMPLIFY: Second argument must be either a string or an Options structure!');
    end
else
    if ~ischar(how)
        error('MPT_SIMPLIFY: Allowed values are ''greedy'', ''optimal'' and ''orm''.');
    end
    if ~isstruct(Options)
        error('MPT_SIMPLIFY: Third argument must be an Options structure!');
    end
end

how = lower(how);
if ~(strcmp(how, 'greedy') | strcmp(how, 'optimal') | strcmp(how, 'orm')),
    error('MPT_SIMPLIFY: Allowed methods are ''greedy'', ''optimal'' and ''orm''.');
end

if ~isfield(Options, 'merging')
    Options.merging = how;
end
if ~isfield(Options, 'verbose')
    Options.verbose = 0;
end

if ~isfield(Options, 'statusbar'),
    Options.statusbar = 0;
end

if ~isfield(Options,'trials')
    % number of attempts to improve solution of greedy merging
    Options.trials = 1;
end
if strcmpi(Options.merging,'greedy')
    Options.greedy = 1;
else
    Options.greedy = 0;
end

if isa(ctrl, 'mptctrl')
    mptctrl_input = 1;
    ctrlStruct = struct(ctrl);
else
    mptctrl_input = 0;
    ctrlStruct = ctrl;
end

if mptctrl_input,
    if ~mpt_isValidCS(ctrlStruct),
        error('mpt_simplify: input argument must be a valid controller structure!');
    elseif ctrlStruct.overlaps & ctrlStruct.probStruct.norm==2,
        error('Can''t use this function for overlapping partitions.');
    elseif ctrlStruct.overlaps,
        disp('Overlaps detected, removing overlaps...');
        % we have to call mpt_removeOverlaps() with the controller object,
        % otherwise it would not return sysStruct/probStruct fields, which
        % we need later
        ctrlStruct = mpt_removeOverlaps(ctrl);
        ctrlStruct = struct(ctrlStruct);
    end
end

if mptctrl_input,
    [nx,nu,ny,ndyn] = mpt_sysStructInfo(ctrlStruct.sysStruct);
else
    % handle partitions which are not valid controller structures
    nx = size(ctrlStruct.Fi{1},2);
    if isfield(Options, 'nu'),
        nu = Options.nu;
    else
        nu = size(ctrlStruct.Fi{1},1);
    end
end

if size(ctrlStruct.Fi{1},1)>nu,
    % if the controller was obtained by solving a CFTOC problem for LTI
    % systems, extract only control law assigned to the first step.
    % after doing so, information about the open-loop trajectory is lost.
    nR = length(ctrlStruct.Fi);
    newFi = cell(1, nR);
    newGi = cell(1, nR);
    for ir = 1:nR
        newFi{ir} = ctrlStruct.Fi{ir}(1:nu,:);
        newGi{ir} = ctrlStruct.Gi{ir}(1:nu);
    end
    ctrlStruct.Fi = newFi;
    ctrlStruct.Gi = newGi;
end

% isloate regions which have the same value of control law
sameregions = sub_uniqueOpt(ctrlStruct.Fi, ctrlStruct.Gi, nu);

newPn = mptOptions.emptypoly;

if Options.statusbar,
    Options.verbose = -1;
end

if Options.verbose > -1,
    fprintf('%d regions with %d different control laws\n',...
        length(ctrlStruct.Pn), length(sameregions.Table));
end

simpleCtrlStruct = ctrlStruct;
simpleCtrlStruct.Pn = mptOptions.emptypoly;
simpleCtrlStruct.Fi = {};
simpleCtrlStruct.Gi = {};
simpleCtrlStruct.Ai = {};
simpleCtrlStruct.Bi = {};
simpleCtrlStruct.Ci = {};
if mptctrl_input,
    simpleCtrlStruct.dynamics = [];
    simpleCtrlStruct.simplified = 1;
end

zAi = zeros(size(ctrlStruct.Ai{1}));
zBi = zeros(size(ctrlStruct.Bi{1}));
zCi = zeros(size(ctrlStruct.Ci{1}));
simpl_time = 0;

% added for [issue91]
wstatus = warning;
warning off

Options.closestatbar = 0;

if Options.greedy,
    if Options.statusbar,
        Options.status_handle = mpt_statusbar('Computing...');
    end
    
    for it = 1:length(sameregions.Table),
        
        if Options.statusbar,
            prog_min = (it - 1) / length(sameregions.Table);
            prog_max = it / length(sameregions.Table);
            Options.status_min = prog_min;
            Options.status_max = prog_max;
            if isempty(mpt_statusbar(Options.status_handle, 0, prog_min, prog_max)),
                mpt_statusbar;
                error('Break...');
            end
        end
        
        % first extract regions with identical control law
        Pn = ctrlStruct.Pn(sameregions.Table{it});

        lenPn = length(Pn);
        if lenPn==1,
            regstr = 'region';
        else
            regstr = 'regions';
        end
        if Options.verbose==0,
            fprintf('control law # %d, %d %s --> ',it,lenPn,regstr);
        end

        % try to merge regions
        [Pm,mdetails] = merge(Pn, Options);
        simpl_time = simpl_time + mdetails.runTime;

        lenPm = length(Pm);
        if lenPm==1,
            regstr = 'region';
        else
            regstr = 'regions';
        end
        if Options.verbose==0,
            fprintf('%d %s\n',lenPm,regstr);
        end

        % write data back to new controller
        simpleCtrlStruct.Pn = [simpleCtrlStruct.Pn Pm];
        for ir = 1:length(Pm),
            simpleCtrlStruct.Fi{end+1} = sameregions.Fi{it};
            simpleCtrlStruct.Gi{end+1} = sameregions.Gi{it};
        end
    end
    if Options.statusbar,
        if isempty(mpt_statusbar(Options.status_handle, 1, prog_min, prog_max)),
            mpt_statusbar;
            error('Break...');
        end
    end

else
    % optimal merging - call mpt_optMerge
    startt = cputime;
    Options.verbose = 1;
    Options.algo = 0;
    Options.color = sameregions.Reg;
    Options.Table = sameregions.Table;
    if length(ctrlStruct.Pfinal)==1,
        % feasible set is convex
        Options.PAdom = ctrlStruct.Pfinal;
        Options.PAcompl = mptOptions.emptypoly;
    else
        % feasible set is (presumably) non-convex
        dimPfinal = dimension(ctrlStruct.Pfinal);
        if dimPfinal <= 3,
            % use hull() for dimensions below 4
            Options.PAdom = hull(ctrlStruct.Pfinal);
        else
            % we have two options here how to compute the complement to a
            % non-convex Pfinal:
            %  A) use the bounding box of Pfinal
            %  B) use the envelope of Pfinal
            %
            % approach (A) is easier to compute, but (B) can actually result
            % into a convex piece (if Pfinal is convex), hence giving
            % pottentially  better run-time
            if 0,
                % use bounding box
                Options.PAdom = bounding_box(ctrlStruct.Pfinal);
            else
                % use envelope() for dimensions above 3. this is due to hull() being
                % quite slow on higher dimensions, and the result obtained by
                % envelope is good enough to obtain the domain
                %
                % NOTE! we could as well use the bounding box of Pfinal here, but
                % the advantage of envelope is that it is quite likely that the
                % envelope and Pfinal will be identical. This simplifies further
                % computation significantly.
                Options.PAdom = envelope(ctrlStruct.Pfinal);
                
                % one problem with using envelope is that the output can be a the
                % whole space, i.e. R^n. therefore we need to bound the output. we
                % do it by intersecting it with the largest box which contains
                % ctrlStruct.Pfinal
                boundP = bounding_box(ctrlStruct.Pfinal);
                
                % finally intersect the envelope with the bounding box
                Options.PAdom = Options.PAdom & boundP;
            end
        end
        Options.PAcompl = Options.PAdom \ ctrlStruct.Pfinal;
    end
    if strcmpi(Options.merging, 'orm')
        [PAmer, colorMer] = mpt_optMerge(ctrlStruct.Pn, Options);
    elseif strcmpi(Options.merging, 'optimal')
        [PAmer, colorMer] = mpt_optMergeDivCon(ctrlStruct.Pn, Options);
    else
        error('Unrecognized option "%s".', Options.merging);
    end
    simpleCtrlStruct.Pn = PAmer;
    for ii=1:length(PAmer),
        ind = colorMer.Reg(ii);
        simpleCtrlStruct.Fi{end+1} = sameregions.Fi{ind};
        simpleCtrlStruct.Gi{end+1} = sameregions.Gi{ind};
    end
    simpl_time = cputime - startt;
end

% added for [issue91]
warning(wstatus);

nRnew = length(simpleCtrlStruct.Fi);
simpleCtrlStruct.Ai = cell(1, nRnew);
simpleCtrlStruct.Bi = cell(1, nRnew);
simpleCtrlStruct.Ci = cell(1, nRnew);
[simpleCtrlStruct.Ai{:}] = deal(zAi);
[simpleCtrlStruct.Bi{:}] = deal(zBi);
[simpleCtrlStruct.Ci{:}] = deal(zCi);

if mptctrl_input,
    simpleCtrlStruct.dynamics = repmat(0,1,length(simpleCtrlStruct.Pn));
    simpleCtrlStruct.details.before_simpl = length(ctrlStruct.Pn);
    simpleCtrlStruct.details.after_simpl = length(simpleCtrlStruct.Pn);
    simpleCtrlStruct.details.simpl_time = simpl_time;
    simpleCtrlStruct.details.alg = Options.merging;
    simpleCtrlStruct.simplified = 1;
    if isfield(simpleCtrlStruct.details,'searchTree'),
        simpleCtrlStruct.details = rmfield(simpleCtrlStruct.details,'searchTree');
        disp('recomputing search tree...');
        simpleCtrlStruct = mpt_searchTree(simpleCtrlStruct);
    end
end

if Options.verbose>=0,
    fprintf('controller partition reduced to %d regions\n',length(simpleCtrlStruct.Pn));
end

details.before = length(ctrlStruct.Pn);
details.after = length(simpleCtrlStruct.Pn);
details.runTime = simpl_time;
details.alg = Options.merging;

if nargout < 2,
    clear details
end

if Options.statusbar,
    mpt_statusbar;
end

if mptctrl_input,
    simpleCtrl = mptctrl(simpleCtrlStruct);
else
    simpleCtrl = simpleCtrlStruct;
end

% assign simplified controller in caller's workspace
if ~isempty(inputname(1)) & nargout==0,
    assignin('caller',inputname(1),simpleCtrl);
end
