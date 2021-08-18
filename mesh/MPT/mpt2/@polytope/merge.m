function [Pm,details] = merge(Pn,Options)
%MERGE merges polytopes together
%
% [Pm, details] = merge(Pn, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Uses either greedy merging or optimal merging to join regions.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn  - polytope array
% Options.greedy    - if set to 1 (default), uses greedy merging. If set to
%                     0, optimal merging will be used (may be slow, requires
%                     espresso solver to be installed)
% Options.Pn_convex - if Pn is convex this should be set to 1 to avoid slowdown
%                     and unnecessary increase of complexity. (default 0)
% Options.verbose   - level of verbosity (0/1/2)
% Options.trials    - for greedy merging, defines number of trials to
%                     improve the solution (default is 1, corresponds to 1 run)
% Options.use_reduceunion - if set to 1 (default is 0) regions which are
%                           completely covered by other regions are kicked
%                           out using the reduceunion() function. this can
%                           have both positive and negative implication on
%                           the speed of merging! usually it makes sense to
%                           activate this option when calculating invariant
%                           subsets where a large number of small regions
%                           are generated.
%
% NOTE!!! if optimal merging is used (Options.greedy=0), then the result might
% contain overlaps! if you do not want them, you must additionaly specify
% Options.algo, see 'help mpt_optMerge' for more details.
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Pm               - merged polytopes
% details.before   - number of polytopes before merging
% details.after    - number of polytopes after merging
% details.runTime  - run time of the algorithm
% details.alg      - string, either 'greedy', 'optimal', or 'reduceunion'
%
% see also MPT_GREEDYMERGING, MPT_OPTMERGE

% Copyright is with the following author(s):
%
% (C) 2007 Michal Kvasnica, Slovak University of Technology in Bratislava
%          michal.kvasnica@stuba.sk
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (C) 2005 Tobias Geyer, Automatic Control Laboratory, ETH Zurich,
%          geyer@control.ee.ethz.ch
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

if ~isa(Pn, 'polytope'),
    error('ISINSIDE: input argument MUST be a polytope');
end

if isempty(Pn.Array),
    % nothing to merge, input is single polytope
    Pm = Pn;
    if nargout > 1,
        details.before = 1;
        details.after = 1;
        details.runTime = 0;
        details.alg = 'merge';
    end
    return
end

if nargin<2,
    Options = [];
end

Options = mpt_defaultOptions(Options, ...
    'greedy', 1, ... % set this option to 0 for Optimal merging 
    'verbose', 0, ...
    'trials', 1, ... % number of attempts to improve solution of greedy merging
    'Pn_convex', 0, ... % set to 1 if you know the array is convex (speeds up computation)
    'use_reduceunion', 0, ...
    'overlaps', false ); 

if Options.use_reduceunion
    % try to use reduceunion() first to kick out regions which are completely
    % covered by other regions.
    % NOTE: this only helps if Pn is overlapping! luckily this is mostly the
    %       case when calculating invariant subsets
    startt = clock;
    Pm = reduceunion(Pn);
    runtime = etime(clock, startt);
    if isempty(Pm.Array)
        if nargout > 1,
            details.before = length(Pn);
            details.after = 1;
            details.runTime = runtime;
            details.alg = 'reduceunion';
        end
        return
    else
        Pn = Pm;
    end
end
    

if Options.greedy,
    [Pm, details] = mpt_greedyMerging(Pn, Options);
else
    if ~Options.Pn_convex
        % compute domain and complement
        % Note: This is slow. Usually, this step should be controlled by the
        % user, since in most cases, the set of polyhedra is convex and thus
        % this computation is not needed. For more details see the help in
        % mpt_optMerge. T. Geyer July 20, 2005
        try
            Options.PAdom = hull(Pn);
        catch
            E = pelemfun(@extreme, Pn);
            V = cat(1, E{:});
            V = round(V*1e10)/1e10;
            Options.PAdom = polytope(V);
        end

        Pn_diff = Options.PAdom \ Pn;
        if ~isfulldim(Pn_diff)
            Options.PAcompl = [];
        else
            Options.PAcompl = Pn_diff;
        end
    else
        Options.PAcompl = [];
    end
    if Options.overlaps
        Options.algo = 0;
    else
        Options.algo = 2; % non-overlapping regions please (set .algo=0 for overlaps)
    end
    startt = clock;
    [Pm] = mpt_optMerge(Pn, Options);
    details.runTime = etime(clock, startt);
    details.before = length(Pn);
    details.after = length(Pm);
    details.alg = 'mpt_optMerge';
end

if nargout < 2,
    clear details
end
