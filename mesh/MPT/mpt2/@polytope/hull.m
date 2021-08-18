function [P,Pn]=hull(Pn,Options)
%HULL Convex hull of n polytopes
% 
% [P,Pn]=hull(Pn,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes a convex hull of n polytopes given by the polytope array Pn
%
% USAGE:
%   P = hull([P1 P2 P3])
%   P = hull(Pn,Options)
%   [P,Pn] = hull(Pn)
%   [P,Pn] = hull(Pn,Options)
%
% NOTE:
%   Convex hull requires computing extreme points, what can be expensive.
%   Thats why also the second output argument should be always retrieved,
%   since once the extreme points are computed for the input argument,
%   they can be stored for future faster calculation
%   (see help extreme)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn                      - Polytope array
% Options.extreme_solver  - Which method to use for vertex enumeration
%                           (0 - brute force enumeration, 3 - CDD)
%                           (see help mpt_init)
% Options.verbose         - Level of verbosity
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% P     - convex hull of polytopes stored in Pn
% Pn    - extreme points have been computed and stored for the input argument
%
% see also ENVELOPE, UNION, EXTREME
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%          baotic@control.ee.ethz.ch

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


global mptOptions;

if nargin<2,
    Options=[];
end

if ~isstruct(mptOptions),
    mpt_error;
end
if ~isfield(Options,'extreme_solver'),
    Options.extreme_solver=mptOptions.extreme_solver;
end
if ~isfield(Options, 'verbose'),
    Options.verbose = mptOptions.verbose;
end

V=[];

% if ~isa(Pn,'polytope'),
%     error('HULL: Input argument MUST be polytopes!');
% end

if length(Pn.Array)>0,           % cycle through all elements of Pn
    for jj=1:length(Pn.Array),
        [Vext,R,Pn.Array{jj}]=extreme(Pn.Array{jj},Options);   % compute extreme points of Pn(jj)
        if ~isempty(R),
            error('HULL: Unbounded polytope detected, polyhedra are not allowed!');
        end
        V=[V; Vext];         % merge vertices together
    end
else
    [Vext,R,Pn]=extreme(Pn,Options);   % Pn is just a single polytope
    if ~isempty(R),
        error('HULL: Unbounded polytope detected, polyhedra are not allowed!');
    end
    V=[V; Vext];
end

extsolvers = mptOptions.solvers.extreme(:);
% remove current solver
extsolvers(find(extsolvers == Options.extreme_solver)) = [];
% put current solver to the front
extsolvers = [Options.extreme_solver; extsolvers]';

Pnbounded = isbounded(Pn);

for extsol = extsolvers,
    Options.extreme_solver = extsol;
    P=hull(V,Options);    % convex hull of all points
    failed = 0;
    if Pnbounded,
        if ~isbounded(P),
            % hull must be bounded if input is bounded
            failed = 1;
        end
    end
    if failed == 0,
        % no problems due to unboundeness, exit
        return
    elseif Options.verbose > 1,
        fprintf('POLYTOPE/HULL: "%s" failed, trying another method...\n', mpt_solverInfo('extreme', extsol));
    end
end
fprintf('POLYTOPE/HULL: All methods failed, hull is unbounded!\n');
