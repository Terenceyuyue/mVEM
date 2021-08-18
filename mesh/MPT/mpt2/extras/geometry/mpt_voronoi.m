function [Pn]=mpt_voronoi(points,Options)
%MPT_VORONOI Computes the voronoi diagram via mpLP
%
% [Pn]=mpt_voronoi(points,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% The voronoi diagram is a partition of the state space; For a given set of
% points pj, each region Pn(j) is defined as 
%           Pn(j)={x \in R^n | d(x,pj)<=d(x,pi), \forall i \neq j}
%  
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% points         -  Optional input:
%                   Matrix p times nx of points: nx is state space dimension and
%                   p is the number of points
%                   The entry is graphical in 2D if no parameters are passed.
% Options.pbound -  A "bounding polytope". If provided, the voronoi cells will
%                   be bounded by this polytope. If not provided, the cells will
%                   be bounded by a hypercube as big as 1.5x the maximum
%                   coordinate of any of the seed points
% Options.plot   -  If set to 1, plots the voronoi diagram (0 is default)
% Options.sortcells - If set to 1, resulting Voronoi partition will be ordered
%                     in a way such that Pn(i) corresponds to seed point i.
%                     (Default is 1)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% Pn            -   Voronoi partition
%
% see also MPT_DELAUNAY
%

% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%     grieder@control.ee.ethz.ch

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

if ~isstruct(mptOptions)
    mpt_error;
end

if nargin < 2
    Options = [];
end
if ~isfield(Options, 'verbose')
    % this is to keep mpt_mplp silent
    Options.verbose = -1;
end
if ~isfield(Options, 'sortcells')
    Options.sortcells = 1;
end

% if(isa(points,'polytope'))
%      points = extreme(points);
% end
if(~isfield(Options,'abs_tol'))
    Options.abs_tol = mptOptions.abs_tol;
end

nx = size(points,2);
npoints = size(points, 1);

% perturb points to obtain general position
% points = points+(rand(nopoints,nx)-0.5)*Options.abs_tol;

%build voronoi constraint matrices
Matrices.G = -ones(npoints,1);
Matrices.W = sum(points.^2,2);
Matrices.E = -2*points;
Matrices.H = 1;

if isfield(Options, 'pbound')
    % user has provided a bounding polytope, check if the input is correct
    pbound = Options.pbound;
    if ~isa(pbound, 'polytope')
        error('''Options.pbound'' must be a polytope object.');
    elseif dimension(pbound) ~= nx
        error(sprintf('''Options.pbounds'' must be a polytope in %dD.', nx)); %#ok<SPERR>
    elseif length(pbound)>1
        error('''Options.pbound'' must be a single polytope.');
    end
else
    % by default we bound the voronoi cells with a hypercube of size 1.5x bigger
    % than the maximum coordinate of seeds
    pbound = unitbox(nx, max(max(points))*1.5);
end
[Matrices.bndA, Matrices.bndb] = double(pbound);

%solve mpLP
Pn = mpt_mplp(Matrices,Options);

% re-order regions such that Pn(i) corresponds to seed "i"
Idx = [];
for i = 1:npoints
    % find which region corresponds to seed "i"
    [isin, inwhich] = isinside(Pn, points(i, :)');
    if isin
        % actually it shouls never happen that a seed does not belong to any
        % polytope, but double-check that
        Idx = [Idx,inwhich(1)];
    else
        warning(sprintf('MPT_VORONOI: point %d does not belong to any polytope!', i)); %#ok<SPWRN>
    end
end
% only re-order polytopes if all regions have an associated seed point
if length(Idx) == npoints
    Pn = Pn(Idx);
else
    warning('MPT_VORONOI: returning unordered partition.');
end
