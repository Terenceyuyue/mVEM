function X0 = mpt_feasibleStates(ctrl,gridpoints,Options)
%MPT_FEASIBLESTATES returns equidistantly spaced data points in feasible set
%
% X0 = mpt_feasibleStates(ctrl, gridpoints, Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Grids the state-space into given number of points and returns those which
% lie in the feasible set of a given explicit controller
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl            - Explicit controller (MPTCTRL object)
% gridpoints      - number of grid points (if not provided, 30 is default)
% Options.Pfinal  - polytope defining part of the state-space which should
%                   be considered for cost computation. (only reasonable
%                   if ctrl.Pfinal is an empty polytope)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% X0     - set of feasible initial states (states contained in ctrl.Pn)
%
% see also MPT_PERFORMANCE, POLYTOPE/GRID
%

% Copyright is with the following author(s):
%
%(C) 2004-2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(1,3,nargin));
if nargin < 3,
    Options = [];
end

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl)
        error('This function supports only explicit controllers!');
    end
    x0format = ctrl.details.x0format;
    nx = x0format.required - x0format.reference - x0format.uprev;
    
elseif isstruct(ctrl),
    if mpt_isValidCS(ctrl)
        nx = mpt_sysStructInfo(ctrl.details.origSysStruct);
    else
        error('First argument has to be a valid controller structure!');
    end
    
else
    error('Unknown type of first input argument.');
    
end

if nargin>1 & ~isa(gridpoints,'double')
    error('Second input argument must be number of grid points!');
end

if nargin < 2
    gridpoints = 30;
end

if isfield(Options, 'Pfinal'),
    if ~isa(Options.Pfinal, 'polytope'),
        error('Options.Pfinal must be a polytope object!');
    end
    if dimension(Options.Pfinal) ~= nx,
        error('Wrong dimension of Options.Pfinal');
    end
    Pfinal = Options.Pfinal;
elseif isfulldim(ctrl.Pfinal)
    Pfinal = ctrl.Pfinal;
else
    Pfinal = ctrl.Pn;
end

if dimension(Pfinal) > nx,
    Pfinal = projection(Pfinal, 1:nx);
end

X0 = grid(Pfinal, gridpoints);
