function [U,feasible,region]=mpt_getOptimizer(Pn,Fi,Gi,x0,Options)
%MPT_GETOPTIMIZER For a given state, extracts the optimizer obtained by a multi-parametric solver
%
% [U,feasible,region]=mpt_getOptimizer(Pn,Fi,Gi,x0,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% For the given state x0, this function extracts the optimizer $U = Fi x + Gi$
% from a solution obtained by a multi-parametric solver
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn                - polyhedral partition
% Fi, Gi            - cells containing the associated control law
% x0                - initial state
% Options.lpsolver  - Solver for LPs (see help mpt_solveLP for details)
% Options.abs_tol   - absolute tolerance
% Options.verbose   - Level of verbosity
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% U         - control input computed as U=F*x0 + G,
% region    - index of a region which contains the optimal control input associated
%             to the given state x0
% feasible  - returns 1 if  the there is at least one control law associated to
%             a given state x0, 0 otherwise
%
% see also MPT_MPLP, MPT_MPQP, MPT_GETINPUT

% Copyright is with the following author(s):
%
%(C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%         kvasnica@control.ee.ethz.ch

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

error(nargchk(4,5,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<4,
    error('mpt_getOptimizer: Wrong number of input arguments!');
end
if nargin<5,
    Options=[];
end

if ~isfield(Options,'lpsolver') % lpsolver to be used
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'abs_tol') % absolute tolerance
    Options.abs_tol=mptOptions.abs_tol;
end

x0=x0(:);

region=0;

Options.fastbreak = 1;

[isin, inwhich] = isinside(Pn,x0,Options);

if ~isin
    % no associated control law
    feasible=0;
    disp(sprintf('NO REGION FOUND FOR STATE x = [%s]',num2str(x0')));
    return
else
    % x0 lies only in one region, return the control law and cost evaluated at x0
    U = Fi{inwhich}*x0 + Gi{inwhich};
    feasible = 1;
    region=inwhich;
end
    
return