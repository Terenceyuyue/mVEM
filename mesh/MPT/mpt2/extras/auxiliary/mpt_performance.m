function [cost,XC,X0,Xinfeas] = mpt_performance(ctrl,gridpoints,Options)
%MPT_PERFORMANCE Computes performance (i.e. sum of closed-loop costs) associated to a given controller
%
% [cost,XC] = mpt_performance(ctrl,gridpoints,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Grids the state-space into given number of points and for each state
% computes the closed-loop cost index. The returning argument is a sum of
% these costs.
%
% This function can be used to compare efficiency of low-complexity controllers
% with respect to optimal control strategy.
%
% USAGE:
%   cost = mpt_performance(ctrl)
%   cost = mpt_performance(ctrl,gridpoints)
%   cost = mpt_performance(ctrl,gridpoints,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% ctrl              - Explicit controller (MPTCTRL object)
% gridpoints        - number of grid points (if not provided, 30 is default)
% Options.verbose   - Level of verbosity
% Options.Pfinal    - polytope defining part of the state-space which should
%                     be considered for cost computation. (only reasonable
%                     if ctrlStruct.Pfinal is an empty polytope)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% cost   - contains "overall cost"
% XC     - a matrix which associates each state to a corresponding cost
%          first n columns represent the state vector (x1,x2,...) and the last
%          column is the cost obtained for this state
% X0     - set of feasible initial states (states which lie in
%          ctrlStruct.Pfinal)
% Xinfeas - matrix of states for which the following holds:
%             a) they feasible initial states (i.e. x \in Pfinal)
%             b) no feasible closed-loop trajectory
%
% see also MPT_COMPUTETRAJECTORY, MPT_PLOTTIMETRAJECTORY
%

% Copyright is with the following author(s):
%
%(C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<1,
    error('mpt_performance: Wrong number of input arguments!');
end

if nargin<3,
    Options = [];
end
if ~isfield(Options,'verbose')
    Options.verbose=mptOptions.verbose;
end

if nargin>1 & ~isa(gridpoints,'double')
    error('Second input argument must be number of grid points!');
end

if isa(ctrl, 'mptctrl')
    if ~isexplicit(ctrl),
        error('This function supports only explicit controllers!');
    end
else
    if ~mpt_isValidCS(ctrl)
        error('mpt_performance: First argument has to be a valid controller structure! See mpt_control for details.');
    end
end
    
ctrlStruct = ctrl;

if nargin<2
    gridpoints = 30;
    if Options.verbose>0,
        disp('mpt_performance: Number of grid points not given, assuming 30');
    end
end

if isfield(Options,'Pfinal') %& ~isfulldim(ctrlStruct.Pfinal),
    if ~isa(Options.Pfinal,'polytope'),
        error('mpt_performance: Options.Pfinal must be a polytope object!');
    end
    if isa(ctrlStruct, 'mptctrl'),
        ctrlStruct = set(ctrlStruct, 'Pfinal', Options.Pfinal);
    elseif isstruct(ctrlStruct),
        ctrlStruct.Pfinal = Options.Pfinal;
    else
        error('mpt_performance: Unknown class of controller.');
    end
end

% first compute bounds on feasible state-space
if ~isfulldim(ctrlStruct.Pfinal),
    error('mpt_performance: Please limit the state-space of interest by Options.Pfinal !');
end

bbOptions = Options;
bbOptions.noPolyOutput = 1; % we don't need the bounding box as a polytope object
[B, lb, ub] = bounding_box(ctrlStruct.Pfinal, bbOptions);

% grid the state-space into equidistantly placed points
%dimB = dimension(B);
dimB = size(lb(:),1);
Xpoints = zeros(gridpoints,dimB);
for ii=1:dimB
    Xpoints(:,ii) = linspace(lb(ii),ub(ii),gridpoints)';
end
    

% generate all possible combinations of states
% one could use kron() here, but that one fails for high number of elements
n_states=dimB;
ZZ=[];
ZZ{n_states}=Xpoints(:,n_states);
for ii=n_states-1:-1:1,
    Zd=[];
    for jj=1:size(Xpoints,1),
        Zd=[Zd; repmat(Xpoints(jj,ii),length(ZZ{ii+1}),1)];
    end
    ZZ{ii}=Zd;
end
for ii=2:n_states,
    ZZ{ii}=repmat(ZZ{ii},length(ZZ{ii-1})/length(ZZ{ii}),1);
end
datapoints=[];
for ii=1:n_states,
    datapoints(:,ii)=ZZ{ii};
end

% datapoints now contains all possible states, but be careful, some of them do
% not belong to feasible area! 

npoints = size(datapoints,1);
locOpt = Options;
locOpt.verbose = -1;
locOpt.fastbreak = 0;
locOpt.openloop = 0;

cost = zeros(npoints,1);

% if iterative solution was computed, regions are already order in such a way,
% that the first region found is in fact the one with least step distance to the
% origin, thus we can break in mpt_getInput as soon as at least one region is
% found (results in much faster run-time!)
locOpt.fastbreak = (ctrlStruct.probStruct.subopt_lev>0);

if iscell(ctrlStruct.sysStruct.A),
    nx = size(ctrlStruct.sysStruct.A{1},2);
else
    nx = size(ctrlStruct.sysStruct.A,2);
end

narg3 = (nargout==3);
X0 = [];
Xinfeas = [];
for ii=1:npoints,
    if ii==1 | ii==npoints | mod(ii*gridpoints,gridpoints^nx)==0,
        if ii>1, ind = (ii/gridpoints)*gridpoints; else ind = 0; end
        if Options.verbose > 0
            fprintf('%d/%d         \r',ind,gridpoints^nx);
        end
    end
    x0 = datapoints(ii,:)';     % one state vector
    if narg3,
        if isinside(ctrlStruct.Pn,x0),
            X0 = [X0; x0'];
        end
    end
    [X,U,Y,D,cost(ii),traj,feasible]=mpt_computeTrajectory(ctrlStruct,x0,[],locOpt);  % compute closed-loop trajectory
    if isinf(cost(ii)) & size(X) > 1,
        % no convergence for this given point
        if Options.verbose > 1,
            fprintf('No feasible trajectory for feasible initial state %s\n', mat2str(x0(:))');
        end
        Xinfeas = [Xinfeas; x0(:)'];
    end
end

if ~isempty(Xinfeas),
    if Options.verbose > -1,
        fprintf('\nPartition is not invariant! Check failed for %d states.\n\n', size(Xinfeas,1));
    end
end

feasible_index = find(~isinf(-cost));
feasibleX = datapoints(feasible_index,:);
feasibleCost = cost(feasible_index,:);
cost = sum(feasibleCost);
XC = [feasibleX feasibleCost];
