function [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp(Matrices,Options)
%MPT_MPLP Explicitly solves the given linear program (LP)
%
% [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp(Matrices,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Multiparametric linear programming
%
% Solves the problem
%   V(x) = min H U + F x
%           U
%   s.t.   G U <= W + E x
%          bndA*x <= bndb
%
% As a solution we get 'nR' regions
%   Pn(i)={x : H x <= K}
% 
% with the optimal control law
%   U = Fi{i} x + Gi{i}
%
% and the corresponding cost function expression
%   V(x) = Bi{i} x + Ci{i}
%  
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Matrices - a struct with all the parameters which are needed.
%            See description above for explanation.
%   Matrices.H=H;
%   Matrices.G=G;   
%   Matrices.W=W;
%   Matrices.E=E;
%   Matrices.bndA=bndA;   Limits on exploration space, i.e. bndA*x<=bndb
%   Matrices.bndb=bndb;
%
% Options.mplpver      - version of mpLP solver to use (3 is default)
% Options.verbose      - level of verbosity
% Options.lpsolver     - which LP solver to use (help mpt_solveLP)
% Options.max_iter     - maximum number of iterations of the algorithm
% Options.step_size    - length of step over a facet
% Options.f_perturb    - Perturbation of the optimization direction
% Options.nu           - How many elements to extract from the optimizer (to
%                        deal with slacks) 
% Options.debug_level  
%           Due to numerical problems tiny regions are sometimes difficult to
%           calculate, i.e. are not identified at all. This may create "gaps"
%           in the computed control law. For the exploration, these will be
%           jumped over and the exploration in the state space will continue.
%           "debug_level" can have three values:
%       
%           0: No debug done
%           1: A tolerance is given to find gap in the region partition,
%              small empty regions inside the region partition will be discarded.
%              Note that this is generally not a problem, since the feedback law 
%              is continuous and can therefore be interpolated easily.
%              Correction to the calculation of the outer hull.
%           2: Zero tolerance to find gap in the region partition, empty regions
%              if they exist, will be detected, i.e. the user will be notified.
%              Correction to the calculation of the outer hull.
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% Pn,Fi,Gi           - for region Pn(i).H*x <= Pn(i).K computed input is
%                      U=Fi{i}*x+Gi{i}   
% activeConstraints  - Cell Array which stores the active constraints 
%                      of the optimizer in each region.
% Phard              - Defines the feasible state space partition (i.e. union of
%                      all regions) as Phard.H*x<=Phard.K
% details            - a structure with the following fields:
%     nR      number of regions
%     Pn      polyhedral partition
%     Fi      control law
%     Gi      control law
%     BC      connection list
%     Bi      value function
%     Ci      value function
%     nHard   number of hard constraints
%     Phard   polytope given by hard constraints
%     nb      number of constraints for each region
%     LISTa   list of active constraints
%
% see also MPT_CONSTRUCTMATRICES, MPT_MPQP, MPT_OPTCONTROL, MPT_OPTCONTROLPWA

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%     kvasnica@control.ee.ethz.ch
% (C) 2003 Mato Baotic, Automatic Control Laboratory, ETH Zurich,
%     baotic@control.ee.ethz.ch
% (C) 2002 Francesco Borrelli, Automatic Control Laboratory, ETH Zurich

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

narginchk(1,2);

global mptOptions;
if ~isstruct(mptOptions)
    mpt_error;
end
if nargin<2
    Options=[];
end

if ~isfield(Options,'mplpver')
    % 1 - original mpLP
    % 2 - improved version (requires QP solver)
    % 3 - mpLP solver introduced in MPT version 1.4.1 (requires QP solver)
    % 4 - improved version of the previous option (requires QP solver)
    % 5 - successor of solver 4, supports parameterized cost (requires QP solver)
    % 6 - latest version
    Options.mplpver=Inf;
end
if ~isfield(Options,'qpsolver')
    Options.qpsolver = mptOptions.qpsolver;
end

if Options.qpsolver==-1
    % cannot use "new" mpLP solvers without a QP solver
    Options.mplpver = 1;
end

% make sure all matrices are in full format. this is important to do, because we
% call mpt_solvelpi which doesn't convert sparse matrices to full for solvers
% which don't support them
Matrices.H  = full(Matrices.H);
Matrices.G  = full(Matrices.G);
Matrices.E  = full(Matrices.E);
Matrices.W  = full(Matrices.W);
if isfield(Matrices, 'F')
    Matrices.F  = full(Matrices.F);
end
if isfield(Matrices, 'Y')
    Matrices.Y  = full(Matrices.Y);
end
if isfield(Matrices, 'Cf')
    Matrices.Cf  = full(Matrices.Cf);
end
if isfield(Matrices, 'Cx')
    Matrices.Cx  = full(Matrices.Cx);
end
if isfield(Matrices, 'Cc')
    Matrices.Cc  = full(Matrices.Cc);
end
if isfield(Matrices, 'bndA')
    Matrices.bndA = full(Matrices.bndA);
end
if isfield(Matrices, 'bndb')
    Matrices.bndb = full(Matrices.bndb);
end
if isfield(Matrices, 'D')
    Matrices.D = full(Matrices.D);
end

switch Options.mplpver
    case 1
        [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver1(Matrices,Options);
        details.feasible = isfulldim(Pn);
    case 2
        [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver2(Matrices,Options);
        details.feasible = isfulldim(Pn);
    case 3
        [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver3(Matrices,Options);
    case 4
        [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver4(Matrices,Options);
    case 5
        [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver5(Matrices,Options);
    case 6
        [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver6(Matrices,Options);
    case {7, Inf}
        [Pn,Fi,Gi,activeConstraints,Phard,details]=mpt_mplp_ver7(Matrices,Options);
    otherwise
        error('mpt_mplp: unknown solver in Options.mplpver!');
end

if Options.mplpver >= 6
    if isfield(Options, 'nu')
        % from mpt_mplp_ver6 on, optimizers corresponding to slack
        % variables are NOT removed automatically, thus we have to do it here
        nu = Options.nu;
        for idxRegions = 1:length(Fi)
            Fi{idxRegions} = Fi{idxRegions}(1:nu,:);
            Gi{idxRegions} = Gi{idxRegions}(1:nu,:);
        end
    end
end

if ~isfield(details, 'Bi')
    details.Bi = {};
end
if ~isfield(details, 'Ci')
    details.Ci = {};
end

% handle non-zero affine cost term
if isfield(Matrices, 'F')
    F = Matrices.F(:)';
    if nnz(F)>0
        Bi = details.Bi;
        for i = 1:length(Bi)
            Bi{i} = Bi{i} + F;
        end
        details.Bi = Bi;
    end
end
