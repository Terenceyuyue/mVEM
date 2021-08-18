function [P]=mldivide(P1,P2,Options)
%MLDIVIDE Set difference
%
% Set difference
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes set difference between P1 and P2
%
% P1 and P2 can be both single polytopes and polytope arrays
%
% USAGE:
%   R=P\Q
%   R=mldivide(P,Q,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P,Q                 - Polytopes or arrays of polytopes
% Options.abs_tol     - absolute tolerance
% Options.lpsolver    - Which LP solver to use (see help mpt_solveLP)
% Options.simplecheck - If set to 1, function aborts as soon as it detects
%                       that P1\P2 is non-empty (answers the question: is P1
%                       fully covered by P2? answer is true if P1\P2 is empty set)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R   - polytope (or polytope array) describing the set difference
%
% see also REGIONDIFF
%

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
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


if ~isa(P1,'polytope') | ~isa(P2,'polytope')
    error('MLDIVIDE: input arguments MUST be polytopes!');
end

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3,
    Options=[];
end

if ~isfield(Options,'lpsolver')
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'simplecheck')
    Options.simplecheck=0;
end

if Options.simplecheck,              
    % for fast coverage detection, will not construct the whole set
    % just answers the question: is P1 fully covered by P2?
    % the answer is true if P1\P2 is empty set
    Options.reduce=0;
    Options.reduce_output=0;
    Options.constructset=0;
else
    % construct the full solution
    Options.reduce=1;
    Options.reduce_output=1;
    Options.constructset=1;
end

P=polytope;
R=[];
if length(P1.Array)>0,
    for ii=1:length(P1.Array),
        % if P1 is a polyarray, regiondiff has to be called for each element of P1
        Pdiff=regiondiff(P1.Array{ii},P2,Options);
        P=[P Pdiff];  
        if ~Options.constructset & isfulldim(Pdiff),   % if P1 is not covered with P2, abort
            return
        end
    end
else
    P=regiondiff(P1,P2,Options);
end
