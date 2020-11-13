function status = le(P,Q,Options)
%LE Checks if polytope P is a subset of polytope Q
%
% status = le(P,Q,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = LE(P,Q) returns TRUE (1) if P<=Q (i.e. P is a subset of Q)
%
% USAGE:
%   P<=Q
%   le(P,Q)
%   le(P,Q,Options)
%
% NOTE:
%   If P is an empty polytope, the statement is always TRUE
%   If Q is an empty polytope, the statement is always FALSE
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P,Q              - Polytopes
% Options.rel_tol  - relative tolerance
% Options.abs_tol  - absolute tolerance
% Options.lpsolver - LP solver to use (see help mpt_solveLP)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status           - Logical statement (true if P<=Q, false otherwise)
%
% see also GE, EQ, NE, LT, GT
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

if ~isa(P, 'polytope') | ~isa(Q, 'polytope')
  error('LE: Arguments MUST be a polytope object');
end

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3
    Options=[];
end

Options.rel_tol=mptOptions.rel_tol;   % relative tolerance
Options.abs_tol=mptOptions.abs_tol;   % absolute tolerance
Options.lpsolver=mptOptions.lpsolver; % default LP solver
if ~isfield(Options,'elementwise')
    Options.elementwise=0;
end

lenP=length(P.Array);
lenQ=length(Q.Array);

if ~isfulldim(P,Options)
    % if P is empty, statement is true
    status = 1;
    return
end

if ~isfulldim(Q,Options)
    % if Q is empty, statement is false
    status = 0;
    return
end

% check dimensions
dimP = dimension(P);
dimQ = dimension(Q);
if dimP ~= dimQ,
    error('LE: Only polytopes of equal dimensionality can be compared.');
end

if Options.elementwise,
    % we compare P and Q elementwise, i.e. P(1)<=Q(1), P(2)<=Q(1), ..., P(n)<=Q(1);
    % P(1)<=Q(2), P(2)<=Q(2), ..., P(n)<=Q(2); ... P(1)<=Q(m), ..., P(n)<=Q(m)
    if lenP>0,
        if lenQ==0,
            lenQ=1;
        end
        status=zeros(lenP,lenQ);
        for ii=1:lenP,
            status(ii,:)=(le(P.Array{ii},Q,Options));   % recursive call
        end
        return
    end
    if lenQ>0,
        status=zeros(1,lenQ);
        for ii=1:lenQ,
            status(ii)=any(le(P,Q.Array{ii},Options));  % recursive call
        end
        return
    end
else
    if lenP>0 | lenQ>0,
        % try to rule out some cases based on bounding boxes
        bboxOpt.noPolyOutput = 1;    % tell bounding_box() not to create a polytope object
        [R, Plow, Pup] = bounding_box(P, bboxOpt);
        [R, Qlow, Qup] = bounding_box(Q, bboxOpt);
        if ~isempty(Plow) & ~isempty(Qlow),
            bbox_tol = 1000*Options.abs_tol;
            if any(Plow + bbox_tol < Qlow) | any(Pup - bbox_tol > Qup),
                % bounding box of P violates bounding box of Q, hence P cannot be a
                % subset of Q
                status = 0;
                return
            end
            % we cannot reach any conclusion based solely on the fact that bounding
            % boxes are identical, therefore we continue...
        end

        Options.simplecheck=1;
        status = (~isfulldim(mldivide(P,Q,Options)));   % P is subset of Q if P\Q is empty polytope
        return
    end
end

[ncP,nxP]=size(P.H);
[ncQ,nxQ]=size(Q.H);

abs_tol = Options.abs_tol;

Pbbox = P.bbox;
Qbbox = Q.bbox;
if ~isempty(Pbbox) & ~isempty(Qbbox),
    bbox_tol = 1000*abs_tol;
    if any(Pbbox(:,1) + bbox_tol < Qbbox(:,1)) | any(Pbbox(:,2) - bbox_tol > Qbbox(:,2)),
        % bounding box of P violates bounding box of Q, hence P cannot be a
        % subset of Q
        status = 0;
        return
    end
end

status=1;
for ii=1:ncQ
    [xc,Rc]=chebyball_f([P.H;-Q.H(ii,:)],[P.K;-Q.K(ii)],Options);
    if (Rc > abs_tol)   % if ii-th border of Q is crossed with P then P~<=Q
        status=0;
        return;
    end
end
