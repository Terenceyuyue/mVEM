function [R,fulldim] = intersect(P1,P2,Options)
%INTERSECT Intersection of 2 polytopes or polytope arrays
%
% R = intersect(P1,P2,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns normalized minimal representation of intersection of input
% arguments
%
% USAGE:
%   R=intersect(P1,P2)
%   R=intersect(P1,P2,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P1,P2   - Polytopes or polyarrays, where a polyarray is considered as 
%           the union of polytopes
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R - Polytope containing the intersection of input polytopes
%
%
% Note: INTERSECT and AND have different functionality.
%
% see also AND
%

% Copyright is with the following author(s):
%
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
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

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin==1 & isa(P1, 'polytope')
    if isempty(P1.Array),
        P1.Array{1} = P1;
    end
    lenP1 = length(P1.Array);
    R = P1.Array{1};
    for i = 2:lenP1
        R = R & P1.Array{i};
    end
    fulldim = isfulldim(R);
    return
end

if nargin<3
    Options = [];
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end
if ~isfield(Options,'verbose')
    Options.verbose=mptOptions.verbose;    % level of verbosity
end
if ~isfield(Options,'reduce_intersection') % if set to 0, the intersection will not be reduced, i.e. no removal of redundant constraints takes place 
    Options.reduce_intersection=1;
end
if ~isfield(Options, 'lpsolver'),
    Options.lpsolver = mptOptions.lpsolver;
end

normal = 0;
if Options.reduce_intersection==0,
    reduceit=2;
    Options.simplecheck=1;         % for fast regiondiff call
    normal=1;
else
    reduceit=0;
end
    
maxdimP=0;


normal=1;


if ~isa(P1,'polytope') | ~isa(P2,'polytope')
    error('INTERSECT: arguments MUST be a polytope object!');
end

lenP1 = length(P1.Array);
lenP2 = length(P2.Array);
fulldim = 0;

if dimension(P1)~=dimension(P2),
    error('INTERSECT: Polytopes must be of same dimension!');
end

if lenP1>0 | lenP2>0,
    R = mldivide(P1,mldivide(P1,P2,Options),Options);
    fulldim = isfulldim(R(1));
    return
else
    havebboxes = 0;
    if ~isempty(P1.bbox) & ~isempty(P2.bbox)
        bbox_tol = Options.abs_tol*1e4;
        havebboxes = 1;
        if any(P1.bbox(:,2) + bbox_tol < P2.bbox(:,1)) | any(P2.bbox(:,2) + bbox_tol < P1.bbox(:,1))
            % bounding boxes do not intersect => polytopes do not intersect
            fulldim = 0;
            R=mptOptions.emptypoly;
            return
        end
    end

    if ~isfulldim(P1),
        % intersection with an empty polytope is an empty polytope
        fulldim = 0;
        R = mptOptions.emptypoly;
        return    
    else
        H1 = P1.H;
        K1 = P1.K;
    end
    if ~isfulldim(P2),
        % intersection with an empty polytope is an empty polytope
        fulldim = 0;
        R = mptOptions.emptypoly;
        return    
    else
        H2 = P2.H;
        K2 = P2.K;
    end
    
    HH = [H1; H2];
    KK = [K1; K2];
    
    if havebboxes,
        l = max([P1.bbox(:,1) P2.bbox(:,1)]')';
        u = min([P1.bbox(:,2) P2.bbox(:,2)]')';      
        cand = find(~((HH>0).*HH*(u-l) - (KK-HH*l) < -bbox_tol));
        if ~all(cand)
            cand = find(cand)
            HH = HH(cand,:);
            KK = KK(cand,:);
        end
    end

    [xcheb, rcheb] = chebyball_f(HH, KK, Options);
    fulldim = (rcheb > Options.abs_tol);
    if fulldim
        % intersection is obtained by eliminating redundant constraints from the system of inequalities
        
        % note that we also provide xcheb and rcheb to polytope(), otherwise we
        % would re-compute these two parameters in the polytope constructor
        R = polytope(HH, KK, normal, reduceit, xcheb, rcheb);
    else
        R = mptOptions.emptypoly;
    end
    
end
