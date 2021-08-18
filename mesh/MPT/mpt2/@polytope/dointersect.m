function [answer, IA, IB, IAB] = dointersect(P1,P2)
%DOINTERSECT Checks if two polytopes / polyarrays intersect
%
% [answer, IA, IB, IAB] = dointersect(P1,P2)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns true (1) if input polytopes / polyarrays intersect.
%
% This function does not return the intersection!
%
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P1,P2   - Polytopes or polyarrays
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% answer  - 1 if polytopes / polyarrays intersect, 0 otherwise
% IA, IB  - index vectors such that P1(IA) & P2(IB) = intersect(P1, P2)
% IAB     - indeces of intersecting polytopes stored as a n x 2 matrix, where
%           every row gives indeces of intersecting polytopes. returned as an
%           empty matrix if no intersections exist. Example:
%             IAB = [1 1; 2 1; 3 1] means that there exists a full-dimensional
%             intersection of (P1(1) and P2(1)), (P1(2) and (P2(1)) and (P1(3)
%             and (P2(1)).
%
% see also AND, INTERSECT
%

% Copyright is with the following author(s):
%
% (C) 2004-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch

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

bbox_tol = 1e4*mptOptions.abs_tol;

if ~isa(P1, 'polytope') | ~isa(P2, 'polytope'),
    error('DOINTERSECT: both input arguments must be polytope objects!');
end

if dimension(P1) ~= dimension(P2)
    error('DOINTERSECT: polytopes must have the same dimension!');
end

IAB = []; IA = []; IB = [];

if isempty(P1.Array) & isempty(P2.Array),
    % special case, both inputs are single polytopes

    if ~isempty(P1.bbox) & ~isempty(P2.bbox)
        if any(P1.bbox(:,2) + bbox_tol < P2.bbox(:,1)) | ...
                any(P1.bbox(:,1) - bbox_tol > P2.bbox(:,2))
            % even bounding boxes of the two polytopes do not intersect, abort
            % quickly
            answer = 0;
            return
        end
    end

    Hint = [P1.H; P2.H];
    Kint = [P1.K; P2.K];
    [xcheb, rcheb] = chebyball_f(Hint, Kint, mptOptions);
    answer =  (rcheb > mptOptions.abs_tol);
    if answer,
        IAB = [1 1];
        IA = 1;
        IB = 1;
    end
    return
    
else
    [H1, K1] = double(P1);
    [H2, K2] = double(P2);
    if ~iscell(H1),
        H1 = {H1};
        K1 = {K1};
    end
    if ~iscell(H2),
        H2 = {H2};
        K2 = {K2};
    end
    lenP1 = length(H1);
    lenP2 = length(H2);
    
    answer = 0;    
    bboxOpt.noPolyOutput = 1;  % do not return the bounding box as a polytope object
    [d, b1min, b1max] = bounding_box(P1, bboxOpt); 
    [d, b2min, b2max] = bounding_box(P2, bboxOpt);
    if any(b1max+bbox_tol < b2min) | any(b2max+bbox_tol < b1min),
        % even bounding boxes of the two polytopes do not intersect, abort
        % quickly
        return
    end

    for i1 = 1:lenP1,
        for i2 = 1:lenP2,
            Hint = [H1{i1}; H2{i2}];
            Kint = [K1{i1}; K2{i2}];
            [xcheb, rcheb] = chebyball_f(Hint, Kint, mptOptions);
            if (rcheb > mptOptions.abs_tol),
                answer = 1;
                if nargout > 1,
                    IAB = [IAB; i1 i2];
                else
                    % exit quickly if we don't care about indeces of
                    % intersecting regions
                    return
                end
            end
        end
    end
end

if nargout > 1 & ~isempty(IAB),
    IA = unique(IAB(:,1));
    IB = unique(IAB(:,2));
end
