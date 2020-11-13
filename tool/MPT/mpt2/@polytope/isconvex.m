function [status, Pconv] = isconvex(P,Options)
%ISCONVEX Checks if a polytope array forms a convex union
%
% [status, Pconv] = isconvex(P)
% [status, Pconv] = isconvex(P,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% STATUS = ISCONVEX(P) returns TRUE (1) if a given polytope array P is convex
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                - polytope array
% Options.usehull  - if set to 1 (default is 0), uses convex hulls instead of
%                    convex envelopes to detect convexity
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% status           - Logical statement
% Pconv            - if P is convex it returns the convex set, otherwise 
%                    an empty polytope
%

% Copyright is with the following author(s):
%
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch

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
    Options = [];
end

if ~isfield(Options, 'usehull'),
    Options.usehull = 0;
end
if ~isfield(Options, 'abs_tol'),
    Options.abs_tol = mptOptions.abs_tol;
end
if ~isfield(Options, 'bbox_tol'),
    % higher tolerance to check bounding boxes
    % it is so high to avoid numerical problems
    Options.bbox_tol = 1e4*mptOptions.abs_tol;
end

if ~isa(P, 'polytope'),
    error('ISCONVEX: First input must be a polytope object.');
end

if ~isfulldim(P),
    % non-fully dimensional polytopes are assumed to be convex
    status = 1;
    Pconv  = P;
    return
end
if isempty(P.Array),
    % single polytope, it is for sure convex
    status = 1;
    Pconv  = P;
    return
end

if Options.usehull==1,
    outer = hull(P, Options);
else
    outer = envelope(P, Options);
end

bboxOpt.noPolyOutput = 1;  % tell bounding_box that we just need vertices
[R, Pl, Pu] = bounding_box(P, bboxOpt);
[R, Ol, Ou] = bounding_box(outer, bboxOpt);
bboxP = [Pl Pu];
bboxO = [Ol Ou];

if any(isnan(Ol)) | any(isnan(Ou)) | isempty(Ol) | isempty(Ou)
    % if envelope returns R^n, bounding_box returns NaNs
    status = 0;
    Pconv = mptOptions.emptypoly;
    return
end    

bbox_tol = Options.bbox_tol;
if any(abs(bboxP(:,1) - bboxO(:,1)) > bbox_tol) | any(abs(bboxP(:,2) - bboxO(:,2)) > bbox_tol),
    % bounding boxes differ by more than bbox_tol => polytopes cannot be equal
    % therefore P is not convex
    status = 0;
    Pconv = mptOptions.emptypoly;
    return
end
% we cannot reach any conclusion based solely on the fact that bounding
% boxes are identical, therefore we continue...

mldivideOpt = Options;
mldivideOpt.simplecheck = 1;   % tell set difference not to construct polytopes
if isfulldim(mldivide(outer, P, mldivideOpt)),
    % if set difference between the envelope (hull) and P is fully dimensional,
    % it means that P is not convex
    status = 0;
    Pconv  = mptOptions.emptypoly;
else
    % set difference is empty => P is convex
    status = 1;
    Pconv  = outer;
end
