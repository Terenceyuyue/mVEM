function [value, regions] = mpt_evalPWA(Pn, Fi, Gi, x0)
%MPT_EVALPWA Evaluates a PWA function at a given location
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Evaluates a PWA function on a given location
%
% USAGE:
%   value = mpt_evalPWA(Pn, Fi, Gi, x0)
%   [value, regions] = mpt_evalPWA(Pn, Fi, Gi, x0)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Pn         - polyhedral partition over which the PWA function is defined
% Fi, Gi     - cell arrays containing parameters of the PWA function
% x0         - value of parameter on which to evaluate the function
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% value      - value of the PWA function at given location (if point x0 is
%              contained in more than one region of Pn, "value" will be returned
%              as a vector, each element corresponding to cost associated to
%              every region containing x0)
% regions    - region(s) of Pn which contain(s) point x0
%

% Copyright is with the following author(s):
%
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

error(nargchk(4,4,nargin));

if ~isa(Pn,'polytope')
    error('mpt_evalPWA: first input argument must be a polytope.');
end

if ~iscell(Fi) | ~iscell(Gi),
    error('mpt_evalPWA: second and third arguments must be cells.');
end

if length(Fi) ~= length(Gi),
    error('mpt_evalPWA: length of Fi must be equal to length of Gi.');
end

if length(Pn) ~= length(Fi),
    error('mpt_evalPWA: length of Pn must be equal to length of Fi.');
end

x0 = x0(:);
if length(x0) ~= dimension(Pn),
    error('mpt_evalPWA: incorrect dimension of x0.');
end

if size(Fi{1}, 2) ~= length(x0),
    error('mpt_evalPWA: dimensions of Fi{1} and x0 must match.');
end
if size(Fi{1}, 1) ~= size(Gi{1}, 1),
    error('mpt_evalPWA: dimensions of Fi{1} and Gi{1} must match.');
end

value = [];
[isin, regions] = isinside(Pn, x0);
if ~isin,
    return
end

for ii = 1:length(regions),
    region = regions(ii);
    value = [value; Fi{region}*x0 + Gi{region}];
end
