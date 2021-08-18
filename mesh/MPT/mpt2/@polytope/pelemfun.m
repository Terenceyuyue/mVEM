function varargout = pelemfun(funname, Pn, varargin)
%PELEMFUN Execute the specified function on each element of a polytope array
%
%   PELEMFUN(@FHANDLE, Pn, [X1, ..., Xn]) evaluates a function defined by the
%   function handle FHANDLE on all _elements_ of the polytope array Pn.
%   Additional input arguments X1, ..., Xn can be specified. Couple of examples:
%
%   E = pelemfun(@extreme, Pn) returns extreme points of all polytopes stored in
%   the polytope array Pn as a cell array.
%
%   [X, R] = pelemfun(@chebyball, Pn) returns centers and radii of the chebychev
%   ball of each element of Pn.
%
%   B = pelemfun(@bounding_box, Pn, struct('Voutput', 1)); V = [B{:}] returns
%   bounding boxes of all elements of Pn converted to a matrix.
%
%   I = pelemfun(@le, Pn, Q) returns a cell array of logical statements with
%   true value at position "i" if Pn(i) is a subset of Q.

% ---------------------------------------------------------------------------
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
% ---------------------------------------------------------------------------

error(nargchk(2,Inf,nargin));

if ~isa(funname, 'function_handle'),
    error('First input must be a function handle.');
end
if ~isa(Pn, 'polytope'),
    error('First input must be a polytope.');
end

nout = nargout;
if nout < 1,
    % always assume at least one output argument
    nout = 1;
end
if isempty(Pn.Array),
    % input is a single polytope
    o = cell(1, nout);
    [o{:}] = feval(funname, Pn, varargin{:});
    for io = 1:nout,
        % convert outputs to cells
        varargout{io}{1} = o{io};
    end
else
    % input is a polyarray
    lenP = length(Pn.Array);
    varargout = cell(1, nout);
    for io = 1:nout,
        varargout{io} = cell(1, lenP);
    end
    for ip = 1:lenP,
        o = cell(1, nout);
        [o{:}] = feval(funname, Pn.Array{ip}, varargin{:});
        for io = 1:nout,
            varargout{io}{ip} = o{io};
        end
    end
end
