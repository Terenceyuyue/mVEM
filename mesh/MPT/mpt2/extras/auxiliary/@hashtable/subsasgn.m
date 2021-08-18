function obj = subsasgn(obj, X, V)
%SUBSREF Subreferencing operator for HASHTABLE objects
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
%

% Copyright is with the following author(s):
%
%(C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

indfound = [];
for ii = 1:length(X.subs),
    key = X.subs{ii};
    if ~ischar(key),
        error('Key must be a char.');
    end
    indfound = [indfound find(ismember(obj.keys, key))];
end
if isempty(indfound),
    for ii = 1:length(X.subs),
        obj.keys{end+1} = X.subs{ii};
        obj.values{end+1} = V;
    end
else
    obj.values(indfound) = V;
    if length(indfound)==1,
        obj.values{indfound} = V;
    else
        for ii = 1:length(indfound),
            obj.values{indfound(ii)} = V;
        end
    end
end
