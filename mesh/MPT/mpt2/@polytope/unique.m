function [P, ind] = unique(P)
%UNIQUE Removes redundant entries from a polytope array
%
% [P, ind] = unique(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% unique(P) removes redundant polytopes from a polytope array P.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P   - polytope array
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% P   - polytope array containing only unique elements
% ind - indicies of unique regions
%
% see also POLYTOPE/REDUCEUNION

% Copyright is with the following author(s):
%
% (C) 2005 Johan Loefberg, Automatic Control Laboratory, ETH Zurich,
%          loefberg@control.ee.ethz.ch

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

if isempty(P.Array),
    % input is a single polytope, nothing to do
    ind = 1;
    return
end

hash = randn(size(P.Array{1}.H,2)+1,1);
HASH = 0;
for i = 1:1:length(P.Array)
    HASH(1:size(P.Array{i}.H,1),i) = sort([P.Array{i}.H P.Array{i}.K]*hash);
end

HASH = HASH'*randn(size(HASH,1),1);

[ii, jj] = unique(HASH,'rows');
ind = sort(jj);
P.Array = {P.Array{ind}};
