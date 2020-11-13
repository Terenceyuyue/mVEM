function X = grid(P, gridpoints)
%GRID Returns equidistantly spaced data points contained in a polytope/polyarray
%
% X = grid(P, gridpoints)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Returns a list of equidistantly spaced points in a given polytope or a
% polytope array.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P               - Polytope object
% gridpoints      - Number of grid points
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% X     - set of equidistantly spaced points
%

% Copyright is with the following author(s):
%
%(C) 2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(2, 2, nargin));

bbOptions.noPolyOutput = 1; % we don't need the bounding box as a polytope object
[B, lb, ub] = bounding_box(P, bbOptions);

% grid the state-space into equidistantly placed points
dimB = size(lb(:),1);
Xpoints = zeros(gridpoints, dimB);
for ii=1:dimB
    Xpoints(:,ii) = linspace(lb(ii),ub(ii),gridpoints)';
end
    
% generate all possible combinations of states
% one could use kron() here, but that one fails for high number of elements
n_states = dimB;
ZZ=[];
ZZ{n_states}=Xpoints(:,n_states);
for ii=n_states-1:-1:1,
    Zd=[];
    for jj=1:size(Xpoints,1),
        Zd=[Zd; repmat(Xpoints(jj,ii),length(ZZ{ii+1}),1)];
    end
    ZZ{ii}=Zd;
end
for ii=2:n_states,
    ZZ{ii}=repmat(ZZ{ii},length(ZZ{ii-1})/length(ZZ{ii}),1);
end
datapoints=[];
for ii=1:n_states,
    datapoints(:,ii) = ZZ{ii};
end
npoints = size(datapoints,1);
isin = zeros(npoints, 1);
locOpt.fastbreak = 1;
for i = 1:npoints,
    isin(i) = isinside(P, datapoints(i, :)', locOpt);
end
X = datapoints(find(isin), :);
