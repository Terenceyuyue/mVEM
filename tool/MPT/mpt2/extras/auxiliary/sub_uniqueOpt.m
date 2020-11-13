function table = sub_uniqueOpt(Fi, Gi, nu)
% SUB_UNIQUEOPT identifies unique optimizers
%
% table = sub_uniqueOpt(Fi, Gi, nu)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Identifies unique optimizers
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Fi, Gi   - cell arrays defining the optimizer as J = Fi{r}*x + Gi{r}
% nu       - number of system inputs (optional)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% table    - linking table

% Copyright is with the following author(s):
%
% (C) 2007 Michal Kvasnica, Slovak University of Technology in Bratislava
%          michal.kvasnica@stuba.sk

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

error(nargchk(2,3,nargin));

if ~iscell(Fi)
    Fi = { Fi };
    Gi = { Gi };
end
if length(Fi) ~= length(Gi)
    error('Length of Fi and Gi must be identical.');
end
if nargin < 3
    nu = size(Fi{1}, 1);
elseif nu > size(Fi{1}, 1) | nu < 1
    error('Index exceeds matrix dimension.');
end

nr = length(Fi);
table.Reg = NaN*ones(1,nr);
table.Table = {}; 
table.Fi = {};
table.Gi = {};
reg_set = 1:nr;
nx = size(Fi{1}, 2);

tolEq = 10e-10;
d=0;

% only keep first "nu" lines in Fi, Gi
if any(cellfun('prodofsize', Fi) > nu*nx) & ...
        nu < size(Fi{1}, 1)
    for i = 1:nr
        Fi{i} = Fi{i}(1:nu, :);
        Gi{i} = Gi{i}(1:nu);
    end
end

% create a hash of each Fi{r}, Gi{r} pair to speed up comparison
hashes = zeros(1, nr);
randv = randn(nx+1, 1);
for i = 1:nr
    hashes(i) = sum([Fi{i} Gi{i}]*randv);
end

% find identical elements
while length(reg_set) > 0
    
    % add the first region in the set to the new entry in the table 
    % and remove it from the set
    r = reg_set(1);
    d = d+1;
    table.Table{d} = r;
    table.Reg(r) = d;
    reg_set(1) = [];
    
    % check which Fi,Gi pairs are identical
    for k = reg_set
        % fast implementation
        if abs(hashes(r) - hashes(k)) <= tolEq 
            if all( abs(Gi{r}-Gi{k}) <= tolEq )
                if all(all( abs(Fi{r}-Fi{k}) <= tolEq ))
                    table.Table{d}(end+1) = k;
                    table.Reg(k) = d;
                    reg_set(find(reg_set==k)) = [];
                end
            end
        end
    end
    
end

% extract unique Fi, Gi elements
n_unique = length(table.Table);
idx = zeros(1, n_unique);
for i = 1:n_unique
    idx(i) = table.Table{i}(1);
end
table.Fi = Fi(idx);
table.Gi = Gi(idx);
