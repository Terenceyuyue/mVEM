function Q = subsref(P, X)
%SUBSREF Indexed referencing for polytope objects
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% P(I) is an array formed from the elements of A specifed by the
%      subscript vector I.
%
% If P is a polytope, P(1) returns that polytope
%
% If P is an array of polytopes (e.g. P=[P1 P2 P3]), P(2) returns P2, etc.
%
% Indexing by logicals is supported as well, e.g.
%
%    P(logical([0 1 0 1]))
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% see also SUBSASGN, END
%

% Copyright is with the following author(s):
%
% (C) 2003-2006 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch
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


if numel(X)>1,
    error('??? Attempt to reference field of non-structure array.');
else
    if (~strcmp(X.type,'()')),
        % only indexes in round brackets are allowed
        if X.type(1)=='.',
            error(['Indexing with ''' X.type ''' not supported! Use [H,K]=double(P) to access the H-representation; type ''help polytope'' for more details']);
        end
        error(['Indexing with ''' X.type ''' not supported!']);
    end
end

indices = X.subs{1};
if islogical(indices)
    indices = double(find(indices));
elseif any(indices<0),
    error('??? Subscript indices must either be real positive integers or logicals.');
end
if isempty(indices),
    Q = polytope;
    return
end

lenP = length(P.Array);
lenInd = length(indices);

if (lenP==0)
    if any(indices>1),
        error('??? Index exceeds array dimension');
    end
    Q = P;
    if lenInd > 1,
        Q.Array = cell(1, lenInd);
        for i = 1:lenInd,
            Q.Array{i} = P;
        end
    end
    
else
    if any(indices>lenP),
        error('??? Index exceeds matrix dimension');
    end
    if (lenInd==1),
        Q = P.Array{indices};
    else
        Q = P;
        Q.Array = P.Array(indices);
    end
end
