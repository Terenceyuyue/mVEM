function P = subsasgn(P, X, Q)
%SUBSASGN Indexed assignments for polytope objects
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% 
% P(I) = Q B assigns the values of B into the elements of A specifed by
%            the subscript vector I.
%
% Note: if P(I) = [], removes polytope at position 'I' from the polytope array
%
% USAGE:
%   P(2) = [];
%   P(4) = Q1;
%   P([1 3 5]) = Q2;
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% see also SUBSREF, HORZCAT
%

% Copyright is with the following author(s):
%
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


if isempty(P),
    P=polytope;
end

oldlen=length(P.Array)+1;

if length(Q)>1,
    error('SUBSASGN: Polyarrays not allowed as right-hand side arguments!');
end

Xs = X.subs{1};
if isempty(Xs),
    % corresponds to P([]) = something
    return
end

%%% if length(Xs)==1 & Xs(1)==1 & length(P.Array)==0,
if length(Xs)==1 & Xs(1)==1 & isempty(P.Array),
    %%% if length(P.Array)==0,
    if isempty(P.Array),
        if isempty(Q),
            P=polytope;
        else
            P=Q;
        end
    else
        P.Array{1}=Q;
    end
    return
end

if isempty(P.Array);
    P1 = P;
    P1.Array = {};
    P.Array{1} = P1;
end
if isempty(Q),
    R = polytope;
    for jj=1:length(P.Array),
        if ~any(jj==Xs),
            R = [R P.Array{jj}];
        end
    end
    P=R;
else
    if ~isa(Q,'polytope'),
        error('SUBSASGN: Argument MUST be a polytope object!');
    end
    if dimension(Q)~=dimension(P),
        error('SUBSASGN: Dimensions of all polytopes must be identical.');
    end
    for ii=1:length(Xs),
        P.Array{Xs(ii)}=Q;
    end
    for ii=oldlen:length(P.Array)
        if isempty(P.Array{ii}),
            Q.H = 1;
            Q.K = -Inf;
            Q.normal = logical(0);
            Q.minrep = logical(0);
            Q.xCheb = [];
            Q.RCheb = -Inf;
            %%Q.keptrows = 0;
            Q.Array = {};            
            Q.bbox = [];
            P.Array{ii}=Q;
        end
    end
end
