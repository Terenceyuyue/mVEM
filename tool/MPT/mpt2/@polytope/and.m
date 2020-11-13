function R = and(varargin)
%AND Intersection of n polytopes
%
% R = and(varargin)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Intersection of n polytopes
%
%   R = AND(P1,P2,P3,...) returns normalized minimal representation
%   of intersection of n polytopes P1, P2, P3,..., Pn
%
% USAGE:
%   R = and(P1,P2,P3)
%   R = P1 & P2
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Polytopes delimited by comma. If P1 (or P2, etc.) are polytope arrays,
% then the elements of P1 are considered as seperate polytopes, i.e
% AND(P1) is the same as AND(P1(1),P1(2),P1(3),...).
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R - Polytope containing the intersection of input arguments
%
%
% Note: INTERSECT and AND have different functionality.
%
% see also INTERSECT
%

% Copyright is with the following author(s):
%
% (C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%          fjc@control.ee.ethz.ch
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

global mptOptions;
Options=[];
if ~isstruct(mptOptions),
    mpt_error;
end

if ~isfield(Options,'rel_tol')
    Options.rel_tol=mptOptions.rel_tol;    % relative tolerance
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end
if ~isfield(Options,'verbose')
    Options.verbose=mptOptions.verbose;
end

ni=nargin;
HH=[];
KK=[];
normal=1;
maxdimP=0;
for ii=1:ni,
    % determine maximum dimension
    if ~isa(varargin{ii}, 'polytope'),
        error('AND: all input arguments must be polytope objects.');
    end
    dimP=dimension(varargin{ii});
    if dimP>maxdimP,
        maxdimP=dimP;
    end
end
for ii=1:ni,   % go through all input arguments
    P = varargin{ii};
    lenP = length(P.Array);  
    if lenP>0, % if the input argument is a polyarray, we cycle through all it's elements
        for jj=1:lenP,
            Q=P.Array{jj};    % pick up a polytope from polyarray
            if ~isfulldim(Q),
                Q.H=[];
                Q.K=[];
            else
                if size(Q.H,2)<maxdimP,
                    error('INTERSECT: Polytopes must be of same dimension!');
                end
            end
            HH=[HH; Q.H];     % concatenate all H and K matrices together
            KK=[KK; Q.K];
            normal = normal & Q.normal;
            if Q.RCheb < Options.abs_tol,   % if the polytope is empty, we return
                if Options.verbose>=2,
                    disp('warning: empty polytope detected!');
                end
                R=polytope;
                return
            end
        end
    else
        if ~isfulldim(P),
            P.H=[];
            P.K=[];
        else
            if size(P.H,2)<maxdimP,
                error('INTERSECT: Polytopes must be of same dimension!');
            end
        end
        HH=[HH; P.H];
        KK=[KK; P.K];
        normal = normal & P.normal;
        if P.RCheb < Options.abs_tol,
            if Options.verbose>=2,
                disp('warning: empty polytope detected!');
            end
            R=mptOptions.emptypoly;
            return
        end
    end
end

Options.lpsolver = mptOptions.lpsolver;
[xcheb, rcheb] = chebyball_f(HH, KK, Options);
if rcheb > Options.abs_tol,
    % intersection is a full-dimensional polytope
    
    % note that we also provide xcheb and rcheb to polytope(), otherwise we
    % would re-compute these two parameters in the polytope constructor
    R = polytope(HH, KK, normal, 0, xcheb, rcheb);
else
    R = mptOptions.emptypoly;
end
