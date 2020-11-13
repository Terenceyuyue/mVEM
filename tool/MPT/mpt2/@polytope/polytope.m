function [P] = polytope(varargin)
%POLYTOPE Default constructor for the POLYTOPE object
%
% POLYTOPE  Defines a new POLYTOPE object
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%
% NOTE: preferred call is P=polytope(H,K)
%
%   P = POLYTOPE     Creates an empty polytope in R^1
%
%   P = POLYTOPE(Q)  Creates polytope P={x | H x <= K} from a structure Q
%                    Structure Q has following fields:
%                    H       - nc x nx matrix, (Note: Required field)
%                    K       - nc x 1  vector, (Note: Required field)
%                    normal  - is polytope in a normalized representation,
%                              {0=no, 1=yes}, default: 0,
%                    minrep  - is polytope in a minimal representation,
%                              {0=no, 1=yes}, default: 0,
%                    xCheb   - Center of a Chebyshev ball (i.e., interior point),
%                    RCheb   - Radius of a Chebyshev ball
%
%                    Note: If Q is of POLYTOPE class then P=Q (simple copy).
%
%   P = POLYTOPE(V)  Creates polytope by making a convex hull of the given vertices V
%
% Following calls may be used as well
%   P = POLYTOPE(H,K)
%   P = POLYTOPE(H,K,normal)
%   P = POLYTOPE(H,K,normal,minrep)
%   P = POLYTOPE(H,K,normal,minrep,xCheb,RCheb)
%
% ---------------------------------------------------------------------------
% ACCESSING INTERNAL DATA OF THE POLYTOPE STRUCTURE
% ---------------------------------------------------------------------------
%
% Each polytope object is internally represented by a structure with the above
% mentioned fields. It is not possible, however, to access these fields directly
% from the Matlab environment using the dot (.) delimiter (e.g. P.H). To
% retrieve these internal data, one has to use one of the following functions:
%
% To access the H-representation (P.H*x <= K):
%   [H,K] = double(P)
%
% To access center and radius of the Chebyshev's ball:
%   [xCheb,RCheb] = chebyball(P)
%
% To get information if a polytope is in normalized representation
%   status = isnormal(P)
%
% To get information if a polytope is in minimal representation
%   status = isminrep(P)
%
% To get information if a polytope is fully dimensional
%   status = isfulldim(P)
%
% To get information if a polytope is bounded
%   status = isbounded(P)
%
% To retrieve extreme points of a polytope (i.e. the V-representation)
%   V = extreme(P)
%
% 
% Consult individual help files and/or the MPT manual for more details
%
% ---------------------------------------------------------------------------
% MERGING OF POLYTOPES INTO POLYARRAYS
% ---------------------------------------------------------------------------
%
% Polytopes can be concatenated into arrays of polytopes (polyarrays). Each
% polyarray is again a polytope object, which allows to use them in any
% overloaded function which accepts individual polytope objects (check the
% manual for exceptions).
%
% Polytopes are concatenated into arrays using the [,] operator as follows:
%
%   Q = [P1, P2, P3]
%
% The above command creates a polyarray with 3 elements (P1, P2 and P3). Each of
% these elements can be a polyarray, e.g.
%
%   R = [P4, Q, [P5, P6], P7]
%
% A polyarray can be indexed using the standard (i) operator, i.e.
%
%   Q(2) will return the second element of polyarray Q
% 
% Check help polytope/subsref and polytope/subsasgn for more details.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Q           - can be either a structure (described uppon) or a matrix with vertices
% H,K         - H-representation of the polytope
% normal      - set to 1 if you know the pair (H,K) is already normalized, 0 otherwise
% minrep      - set to 1 if you knot the pair (H,K) is already reduced, 0 otherwise
% xCheb,RCheb - user provided Chebyshev's ball parameters
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% P   - polytope
%
%
% see also POLYTOPE/DOUBLE, CHEBYBALL, ISFULLDIM, ISMINREP, ISNORMAL, EXTREME
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


global mptOptions;
if ~isstruct(mptOptions),
    mpt_error;
end

% fast call
%----------
if nargin==1 & isa(varargin{1},'polytope')
    P = varargin{1};
    return;
end
    
Options=[];

% default values
%---------------
P.H=1;
P.K=-Inf;
P.normal  =  logical(0);     % is polytope in a normalized representation, {0=no, 1=yes}
P.minrep  =  logical(0);     % is polytope in a minimum representation, {0=no, 1=yes}

P.xCheb = [];
P.RCheb = [];
P.Array = {};
P.vertices = [];
P.bbox = [];
normalit=1;
reduceit=1;


if nargin==1 & iscell(varargin{1}),
    % converts a cell array of polytopes to a polytope array
    Q = varargin{1};
    for ii=1:length(Q),
        if ~isa(Q{ii},'polytope'),
            error('POLYTOPE: all elements in the cell array must be polytope objects!');
        end
        if length(Q{ii}.Array)>0,
            for jj=1:length(Q{ii}.Array),
                if isfulldim(Q{ii}.Array{jj}),
                    P.Array{end+1} = Q{ii}.Array{jj};
                end
            end
        else
            if length(Q)==1,
                P = struct(Q{ii});
            elseif isfulldim(Q{ii}),
                P.Array{end+1} = Q{ii};
            end
        end
    end
    P = class(P,'polytope');
    return
end


if nargin==0
    % empty 1-dimensional polytope
    P.H=1;
    P.K=-Inf;
    P.normal=logical(1);
    P.minrep=logical(1);
    P.xCheb=0;
    P.RCheb=-Inf;
    P.Array={};
    P.vertices=[];
    P.bbox=[];
    reduceit=0;
    normalit=0;
elseif nargin==1 
    if isa(varargin{1},'struct')
        % input is a structure; create POLYTOPE object from it
        if isfield(varargin{1},'H')
            P.H=varargin{1}.H;
        else
            error('POLYTOPE: Input structure MUST have ''H'' field');
        end
        if isfield(varargin{1},'K')
            P.K=varargin{1}.K;
        else
            error('POLYTOPE: Input structure MUST have ''K'' field');
        end
        if isfield(varargin{1},'normal')
            P.normal=logical(varargin{1}.normal);
        end
        if isfield(varargin{1},'minrep')
            P.minrep=varargin{1}.minrep;
        end
        if isfield(varargin{1},'xCheb')
            P.xCheb=varargin{1}.xCheb;
        end
        if isfield(varargin{1},'RCheb')
            P.RCheb=varargin{1}.RCheb;
        end
        %if isfield(varargin{1},'keptrows')
        %    P.keptorws=varargin{1}.keptrows;
        %end
        if isfield(varargin{1},'Array')
            P.Array=varargin{1}.Array;
        end
        if isfield(varargin{1},'vertices')
            P.vertices=varargin{1}.vertices;
        end
        if isfield(varargin{1},'bbox')
            P.bbox=varargin{1}.bbox;
        end
    elseif isa(varargin{1},'double'),
        [P,P.vertices]=hull(varargin{1}); % compute the convex hull and remove points which are not extreme
        return
    else
        error('POLYTOPE: Input should be a POLYTOPE, or a STRUCT or a matrix with vertices');
    end
else
    P.H=varargin{1};
    P.K=varargin{2};
    if nargin>=3
        normalarg=varargin{3};
        normalit=(normalarg==0);
        P.normal=(normalarg==1);
        if nargin>=4
            reducearg=varargin{4};
            reduceit=(reducearg==0);
            P.minrep=(reducearg==1);
            if nargin>=5
                P.xCheb=varargin{5};
                if nargin>=6,
                    P.RCheb=varargin{6};
                    if nargin<6
                        error('POLYTOPE: Wrong number of input arguments');
                    elseif nargin==7,
                        Options=varargin{7};
                    end
                end
            end
        end
    end
end

if issparse(P.H),
    P.H = full(P.H);
end
if issparse(P.K),
    P.K = full(P.K);
end

[nc,nx]=size(P.H);

% check if arguments are correct
%-------------------------------
if nx<1
    error('POLYTOPE: Dimension of a polytope has to be nx >= 1');
end
if size(P.K,1)~=nc | size(P.K,2)~=1
    error(sprintf('POLYTOPE: K should be a row vector of length %d',nc));
end
if all(P.normal~=[0 1 2])
    error('POLYTOPE: Wrong value for NORMAL');
end
if all(P.minrep~=[0 1 2])
    error('POLYTOPE: Wrong value for MINREP');
end

P = class(P,'polytope');

if normalit,
    P=normalize(P);
end

if isempty(P.xCheb) | isempty(P.RCheb)
   [P.xCheb,P.RCheb,lambda]=chebyball(P);
end

if(P.RCheb<mptOptions.abs_tol)
    %if chebyball is nonexistent => return emptypoly
    P.H=1;
    P.K=-Inf;
    P.normal=logical(0);
    P.minrep=logical(1);
    P.xCheb=0;
    P.RCheb=-Inf;
    %%P.keptrows = [];
    P.Array={};
    P.vertices=[];
    P.bbox=[];
    return
end

if reduceit,
    P = reduce(P,Options);
end

return;
