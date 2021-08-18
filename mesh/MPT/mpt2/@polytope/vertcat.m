function P = vertcat(varargin)
%VERTCAT Concatenates polytopes into a polytope array
%
% P = vertcat(varargin),
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Concatenates input arguments and creates a polytope array
%
% USAGE:
%   P=[P1; P2; P3]    
%   P=[[P1; P2]; P2];
%
% Note that only 1-dimensional arrays are supported for polytopes!
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% Polytopes or polyarrays
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% P   - concatenated array of polytopes
%
% see also HORZCAT, SUBSREF, SUBSASGN
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


nargs=nargin;

for ii=1:nargs
    if isa(varargin{ii},'double') & isempty(varargin{ii}),
        varargin{ii} = polytope;
    end
    if ~isa(varargin{ii},'polytope')
        error('HORZCAT: input argument MUST be a polytope!');
    end
end

P=varargin{1};
    
if nargs==1,
    return
end

if length(P.Array)>0,
    dimP=size(P.Array{1}.H,2);
elseif isfulldim(P)
    dimP=size(P.H,2);    % dimension of the first polytope
else
    for ii=1:nargs
        Q=varargin{ii};
        if length(Q.Array)>0,
            dimP=size(Q.Array{1}.H,2);
            break
        else
            dimP=size(Q.H,2);
        end
    end
end

if isempty(P.Array),
    if P.K(1)~=-Inf,
        Q=P;                    % all remaining arguments will be stored in the Array field of the first input argument
        P.H = 1;
        P.K = -Inf;
        P.normal = logical(1);
        P.minrep = logical(1);
        P.xCheb = 1;
        P.RCheb = -Inf;
        %%P.keptrows = 0;
        P.Array{1}=Q;           % store the first polytope into P.Array{1}
        P.vertices=[];
        P.bbox = [];
    end
end

for jj=2:nargs,
    Q = varargin{jj};
    if ~isempty(Q.Array),         % in case Q is a polyarray
        for ii=1:length(Q.Array)  % cycle through all it's elements
            if Q.Array{ii}.H==1 & Q.Array{ii}.K==-Inf & Q.Array{ii}.RCheb==-Inf,
                % if polytope is empty, skip it
                continue
            elseif size(Q.Array{ii}.H,2)~=dimP,
                % we merge only polytopes of the same dimension
                error('When concatenating polytopes, all input arguments must have same dimension!');
            else
                % otherwise add it to the structure
                P.Array{length(P.Array)+1} = Q.Array{ii};
            end
        end
    else
        if Q.H==1 & Q.K==-Inf & Q.RCheb==-Inf,
            continue
        elseif size(Q.H,2)~=dimP,
            error('When concatenating polytopes, all input arguments must have same dimension!');
        else
            P.Array{length(P.Array)+1} = Q;
        end
    end
end

if length(P.Array)==1,
    P = P.Array{1};
    P.Array = {};
end
