function P = horzcat(varargin)
%HORZCAT Concatenates polytopes into a polytope array
%
% P = horzcat(varargin),
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Concatenates input arguments and creates a polytope array
%
% USAGE:
%   P=[P1 P2 P3]    
%   P=[P1 P2]; Q=[P P3]
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
% see also VERTCAT, SUBSREF, SUBSASGN
%

% Copyright is with the following author(s):
%
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
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


nargs=nargin;

for ii=1:nargs
    argii = varargin{ii};
    if isa(argii,'double')
        if isempty(argii),
            emptypoly = polytope;
            varargin{ii} = emptypoly;
            argii = emptypoly;
        end
    end
    if ~isa(argii,'polytope')
        error('HORZCAT: input argument MUST be a polytope!');
    end
end

P=varargin{1};
    
if nargs==1,
    return
end

PArray = P.Array;

dimP = 0;
if ~isempty(PArray),
    dimP=size(PArray{1}.H,2);
elseif isfulldim(P)
    dimP=size(P.H,2);    % dimension of the first polytope
else
    for ii=1:nargs
        Q=varargin{ii};
        QArray = Q.Array;
        if ~isempty(QArray),
            dimP=size(QArray{1}.H,2);
            break
        elseif isfulldim(Q)
            dimP=size(Q.H,2);
        end
    end
end
if dimP == 0
    P = polytope;
    return
end

if isempty(PArray),
    if P.K(1)~=-Inf,
        Q=P;                    % all remaining arguments will be stored in the Array field of the first input argument
        P.H = 1;
        P.K = -Inf;
        P.normal = logical(1);
        P.minrep = logical(1);
        P.xCheb = 1;
        P.RCheb = -Inf;
        P.Array{1} = Q;           % store the first polytope into P.Array{1}
        P.vertices = [];
        P.bbox = [];
    end
end

for jj=2:nargs,
    Q = varargin{jj};
    QArray = Q.Array;
    if ~isempty(QArray),         % in case Q is a polyarray
        for ii=1:length(QArray)  % cycle through all it's elements
            QAHii = QArray{ii}.H;
            if QAHii==1 & QArray{ii}.K==-Inf & Q.Array{ii}.RCheb==-Inf,
                % if polytope is empty, skip it
                continue
            elseif size(QAHii,2)~=dimP,
                % we merge only polytopes of the same dimension
                error('When concatenating polytopes, all input arguments must have same dimension!');
            else
                % otherwise add it to the structure
                %P.Array{length(P.Array)+1} = QArray{ii};
                P.Array{end+1} = QArray{ii};
            end
        end
    else
        H = Q.H;
        if H==1 & Q.K==-Inf & Q.RCheb==-Inf,
            continue
        elseif size(H,2)~=dimP,
            error('When concatenating polytopes, all input arguments must have same dimension!');
        else
            P.Array{end+1} = Q;
        end
    end
end


% for ii=2:nargs,
%     Q = varargin{ii};
%     if ~isempty(Q.Array),         % in case Q is a polyarray
%         for ii=1:length(Q.Array)  % cycle through all it's elements
%             if Q.Array{ii}.H==1 & Q.Array{ii}.K==-Inf & Q.Array{ii}.RCheb==-Inf,
%                 % if polytope is empty, skip it
%                 continue
%             elseif size(Q.Array{ii}.H,2)~=dimP,
%                 % we merge only polytopes of the same dimension
%                 error('When concatenating polytopes, all input arguments must have same dimension!');
%             else
%                 % otherwise add it to the structure
%                 P.Array{length(P.Array)+1} = Q.Array{ii};
%             end
%         end
%     else
%         if Q.H==1 & Q.K==-Inf & Q.RCheb==-Inf,
%             continue
%         elseif size(Q.H,2)~=dimP,
%             error('When concatenating polytopes, all input arguments must have same dimension!');
%         else
%             %P.Array{length(P.Array)+1} = Q;
%             P.Array{end+1} = Q;
%         end
%     end
% end

if length(P.Array)==1,
    P = P.Array{1};
    P.Array = {};
end
