function [hhPoly,adj,vvPoly,vol,P]=triangulate(P,Options)
%TRIANGULATE Calculates triangulation of arbitrary polytopes
%
% function [hPoly,adj,vPoly,vol,P]=triangulate(P,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Partitions any polytope (given in H or V representation) into simplical cells.
% The function returns the simplical cells in V (and H if flag set) representation. 
% An adjecency list is returned as well.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                       - Polytope object
% Options.noHPoly         - Set to "1" if you do not wish to get the H-representation of the 
%                           output simplices.
% Options.extreme_solver  - Which method to use for vertex enumeration
% Options.abs_tol         - Absolute tolerance
% Options.rel_tol         - Relative tolerance
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%
% vPoly                   - Cell array, where vPoly{i} is a matrix containing the 
%                           vertices of simplex i
% adj                     - adjecency list matrix, indexing the vertices of P which form 
%                           a simplex.
% hPoly                   - polytope array containing all simplices obtained for P
%
% vol                     - volume of the polytope P 
%
% P                       - Return polytope P, which now contains the vertices of P
%
%
% (C) 2004 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch

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

if nargin<2,
    Options=[];
end

if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end
if ~isfield(Options,'rel_tol')
    Options.rel_tol=mptOptions.rel_tol;    % relative tolerance
end
if ~isfield(Options,'extreme_solver')
    Options.extreme_solver=mptOptions.extreme_solver;
end
if ~isfield(Options,'debug_level')
    Options.debug_level=mptOptions.debug_level;
end
if ~isfield(Options,'noHPoly') 
    Options.noHPoly=0;  %compute H triangulation
end
if length(P.Array)>0,
    error('triangulate: This function does not support polytope arrays!');
end

% throw an error if the input polytope is unbounded
if ~isbounded(P),
    error('triangulate: Input polytope is unbounded');
end
    
minRad=Options.rel_tol;         %simplices with a chebychev radius smaller than this are considered empty
minVol=Options.rel_tol;         %simplices with a volume smaller than this are considered empty


vert=extreme(P,Options);        %compute/extract extreme vertices  
P.vertices = vert;              % update extreme points
%compute delaunay triangulation using qhull
if mpt_matlabrelease >= 14,
    % matlab R14 should use additional options
    adj = delaunayn(vert, {'Qt', 'Qbb', 'Qc', 'Qz'});
else
    % however matlab R13 does not support callinf delaunayn() with two
    % input arguments
    adj = delaunayn(vert);
end

nx=size(vert(adj(1,:)',:),2);   %number of states
%compute V representation of simplices
vPoly = cell(1,size(adj,1));
for i=1:size(adj,1)
    vPoly{i}=vert(adj(i,:)',:);
end 
if(Options.debug_level>0)
    check=0;
    for i=1:length(vPoly)
        for j=1:size(vPoly{i},1)
            check=check+~isinside(P,vPoly{i}(j,:)');
        end
    end
    if(check>0)
        error(['Error in triangulation:  ' num2str(check) ' vertices is not in feasible set'])
    end
end


hhPoly=polytope;
if(Options.noHPoly==0) | nargout > 2,
    ctr=0;
    %compute H representation of simplices
    index=[];
    for i=1:size(adj,1)
        tmpPoly=hull(vPoly{i},Options);
        if(isfulldim(tmpPoly))
            [H,K]=double(tmpPoly);
            addpoly=1;
            for j=1:size(vPoly{i},1)
                tmp=find(abs(H*vPoly{i}(j,:)'-K)<=Options.abs_tol); %find intersections with hyperplanes
                if length(tmp)~=nx
                    addpoly=0;
                    %check if each point is really a vertex
                end
            end
            [xc,rad]=chebyball(tmpPoly); %only consider simplices with "reasonable" size
            if(addpoly & rad>=minRad & sub_computeSimplexVol(vPoly{i})>minVol)
                index=[index i];
                ctr=ctr+1;
                hhPoly=[hhPoly tmpPoly]; %build array of simplices
                vvPoly{ctr}=vPoly{i};
            end
        end
    end
    adj=adj(index,:);%remove tiny combinations
end 


if(nargout>3)
    vol=0;
    for i=1:length(index)
        Vnew=[];
        vol=vol+sub_computeSimplexVol(vvPoly{i});
    end
end


return


function [vol]=sub_computeSimplexVol(vert)   
    %compute volume of polytope P
    nx=size(vert,2);
    vol=0;
    Vnew=[vert(1:end-1,:)-repmat(vert(end,:),nx,1)]; %shift simplex to origin  
    D = det(Vnew*Vnew);
    if D > 0,
        vol=D^0.5/factorial(nx);
    end
    
    return
 
