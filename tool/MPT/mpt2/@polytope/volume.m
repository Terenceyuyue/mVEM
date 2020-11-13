function V=volume(P,Options)
%VOLUME Calculates volume of a polytope
%
% V=volume(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Computes volume of a given polytope
%
%   V = volume(P,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                      - Polytope 
% Options.extreme_solver - Which method to use for extreme point enumeration
%                          (see help extreme)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% V     - Volume of the polytope. If input argument is a polytope array,
%         the output will be an array of volumes of corresponding polytopes.
%
% see also EXTREME
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

if ~isa(P, 'polytope')
  error('VOLUME: Argument MUST be a polytope object');
end

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if(nargin<2)
    Options=[];
end

if ~isfield(Options,'extreme_solver'),
    Options.extreme_solver=mptOptions.extreme_solver;
end
Options.noHPoly = 1;    

lenP=length(P.Array);

% check if the polytope is bounded
if ~isbounded(P),
    V = Inf;
    return
end
    

if lenP>0
    % P is a polyarray
    V=zeros(lenP,1);
    for ii=1:lenP,
        
        vert = extreme(P.Array{ii},Options);     % compute vertices of P(ii)
        try
            [nconv, V(ii)] = mpt_convhulln(vert);        % compute the volume
        catch
            % in Matlab < 6.1, convhulln does not return the volume, therefore we call 'triangulate' which computes it
            [hPoly,adj,vPoly,V(ii)]=triangulate(P.Array{ii},Options);
        end
    end
else
    % P is a polytope
    
    vertices = extreme(P,Options);        % compute vertices of P
    if(size(vertices,1)+1==size(vertices,2))
        %polytope is simplex
        nconv= 1:(size(size,1)+1);
    else
        % compute the volume        
        if mpt_matlabrelease >= 14,
            % matlab R14 should use additional options
            nconv= delaunayn(vertices, {'Qt','Qbb','Qc','Qz'});
        else
            % however matlab R13 does not support callinf delaunayn() with two
            % input arguments
            nconv= delaunayn(vertices);
        end
            
    end
    nx=size(vertices,2);
    V=0;
    for i=1:size(nconv,1)
        vert=vertices(nconv(i,:),:);
        Vnew=[vert(1:end-1,:)-repmat(vert(end,:),nx,1)]; %shift simplex to origin  
        D = det(Vnew*Vnew);
        if D > 0,
            V=V+D^0.5/factorial(nx);
        end
    end
end
