function [H,K] = double(P)
%DOUBLE Function used to access internal properties of the given polytope
%
% [H,K] = double(P)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Function used to access internal properties of the POLYTOPE object
%
% USAGE:
%   HK = DOUBLE(P) returns matrix [H K] of a polyhedron P={x|H x<=K}
%   [H,K] = DOUBLE(P)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P   - Polytope
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% H,K          - Matrices describing the H-representation of the polytope P
%                  (i.e. P={x|H x<=K})
%
% see also polytope, chebyball, isnormal, isminrep
%

% Copyright is with the following author(s):
%
% (C) 2004 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
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

% if there are no output arguments [H K] matrix is returned
lenP=length(P.Array);

if nargout<=1
    % in case of one output argument, the matrix [H K] is returned
    if lenP>0,
        H = cell(1,lenP);
        for ii=1:lenP,
            H{ii} = [P.Array{ii}.H P.Array{ii}.K];
        end
    else
        if ~isfulldim(P),
            H = [];
        else
            H=[P.H P.K];
        end
    end
else
    % otherwise we return all fields of the polytope structure
    if lenP>0,
        H = cell(1,lenP);
        K = cell(1,lenP);
        for ii=1:lenP,
            H{ii} = P.Array{ii}.H;
            K{ii} = P.Array{ii}.K;
        end
    else
        H=P.H;
        K=P.K;
        if ~isfulldim(P),
            H = [];
            K = [];
        end
    end
end