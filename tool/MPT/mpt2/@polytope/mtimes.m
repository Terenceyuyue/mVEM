function R=mtimes(P,Q,Options)
%MTIMES Polytope multiplication
%
% MTIMES Polytop multiplication
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%  Three cases could happen:
%   1. Multiplication of two polytopes P1*P2 leads to a higher dimensional
%      polytope [H1 0; 0 H2] <= [K1; K2]
%   2. Multiplication of a polytope by a matrix -- this leads an affine
%      transformation performed by range.m
%   3. Multiplication of a polytope by a scalar -- if the origin is included
%      in the polytope, the polytope will be shrinked
%
% USAGE:
%   R = P1*P2
%   R = P1*[1 0; 1 1]
%   R = P1*alfa
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P     - polytope / scalar / matrix
% Q     - polytope / scalar / matrix
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R     - resulting polytope
%

% Copyright is with the following author(s):
%
% (C) 2005 Mato Baotic, FER Zagreb
%          mato.baotic@fer.hr
% (C) 2003-2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%               kvasnica@control.ee.ethz.ch

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

if ~(isa(P,'polytope') | isa(Q, 'polytope'))
  error('MTIMES: Argument MUST be a polytope object');
end

if isa(P,'polytope') & isa(Q,'polytope')
    if length(P.Array)>0 | length(Q.Array)>0,
        error('MTIMES: arrays of polytope are not supported by this function');
    end
    [Pnc, Pnx] = size(P.H);
    [Qnc, Qnx] = size(Q.H);
    % form polytope in higher dimension
    H = [P.H zeros(Pnc,Qnx); zeros(Qnc,Pnx) Q.H];
    K = [P.K; Q.K];
    if P.minrep & Q.minrep,
        R = polytope(H,K,0,1);
    else
        R = polytope(H,K);
    end
    return
end

% P is a double, Q is a polytope
if isa(P,'double') & isa(Q,'polytope'),
    dummy = P;
    P = Q;
    Q = dummy;
end

% P is a polytope, Q is a double
if isa(Q,'double') & isa(P,'polytope'),
    if all(size(Q)==[1 1]),
        % scaling with a scalar

        if Q==1,
            % special case - scaling with 1 = no scaling [issue112]
            R = P;
            return
        elseif Q==0,
            % special case - scaling with zero = empty polytope
            R = polytope;
            return
        end

        % handle polyarrays:
        if length(P.Array)>0,
            R = polytope;
            for ii=1:length(P.Array),
                % make sure the polytope is in minimal and normalized
                % representation 
                A = P.Array{ii};
                if ~isminrep(A),
                    A = reduce(A);
                end
                if ~isnormal(A),
                    A = normalize(A);
                end
                R = [R polytope(A.H*sign(Q), A.K*abs(Q), 1, 1);];
            end
            return
        end
                
        % make sure the polytope is in minimal and normalized representation
        if ~isminrep(P),
            P = reduce(P);
        end
        if ~isnormal(P),
            P = normalize(P);
        end
        R = polytope(P.H*sign(Q), P.K*abs(Q), 1, 1);
        return
    else
        % P is a matrix, perform affine transformation
        R = range(P,Q);
        return
    end
end

error('MTIMES: unhandled case! see help polytope/mtimes for details...');
