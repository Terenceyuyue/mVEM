function [R,keptrows,feasible]=domain(P,A,f,Q,horizon,Options)
%DOMAIN Computes polytope that is mapped to an another polytope using affine map
%
% [R,keptrows,feasible]=domain(P,A,f,Q,horizon)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%   Compute polytope contained in "Q" that is mapped to "P" in "horizon" steps.
%   R=domain(P,A,f,Q,horizon) computes polytope R that is under an
%   affine mapping m: 
%   R={x(0) \in Q | x(horizon) \in P, x(k+1)=Ax(k)+f}
%
%   NOTE: Q is R^n as default
%
% This function is almost identical to "mpt_getReachSubset". However, here we
% assume the feedback to be time invariant while "mpt_getReachSubset" assumes
% the feedback to be time-variant.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P       - Polytope (or union thereof); 
% A, f    - Dynamics matrices; x(k+1)=Ax(k)+f
% Q       - Polytope (or union thereof); 
% horizon - length of the mapping
% Options.noReduce - if set to 1, resulting polytope will not be in minimal
%                    representation and 'keptrows' will be an empty matrix
%                    (off by default)
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% R         - Polytope (or union thereof); See definition in DESCRIPTION;
% keptrows   - is a vector containing the origin of each facet of R
%                1 means that the facet originated from active constraints in Q
%                0 means that the facet originated from reachability restrictions P
%              If P and Q are polygons "keptrows" returns the indices
%	           of the original polytope P for which transition R(i) was computed.         
% feasible  - 1 if a feasible transition was found (i.e., R is full dimensional);
%             0 otherwise;
%
% see also MPT_GETREACHSUBSET, RANGE
%

% ---------------------------------------------------------------------------
% Copyright is with the following author(s):
%
% (C) 2005 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
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
% ---------------------------------------------------------------------------

if ~isa(P, 'polytope')
  error('DOMAIN: Argument MUST be a polytope object');
end

error(nargchk(2,6,nargin));

global mptOptions

if nargin<6,
    Options=[];
end
if(nargin<5)
    horizon=1;
end
if(nargin<4)
    Q=polytope;
end
if nargin<3 | isempty(f)
    f=zeros(size(P.H,2),1);
end

if ~isfield(Options,'noReduce'),
    Options.noReduce = 0;
end
if ~isfield(Options, 'lpsolver'),
    Options.lpsolver = mptOptions.lpsolver;
end
% uncomment the following code to activate pre-processing using bounding boxes
% if ~isfield(Options, 'checkBboxes'),
%     checkBboxes = 1;
% else
%     checkBboxes = Options.checkBboxes;
% end

lenP=length(P.Array);
lenQ=length(Q.Array);
if lenP>0 | lenQ>0,
    %target set is union of polytopes (=polygon); check all transitions for the two polygons
    keptrows=[];
    R=polytope;
    feasible=0;
    for ii=1:max(1,lenP),
        for jj=1:max(1,lenQ)
            if(lenP>0 & lenQ==0)
                dR=domain(P.Array{ii},A,f,Q,horizon);
            elseif(lenP==0 & lenQ>0)
                dR=domain(P,A,f,Q.Array{jj},horizon);
            else
                dR=domain(P.Array{ii},A,f,Q.Array{jj},horizon);
            end
            R=[R dR];
            keptrows=[keptrows ii]; %store dynamics of transition
            feasible=feasible+isfulldim(dR);
        end
    end
    if(feasible>0)
        feasible=1;
    else
        feasible=0;
    end
    return
else
    %target set is single polytope
    Af = f;
    for i=2:horizon
        Af=Af + A^(i-1)*f;   %x(horizon)=A^horizon+Af
    end
    if ~isfulldim(Q),
        Q.H = [];
        Q.K = [];
    end
    Ahorizon = A^horizon;

    % uncomment the following code to activate pre-processing using bounding boxes
    %     if checkBboxes,
    %         % try to rule out certain transitions based on bounding boxes
    %         qbbox = Q.bbox;
    %         pbbox = P.bbox;
    %         if ~isempty(qbbox) & ~isempty(pbbox),
    %             % compute bounding box of range(Q, A, f)
    %             qrbbox = Ahorizon*qbbox + [Af Af];
    %             
    %             % check if bounding box of range(Q, A, f) intersects with P. if not,
    %             % domain will be empty
    %             if (all(pbbox(:,2) < qrbbox(:,1)) | all(pbbox(:,1) > qrbbox(:,2)))
    %                 % even bounding boxes of the two polytopes do not intersect, abort
    %                 % quickly
    %                 feasible = 0;
    %                 R = mptOptions.emptypoly;
    %                 keptrows = [];
    %                 return
    %             end
    %             % we cannot reach any conclusion based solely on the fact that bounding
    %             % boxes do intersect, therefore we continue...
    %         end
    %     end
    
    HH = [P.H*Ahorizon; Q.H];
    KK = [P.K-P.H*Af; Q.K];
    
    % compute center and radius of chebyshev's ball
    [xcheb,rcheb] = chebyball_f(HH, KK, Options);
    
    % if the ball has a non-zero radius, intersection of Q with affine
    % transformation of P exists
    feasible = (rcheb >= mptOptions.abs_tol);
    
    if feasible,
        % intersection exists, compute the domain polytope
        R = polytope(HH, KK, 0, 2, xcheb, rcheb);
        if Options.noReduce,
            % the polytope does not need reduction if this flag is true
            keptrows = [];
        else
            % reduce the polytope
            [R, keptrows] = reduce(R);
            tmp=find(keptrows<=nconstr(P));
            keptrows=ones(1,length(keptrows));
            keptrows(tmp)=0;
        end
    else
        % intersection does not exists, domain is empty
        R = mptOptions.emptypoly;
        keptrows = [];
    end
end
