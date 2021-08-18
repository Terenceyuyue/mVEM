function R=minus(P,Q,Options)
%MINUS Pontryagin difference
%
% MINUS  Minkowski (Pontryagin) difference
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
%  The algorithm for efficiently computing the Minkowski difference between a
%  union of polytopes and a polytope is based on a talk by 
%  S. Rakovic and D. Mayne entitled "Constrained Control Computations"
%  It was the keynote address at the GBT Meeting, London, November, 2002.
%
% NOTE: If Q ia a polyarray, Q=[Q1 Q2 ... Qn], then Q is considered
% as the space covered with polytopes Qi, and the minus.m computes P-Q according
% to the definition 
%   P-Q=\{x\in P| x+w\in P \forall w\in Q\}=(P-Q1)&(P-Q2)&...&(P-Qn) 
%
% USAGE:
%   R=P-Q
%   R=minus(P,Q,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P,Q                     - Polytopes or arrays of polytopes
% Options.extreme_solver  - Which method to use for extreme points computation
%                           (see help extreme)
% Options.lpsolver        - Which LP solver to use (see help mpt_solveLP)
% Options.merge           - If set to true (default), simplifies the final
%                           polytope description using merge.m
% Options.merge_method    - Valid only if .merge is set to 1, if 'greedy'
%                           (default) greedy merging will be used,
%                           otherwise optimal merging is used
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R     - polytope (or polytope array) describing the Pontryagin difference
%
% see also PLUS, MLDIVIDE
%

% Copyright is with the following author(s):
%
% (C) 2005 Mario Vasak, FER, Zagreb
%          mario.vasak@fer.hr  
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2003 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
%          grieder@control.ee.ethz.ch
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

if ~isa(P, 'polytope'),
    error('MINUS: First input arguement must be a polytope object!');
end

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3,
    Options=[];
end

if ~isfield(Options,'lpsolver')
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if ~isfield(Options,'merge')
    Options.merge=1;
end
if Options.merge & ~isfield(Options,'merge_method')
    Options.merge_method='greedy';
end


% MINKOWSKI DIFFERENCE
%=====================

if isa(P, 'polytope') & isa(Q, 'double')
    % special case - second argument is a matrix (possible lower-dimensional
    % polytope given as a set of it's extreme vertices) 
    
    if ~isempty(P.Array),
        
        Phull=envelope(P,Options);          %compute envelope of set (Note: it is also possible to exchange "envelope" with "hull")
        %Phull=hull(P,Options);
        
        PM=minus(Phull,Q,Options);          %Compute minkowski difference
        if(~isfulldim(PM))
            R=polytope;
            return
        end
        E=mldivide(Phull,P,Options);        %E is union of polytopes; sets inside Phull which are not covered by PA
        QM = -Q;
        if ~isfulldim(E),
            Emin = polytope(QM);
        else
            Emin=plus(E,QM,Options);            %Minkowski addition on those sets
        end
        R=mldivide(PM,Emin,Options);        %Compute final polytopes
        if ~Options.merge       %M.V.
            return;
        end
    else    
    
        nx = dimension(P);
        nc = nconstr(P);
        if size(Q, 1)~=nx,
            error(sprintf('The matrix must have %d rows!', nx));
        end
        A = P.H;
        B = P.K;
        for ii = 1:nc,
            a = A(ii, :);
            b = B(ii);
            delta = min(-a*Q);
            B(ii) = B(ii) + delta;
        end
        R = polytope(A, B);
        return
    end
elseif isa(P,'polytope') & isa(Q,'polytope')
    [cx,cr]=chebyball(P);
    if isinf(cr) & all(cx==0),
        R=P;
        return
    end
    lenP=length(P.Array);
    if lenP>0,
        Phull=envelope(P,Options);          %compute envelope of set (Note: it is also possible to exchange "envelope" with "hull")
        %Phull=hull(P,Options);
        
        PM=minus(Phull,Q,Options);          %Compute minkowski difference
        if(~isfulldim(PM))
            R=polytope;
            return
        end
        E=mldivide(Phull,P,Options);        %E is union of polytopes; sets inside Phull which are not covered by PA
        QM=uminus(Q);                       %Flip Polytope 
        try
            Emin=plus(E,QM,Options);            %Minkowski addition on those sets
        catch
            disp('polytope/minus: extreme point enumeration failed, trying projection...');
            Emin=plus(E,QM,struct('msolver',1));
        end
        R=mldivide(PM,Emin,Options);        %Compute final polytopes
        if ~Options.merge       %M.V.
            return;
        end
    else    %M.V. if lenP>0 we don't enter this part, but we enter it only then when lenP=0
        lenQ=length(Q.Array);
        if lenQ>0, %M.V. we do R=(P-Q1)&(P-Q2)&...&(P-Qn)
            R = minus(P,Q.Array{1},Options);
            for ii=2:lenQ,
                R = R & minus(P,Q.Array{ii},Options);
            end
            return %result is a single polytope, no need to merge
        end 
        if ~isfulldim(Q),
            R=P;
            return
        end
        if ~isfulldim(P),
            R=polytope;
            return
        end
        if ~P.minrep
            P=reduce(P);
        end
        if ~Q.minrep
            Q=reduce(Q);
        end
        [ncP,dP]=size(P.H);
        [ncQ,dQ]=size(Q.H);
        if dP~=dQ
            error('MINUS: Polytopes MUST have the same dimension in Minkowski difference');
        end
        H=P.H;
        K=P.K;
        for ii=1:ncP
            [xopt,fval,lambda,exitflag,how]=mpt_solveLPi(-P.H(ii,:),Q.H,Q.K,[],[],[],Options.lpsolver);
            K(ii)=K(ii)+fval;
        end

        R=polytope(H, K, 1);
        return;
    end
else
    error('P must be a polytope object, and Q can either be a polytope object or a matrix of vertices');
end

if Options.merge & ~isempty(R.Array),
    mergeOpt.greedy=strcmpi(Options.merge_method,'greedy');
    R=merge(R,mergeOpt);
end
