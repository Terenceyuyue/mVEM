function [R,P,Q]=plus(P,Q,Options)
%PLUS Minkowski sum of two polytopes
%
% [R,P,Q] = plus(P,Q,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% R=P+Q, or R=SUM(P,Q) computes Minkowski sum of two polytopes.
% One of the arguments can be a matrix in which case it is interpreted
% as a set of vertices of a polytope (each column is a vertex).
% This is useful for computing Minkowski sum of two polytopes that are
% not of the same dimension.
%
% NOTE: If Q ia a polyarray, Q=[Q1 Q2 ... Qn], then Q is considered
% as the space covered with polytopes Qi, and the plus.m computes P+Q according
% to the definition P+Q=\{x+y| x\in P and y\in Q\}=[(P+Q1) (P+Q2) ... (P+Qn)] 
%
% USAGE:
%   R=P+Q
%   R=plus(P,Q,Options)
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P,Q                     - Polytopes or arrays of polytopes
% Options.msolver         - 1=use projection, 2=use vertex based method
% Options.extreme_solver  - Which method to use for extreme points computation
%                           (see help extreme)
% Options.lpsolver        - Which LP solver to use (see help mpt_solveLP)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% R   - polytope (or polytope array) describing the Minkowski sum
% P,Q - since this algorithm computes extreme points, they can be stored and
%       returned back to the user (see help extreme)
%
% see also MINUS, EXTREME, MLDIVIDE
%

% Copyright is with the following author(s):
%
% (C) 2005 Mario Vasak, FER, Zagreb
%          mario.vasak@fer.hr  
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

if ~(isa(P, 'polytope') | isa(Q, 'polytope'))
  error('PLUS: Argument MUST be a polytope object');
end

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if nargin<3
   Options=[];
end

if ~isfield(Options,'extreme_solver')
    Options.extreme_solver=mptOptions.extreme_solver;
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;    % absolute tolerance
end
if ~isfield(Options,'msolver')
    % msolver=1 - use projection (use together with Options.projection / 
    % see help projection)
    % 
    % msolver=2 - minkowski sum based on extreme points and convex hull
    
    Options.msolver = 2;
end


% MINKOWSKI SUM
%==============

if isa(P,'polytope') & isa(Q,'polytope')
    lenP=length(P.Array);
    lenQ=length(Q.Array);
    if lenP>0,
        R = polytope;
        for ii=1:lenP,
            [RP,P.Array{ii},Q]=plus(P.Array{ii},Q,Options);
            if ~isfulldim(RP),
                disp('polytope/plus: extreme point enumeration failed, trying projection...');
                Opt = Options;
                Opt.msolver = 1;
                [RP,P.Array{ii},Q]=plus(P.Array{ii},Q,Opt);
                if ~isfulldim(RP)
                    Opt = Options;
                    Opt.msolver = 1;
                    [RP,P.Array{ii},Q]=plus(P.Array{ii},Q,Opt);
                end
            end
            R = [R RP];
        end
        return
    end
    if lenQ>0,
        R=mptOptions.emptypoly; %M.V. R=(P+Q1)\bigcup(P+Q2)\bigcup...\bigcup(P+Qn)
        for ii=1:lenQ,
            [P_aux,Rr,Q.Array{ii}] = plus(P,Q.Array{ii},Options);
            R=[R P_aux];
        end
        return
    end
    if ~isfulldim(P),
        R=polytope;
        return
    end
    if ~isfulldim(Q),
        R=P;
        return
    end
    if ~P.minrep
        P=reduce(P);
    end
    if ~Q.minrep
        Q=reduce(Q);
    end
    d=dimension(P);
    if dimension(Q)~=d
        error('PLUS: Polytopes MUST have the same dimension in Minkowski sum');
    end
    
    if Options.msolver==1
        % use projection
        
        A = [zeros(nconstr(P),size(Q.H,2)), P.H; Q.H, -Q.H];
        b = [P.K; Q.K];
        
        %dim = (size(P.H,2)+1):(size(P.H,2)+size(Q.H,2));
        dim = 1:size(P.H,2);
        S = polytope(A,b,0,2); 
        
        R = projection(S,dim,Options);
        return;
        
    else
        % use vertex enumeration and conex hull
        try
            [VP.V,VP.R,P]=extreme(P,Options);
            [VQ.V,VQ.R,Q]=extreme(Q,Options);
        catch
            % catch any errors in extreme point computation and if any, use
            % projection
            Options.msolver = 1;
            [R,P,Q] = plus(P,Q,Options);
            return
        end
        if size(VP.R,1)>0 | size(VQ.R,1)>0
            error('PLUS: Rays detected. Polyhedra are not allowed');
        end
        nVP=size(VP.V,1);
        nVQ=size(VQ.V,1);
        V.V=zeros(nVP*nVQ,d);
        for ii=1:nVQ
            V.V((ii-1)*nVP+1:ii*nVP,:)=VP.V+repmat(VQ.V(ii,:),nVP,1);
        end
        V=V.V;
        try
            R=hull(V,Options);
        catch
            % catch any errors in hull computation and if any, use
            % projection
            Options.msolver = 1;
            [R,P,Q] = plus(P,Q,Options);
        end
        return;
    end
end


% SECOND OPTION, one entry is a matrix
if isa(Q,'polytope')
    R=Q;
    Q=P;
else
    R=P;
end

if ~isa(Q,'double')
    error('Can handle only POLYTOPE + VECTOR(MATRIX) case');
end

% special case: Q is a vector
if size(Q,2)==1
    if dimension(R)~=length(Q),
        error('PLUS: Incompatible dimensions!');
    end
    if ~R.minrep
        R=reduce(R);
    end
    lenR = length(R.Array);
    if lenR > 0,
        % handle polyarrays
        Result = mptOptions.emptypoly;
        for ii = 1:lenR,
            Result = [Result plus(R.Array{ii}, Q, Options)];
        end
        R = Result;
    else
        % P is a single polytope
        R.vertices = [];
        R.K=R.K+R.H*Q;
        R.xCheb=R.xCheb+Q;
        R.bbox = [];
    end
else
    d=dimension(R);
    if size(Q,1)~=d
        error('Polytopes MUST have the same dimension in Minkowski sum');
    end
    
    lenR = length(R.Array);
    if lenR>0,
        Result = mptOptions.emptypoly;
        for ii=1:lenR,
            Result = [Result plus(R.Array{ii}, Q, Options)];
        end
        R = Result;
        return
    end
    
    [VR.V,VR.R,R]=extreme(R,Options);
    VQ.V=Q';
    if size(VR.R,1)>0
        error('Rays detected. Polyhedra are not allowed');
    end
    nVR=size(VR.V,1);
    nVQ=size(VQ.V,1);
    V.V=zeros(nVR*nVQ,d);
    for ii=1:nVQ
        V.V((ii-1)*nVR+1:ii*nVR,:)=VR.V+repmat(VQ.V(ii,:),nVR,1);
    end
    V=V.V;
    R=hull(V,Options);
    return;
end
