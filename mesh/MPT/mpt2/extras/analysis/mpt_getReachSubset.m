function [Pret,keptrows,feasible] = mpt_getReachSubset(P,Pfin,A,B,F,G,horizon,Options)
%MPT_GETREACHSUBSET Computes a subset of P which enters Pfin
%
% [Pret,keptrows,feasible] = mpt_getReachSubset(P,Pfin,A,B,F,G,horizon,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% this function returns a subset Pret={x | H x <= K} of a polyhedron defined by P
% from which all points enter a final polytope Pfin in N steps. This function is 
% almost identical to "domain". However, here we assume the feedback 
% to be time-varying while "domain" assumes the feedback to be time-invariant.
%
% ---------------------------------------------------------------------------
% INPUT
% --------------------------------------------------------------------------- 
% P                 - polyhedral partition given by the polytopic array
% Pfin              - target polytope which should be reached in "horizon" steps
% Fi,Gi             - optimal feedback law given by u=Fi{i}*x+Gi{i}
% A, B              - system dynamics are given by x(k+1)=Ax(k)+Bu(k)
% horizon           - number of steps in which the target poly should be reached
% Options.abs_tol   - if abs(x)<tolerance x is considered to be zero
% Options.lpsolver  - LP solver to be used (see help mpt_solveLP for details)
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
%   OUTPUT:
% ---------------------------------------------------------------------------
% Pret       - subset Pret.H*x<=Pret.K of the input poly P.H*x<=P.K whose
%              points enter (P) in "horizon" steps
% keptrows   - is a vector containing the origin of each facet of Pret.H*x<=Pret.K
%               1 means that the facet originated from active constraints in P
%               0 means that the facet originated from reachability restrictions Pfin
% feasible   - is 1 if the resulting region is not empty
%
% see also DOMAIN, RANGE

% Copyright is with the following author(s):
%
% (C) 2003 Michal Kvasnica, Automatic Control Laboratory, ETH Zurich,
%          kvasnica@control.ee.ethz.ch
% (C) 2002 Pascal Grieder, Automatic Control Laboratory, ETH Zurich,
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

error(nargchk(7,8,nargin));

global mptOptions;

if ~isstruct(mptOptions),
    mpt_error;
end

if(nargin<8)
    Options=[];
end
if ~isfield(Options,'lpsolver'),
    Options.lpsolver=mptOptions.lpsolver;
end
if ~isfield(Options,'abs_tol'),
    Options.abs_tol=mptOptions.abs_tol;
end
if(horizon > length(G))
    error('mpt_getReachSubset: Cannot compute more steps than feedback law has entries')
end
if(~isfulldim(P))
    Pret=polytope;
    keptrows=[];
    feasible=0;
    return
end

%compute poly
lengthA =   length(A);
AnB     =   zeros(lengthA,1);
sumBF   =   zeros(lengthA,lengthA);
sumBG   =   zeros(lengthA,1);
%xFinal = [A^horizon + A^(horizon-1) BF + ... + BF]x0+ [ A^(horizon-1)BG + ... + BG]
for i=1:horizon
    j=i*size(B,2);
    AnB     = A^(i-1)*B; 
    sumBG   = sumBG+AnB*G((horizon*size(B,2)-j+1):(horizon*size(B,2)-j+size(B,2)),:);
    sumBF   = sumBF+AnB*F((horizon*size(B,2)-j+1):(horizon*size(B,2)-j+size(B,2)),:);   
end
sumABF  =   A^horizon + sumBF;

% xFinal = sumABF * x0 + sumBG
% xFinal is in the target region:     Hfin*xFinal<=Kfin
%                                   (Hfin*sumABF)*x0    <=(Kfin - Hfin*sumBG)
[Hfin, Kfin]=double(Pfin);
Kret = Kfin-Hfin*sumBG;
Hret = Hfin*sumABF;
if ~isnormal(P),
    P=normalize(P);
end

[H,K]=double(P);
Kret = [Kret;K];
Hret = [Hret;H];
reachLength=size(Kret,1)-size(K,1);
Pret=polytope(Hret,Kret,0,2);
kept = [];
if isfulldim(Pret),
    [Pret,kept] = reduce(Pret);
end
%check if the region obtained has full dimension, i.e. chebychev radius >0

[x,R]=chebyball(Pret,Options);
if (R<0 | ~isinside(Pret, x, Options)),
    feasible = 0;
    keptrows = [];
    Pret=polytope;
    return
end

%%kept = getkeptrows(Pret);

for i=1:length(kept)
    if(kept(i)<=reachLength)
        %reachability
        keptrows(i) = 0;
    else
        %constraint
        keptrows(i) = 1;
    end
end


if(~isfulldim(Pret))
    feasible = 0;
    Pret=polytope;
else
    feasible = 1;
end

return
