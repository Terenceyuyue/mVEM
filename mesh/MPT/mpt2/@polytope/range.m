function R=range(P,A,f,Options)
%RANGE Affine transformation of a polytope
%
% [R]=range(P,A,f,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% R=range(P,A,f) computes affine transformation of a polytope
% affine mapping m: x |-> A*x + f maps polytope P to R, m: P -> R;
% i.e.,   R={x(1) | x(0) \in P, x(1)=Ax(k)+f}
%
% ---------------------------------------------------------------------------
% INPUT:
% --------------------------------------------------------------------------- 
% P                - Polytope
% A,f              - system dynamics matrices; x(k+1)=Ax+f; A needs to be 
%		             invertible.
% Options.abs_tol  - absolute tolerance
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT
% ---------------------------------------------------------------------------
% R   - polytope; See DESCRIPTION for details.
%
% see also DOMAIN
                                                                               
% ---------------------------------------------------------------------------              
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
global mptOptions;

if ~isa(P, 'polytope')
  error('RANGE: Argument MUST be a polytope object');
end
if ~isstruct(mptOptions),
    mpt_error;
end
if nargin<2
    error('RANGE: Not enough input arguments');
end
if nargin<3 | isempty(f)
    if(length(P.Array)>0)
        f=zeros(size(P.Array{1}.H,2),1);
    else
        f=zeros(size(P.H,2),1);
    end
end
if nargin<4
    Options=[];
end
if ~isfield(Options,'abs_tol')
    Options.abs_tol=mptOptions.abs_tol;  % absolute tolerance
end


lenP=length(P.Array);
if lenP>0,
    R=polytope;
    for ii=1:lenP,
        R=[R range(P.Array{ii},A,f,Options)];
    end
    return
end
        
if abs(det(A))<Options.abs_tol
    error('RANGE: For now only invertible mappings can be computed');
end
invA=inv(A);
R=polytope(P.H*invA, P.K+P.H*invA*f);
