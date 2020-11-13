function [P] = normalize(P,Options)
%NORMALIZE Normalizes a given polytope
%
% P = normalize(P,Options)
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% PN = NORMALIZE(P) returns normalized representation PN of a
% polytope P. Polytope PN={x | H_ix <= K_i, i=1,...,nc} has the
% following property: H_i' * H_i=1, for all i.
%
% USAGE:
%   P=normalize(P)
%   P=normalize(P,Options)
%
% NOTE:
%   By default, the POLYTOPE constructor performs the normalization automatically, 
%   so there no need to use this function unless you switch the normalization
%   off while creating the polytope object.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
% P                - Polytope
% Options.rel_tol  - relative tolerance
% Options.abs_tol  - absolute tolerance
%
% Note: If Options is missing or some of the fields are not defined, the default
%       values from mptOptions will be used
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
% P   - polytope in normalized description
%
% see also REDUCE, POLYTOPE
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

global mptOptions;

rel_tol = mptOptions.rel_tol;
abs_tol = mptOptions.abs_tol;


if ~isempty(P.Array),
    lenP=length(P.Array);
    for ii=1:lenP,
        P.Array{ii}=normalize(P.Array{ii});
    end
    return
end

if P.normal
    disp('NORMALIZE Warning: Polytope is already in a normalized representation!');
    return
end

H = P.H;
K = P.K;
[nc,nx]=size(H);
Anorm=sqrt(sum(H .* H,2));


% consistency check
%------------------
if nc<1 | nx<1
    error('NORMALIZE: Polytope has to have nc>0 and nx>0');
end
if any(isinf(Anorm))
   error('NORMALIZE: No Inf terms are allowed in the matrix P.H');
end

% deal with "zero" rows in P.K
%-----------------------------
% use 2-norm combined with absolute and relative tolerance
ii=find(Anorm<=abs_tol | Anorm<=abs(rel_tol*K));

nii=length(ii);
if nii>0
    Anorm(ii)=1;
    H(ii,:)=repmat([1 zeros(1,nx-1)],nii,1);
    
    % decide if some constraint is always true
    jj=(K(ii)>=-abs_tol);
    K(ii(jj))=Inf;
    % or is it always false
    K(ii(~jj))=-Inf;
end

temp=1./Anorm;

P.normal = logical(1);
%P.normal = 1;

P.H = H .* temp(:,ones(1,nx));   % replaces P.H=P.H .* repmat(temp,1,nx); - call to repmat is rather slow
P.K = K .* temp;


return