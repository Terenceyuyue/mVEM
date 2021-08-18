function [how,maxD,maxX,minD,minX] = mpt_isPWAbigger(J1,J2,Options);
%MPT_ISPWABIGGER Test if one PWA function is bigger than a second PWA function
%
% [how,maxD,maxX,minD,minX] = fjc_isPWAbigger(J1,J2,Options);
%
% ---------------------------------------------------------------------------
% DESCRIPTION
% ---------------------------------------------------------------------------
% Test if one PWA function J1(x) is bigger than a second PWA function J2(x)
% defined over the same compact set.
%
% ---------------------------------------------------------------------------
% INPUT
% ---------------------------------------------------------------------------
%  J1, J2      - PWA functions described by Ji(x) = Ji.Bi{k} * x + Ji.Ci{k}
%                where x is in the polytopic region Ji.Pn(k)
%                Note: the total domains of both functions have to be
%                      identical!!
%  Options
%   .Vdiff_tol - tolerance for which two vertices of the polyhedral partition 
%                are considered equal (default: 1e-8)
%   .Jdiff_tol - tolerance for which two function values are considered equal 
%                (default: 1e-6)
%
% ---------------------------------------------------------------------------
% OUTPUT                                                                                                    
% ---------------------------------------------------------------------------
%  how         - =0, if J1(x)= J2(x) for all x in intersect(J1.Pn,J2.Pn)
%                =1, if J1(x)>= J2(x) for all x in intersect(J1.Pn,J2.Pn)
%                =2, if J1(x)<= J2(x) for all x in intersect(J1.Pn,J2.Pn)
%                =3, if J1(x)>= J2(x) for some x in intersect(J1.Pn,J2.Pn)
%                    and J1(x)<= J2(x) for some x in intersect(J1.Pn,J2.Pn)
%  maxD        - max_x |J1(x) - J2(x)|
%  maxX        - argmax_x |J1(x) - J2(x)|
%  minD        - min_x |J1(x) - J2(x)|
%  minX        - argmin_x |J1(x) - J2(x)|
%

% Copyright is with the following author(s):
%
%(C) 2005 Frank J. Christophersen, Automatic Control Laboratory, ETH Zurich,
%         fjc@control.ee.ethz.ch

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

global mptOptions
if ~isstruct(mptOptions)
    mpt_error;
end

if nargin < 2
  error('not enough input arguments!')
end

if nargin<3
    Options=[];
end

if isfield(Options,'Vdiff_tol')
  Vdiff_tol = Options.Vdiff_tol;
else
  Vdiff_tol = 1e-8;
end
if isfield(Options,'Jdiff_tol')
  Jdiff_tol = Options.Jdiff_tol;
else
  Jdiff_tol = 1e-6;
end


P1 = J1.Pn;
P2 = J2.Pn;


% P1 and P2 have to describe identical total domains.
DELTA1 = P1\P2;
DELTA2 = P2\P1;

if isfulldim(DELTA1) | isfulldim(DELTA2)
  error('The total domains of both functions have to be identical!')
end

nP1 = length(P1);
nP2 = length(P2);
nx  = dimension(P1(1));

if nx<=2
  Options.extreme_solver = 0; % analytical computation of the vertices
else
  Options.extreme_solver = 3; % use of CDD to compute the vertices
end


% collect unique (up to difference of Vdiff_tol) Vertices
cnt = 0;
for ii = 1:nP1
  V = [];
  V = extreme(P1(ii),Options);
  
  
  % kick out duplicate vertices if closer than Vdiff_tol (inf-norm)
  if ii>1
    for kk = 1:size(V,1)
      ind = find(sum(abs(V_total-repmat(V(kk,:),size(V_total,1),1)),2)<Vdiff_tol);

      if isempty(ind)
        cnt = cnt + 1;
        V_total = [V_total; V(kk,:)];
      end
    end
    
  else
    V_total = V;
    cnt = size(V,1);
  end

end%ii
    

for ii = 1:nP2
  V = [];
  V = extreme(P2(ii),Options);
  
  
  % kick out duplicate vertices if closer than Vdiff_tol (inf-norm)
  for kk = 1:size(V,1)
    ind = find(sum(abs(V_total-repmat(V(kk,:),size(V_total,1),1)),2)<Vdiff_tol);

    if isempty(ind)
      cnt = cnt + 1;
      V_total = [V_total; V(kk,:)];
    end
  end
end



% computation of the differences
nV    = cnt;
Jdiff = [];
Options.abs_tol = Jdiff_tol;

J1 = struct(J1);
J2 = struct(J2);
if ~isfield(J1,'Bi')
    J1.Bi = J1.details.Bi;
    J1.Ci = J1.details.Ci;
end
if ~isfield(J2,'Bi')    
    J2.Bi = J2.details.Bi;
    J2.Ci = J2.details.Ci;    
end

for ii=1:nV
  V  = [V_total(ii,:)]';
  vect = [];
  
  [isin1, ind1, closest] = isinside(P1,V,Options);
  [isin2, ind2, closest] = isinside(P2,V,Options);
  
  if isin1 & isin2
    
    hh = [];
    for i1=1:length(ind1)
      for i2=1:length(ind2)
        vect = [J1.Bi{ind1(i1)} J1.Ci{ind1(i1)}]-[J2.Bi{ind2(i2)} J2.Ci{ind2(i2)}];
        hh   = [hh; [abs(vect*[V; 1]),i1,i2]];   % collect all possible eps-values
      end
    end

    hh_sort=sortrows(hh,1);
    vect = [J1.Bi{ind1(hh_sort(1,2))} J1.Ci{ind1(hh_sort(1,2))}]-...
      [J2.Bi{ind2(hh_sort(1,3))} J2.Ci{ind2(hh_sort(1,3))}];

    Jdiff = [Jdiff; [vect*[V; 1],abs(vect*[V; 1]),V']];
  end%if
end%ii
  

Jdiff_sort=sortrows(Jdiff,2);

minD = Jdiff_sort(1,2);
maxD = Jdiff_sort(end,2);

ind_minX = find(abs(minD-Jdiff_sort(:,2))<=Jdiff_tol);
ind_maxX = find(abs(maxD-Jdiff_sort(:,2))<=Jdiff_tol);

minX = Jdiff_sort(ind_minX,3:end);
maxX = Jdiff_sort(ind_maxX,3:end);

h0 = Jdiff_sort(:,2) < Jdiff_tol;
h1 = Jdiff_sort(:,1) >= -Jdiff_tol;
h2 = Jdiff_sort(:,1) <= Jdiff_tol;
nh = length(h0);


if sum(h0) == nh  % J1 == J2
  how = 0;
elseif sum(h1) == nh % J1 >= J2
  how = 1;
elseif sum(h2) == nh % J1 <= J2
  how = 2;
else
  how = 3;
end

return
