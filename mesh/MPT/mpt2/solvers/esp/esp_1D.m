function [Q,E] = esp1D(P,ax)
%
% Project P to 1D
%

global zerotol
global verbose

H = H2normalize(P);
C = H(:,ax);
D = H(:,esp_setdiff([1:polydim(P)],ax));
b = H(:,end);

% Remove zero rows/columns
zerorows  = find(sum(abs([C D]),2) < zerotol);
zeroxcols = find(sum(abs(C),1) < zerotol);
zeroycols = find(sum(abs(D),1) < zerotol);

nonzerorows  = esp_setdiff([1:size(C,1)],zerorows);
nonzeroxcols = esp_setdiff([1:size(C,2)],zeroxcols);
nonzeroycols = esp_setdiff([1:size(D,2)],zeroycols);

%%%%% If row i is all zero and so is b(i) should this constraint be in the eset?

if(~isempty(zerorows) | ~isempty(zeroxcols) | ~isempty(zeroycols))
  if(verbose > 2)
    fprintf('All zero rows/cols present\n');
  end;
end;

% Drop the zero columns/rows
C = C(nonzerorows,nonzeroxcols);
D = D(nonzerorows,nonzeroycols);
b = b(nonzerorows);

d = length(ax);
k = size(D,2);
N = size(C,1);

c = [1 zeros(1,k)];
%[zmin,lmin,how,val] = tom_lp(c,[C D b]);
[zmin,lmin,how,val] = esp_LP(c,[C D b]);
if(how ~= 1) error('LP error'); end;

% Handle dual-degeneracy
[prim,dual] = esp_isdegen(c,[C D b],[],zmin,lmin);
if(dual == 1)
  if(verbose >= 2)
    fprintf('Dual degeneracy in esp1D\n');
  end;
  E0min = esp_dual_make_unique(zmin,C,D,b);
else
  E0min = find(abs([C D]*zmin-b)<zerotol)';
end;

e.E = sort(E0min);
e.a = -1;
e.b = -zmin(1);
E = e;
Q = [-1 -zmin(1)];

%[zmax,lmax,how,val] = tom_lp(-c,[C D b]);
[zmax,lmax,how,val] = esp_LP(-c,[C D b]);
if(how ~= 1) error('LP error'); end;

% Handle dual-degeneracy
[prim,dual] = esp_isdegen(-c,[C D b],[],zmax,lmax);
if(dual == 1)
  if(verbose >= 2)
    fprintf('Dual degeneracy in esp1D\n');
  end;
  E0max = esp_dual_make_unique(zmax,C,D,b);
else
  E0max = find(abs([C D]*zmax-b)<zerotol)';    
end;

e.E = sort(E0max);
e.a = 1;
e.b = zmax(1);
E = [E;e];
Q = [Q;1 zmax(1)];

% Return the zero rows/cols
d = length(nonzeroxcols) + length(zeroxcols);
for i=[1:length(E)]
  E(i).E  = nonzerorows(E(i).E);
  t = zeros(d,1);
  t(nonzeroxcols) = E(i).a;
  E(i).a = t;
  
  G(i,:) = E(i).a';
  g(i,1)   = E(i).b;
end;

Q = H2normalize([G g]);
