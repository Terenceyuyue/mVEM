% Compute a random facet
function [E0,af,bf] = esp_shoot(C,D,b)

global zerotol
global verbose

d = size(C,2);
k = size(D,2);
N = size(C,1);

count = 0;
r = 2;
while(r ~= 1)
  count = count + 1;
  if count > 100
    error('Tried 100 random directions in shoot without finding the polyhedra!');
  end;
  
  % Choose a random direction
  q = randn(d,1); q = q./norm(q);

  % Shoot in that direction
  Hi = [C*q D b];
%  [star,lambda,how,val] = tom_lp(-[1 zeros(1,k)],Hi);
  [star,lambda,how,val] = esp_LP(-[1 zeros(1,k)],Hi);
	if(how ~= 1) error('LP error'); end;

  % Compute the affine hull of the projection
  E0 = find(abs(([C*q D]*star)-b)<zerotol)';
  [af,bf] = esp_projaff(C(E0,:),D(E0,:),b(E0));
  
  r = rank(af);
end;

% Handle dual-degeneracy
[prim,dual] = esp_isdegen(-[1 zeros(1,k)],[C*q D b],[],star,lambda);
if(dual == 1)
  if(verbose >= 2)
	fprintf('Dual degeneracy in shoot\n');
  end;
  
  E0 = esp_dual_make_unique([q*star(1);star(2:end)],C,D,b);
end;

% Compute the facet
[af,bf] = esp_projaff(C(E0,:),D(E0,:),b(E0),1);
