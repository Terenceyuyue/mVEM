function [Eadj,aadj,badj] = adj(Er,ar,br,Ef,af,bf,C,D,b,testdegen)

global zerotol
global verbose

if(nargin < 10)
  testdegen = 1;
end;

d = size(C,2);
k = size(D,2);
N = size(C,1);

Cer = C(Er,:);
Der = D(Er,:);
ber = b(Er,:);

% Make af and ar orthogonal
af = af./norm(af);
ar = ar - (af'*ar)*af;
ar = ar./norm(ar);

% Compute a point on the affine hull of the adjacent facet
Hi = [Cer Der ber];
He = [af' zeros(1,k) bf*(1-1e-2)];
c  = [-ar' zeros(1,k)];

nz = find(sum(abs([Hi;He]),1) > zerotol);
if(isempty(find(nz==size(Hi,2))))
  error('b is zero in adj');
end;

Hi = Hi(:,nz);
He = He(:,nz);
c  = c(:,nz(1:end-1));

%[star,lambda,how,val] = tom_lp(c,Hi,He);
[star,lambda,how,val] = esp_LP(c,Hi,He);
if(how ~= 1) 
    error('LP error'); 
end;

% Test for dual degeneracy and find the equality set
[prim,dual] = esp_isdegen(c,Hi,He,star,lambda);
tmp = zeros(d+k,1);
tmp(nz(1:end-1)) = star;
star = tmp;
if(dual == 1)
	if(verbose >= 2)
		fprintf('Dual degeneracy in ADJ\n');
	end;

  Eadj = Er(esp_dual_make_unique(star,Cer,Der,ber));
else
  Eadj = Er(find(abs([Cer Der]*star - ber)<zerotol));
end;

% Compute the affine hull of the adjacent facet
[aadj,badj] = esp_projaff(C(Eadj,:),D(Eadj,:),b(Eadj),1);

