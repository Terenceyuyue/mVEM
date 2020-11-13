function [Ae,be,Ai,bi,E] = affhull(P)
%
% [Ae,be,Ai,bi,E] = aff(P)
%
% Compute the affine hull of P
%

global zerotol

%[A,b] = a2s(H2normalize(Hrep(P)));
[A,b] = a2s(H2normalize(P));
[m,n] = size(A);

%b = b*1e3;

e = 1e-3;
v = [zeros(1,n);eye(n)]';
v = v-repmat(mean(v')',1,n+1);
v = v*e;

bt   = b+min(-A*v,[],2);

Hi = [A -eye(m) bt];
Hi = [Hi;A zeros(m) b];
Hi = [Hi;zeros(m,n) -eye(m) zeros(m,1)];
[z,l,how,val] = esp_LP([zeros(1,n) ones(1,m)],Hi);
%[xopt,lambda,flag,fval]=esp_LP(f,H);
if(how ~= 1) 
	error('LP error'); 
end;

x = z(1:n);
s = z(n+1:end);

E = find(s > zerotol)';
%%H = Hrep(P);
H = P;
He = H(E,:);
Hi = H(esp_setdiff([1:m],E),:);
[Ae,be] = a2s(He);
[Ai,bi] = a2s(Hi);

if(isempty(E))
	E = [];
end;