function [Ae,be,Ai,bi,E]=affhull_checkall(P)

zerotol = globalzerotol;

[A,b] = a2s(Hrep(P));
N = size(A,1);

E = [];
for i=[1:N]
  Hi = [[A(1:i-1,:);A(i+1:end,:)] zeros(N-1,1) [b(1:i-1);b(i+1:end)]];
  Hi = [Hi;A(i,:) -1 b(i)];
  [z0,l,how,val] = tom_lp([zeros(1,size(A,2)) 1],Hi);
  if(how ~= 1) error('LP error'); end;
  if(abs(val) < zerotol)
	E = [E i];
  end;
end;
Ae = A(E,:);
be = b(E,:);
Ec = cnjsetdiff([1:N],E);
Ai = A(Ec,:);
bi = b(Ec,:);
