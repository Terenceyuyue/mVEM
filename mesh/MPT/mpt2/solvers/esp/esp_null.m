function Z = esp_null(A,nonempty)
%
% N = esp_null(A,nonempty)
%
% Compute the nullspace of A with a tolerance
%
% Deal with bug where null(0) = []
%
% If nonempty == 1 then the tolerance will be decreased until the
% nullspace is nonempty
%

if(nargin < 2)
  nonempty = 0;
end;

[m,n] = size(A);
[U,S,V] = svd(A,0);
if m > 1
  s = diag(S);
elseif m == 1
  s = S(1);
else
  s = 0;
end;

tol = max(m,n) * max(s) * eps;
r = sum(s > tol);
Z = V(:,r+1:n);

if(nonempty == 1 & isempty(Z))
  Z = V(:,max(n-1,1):n);
end;
