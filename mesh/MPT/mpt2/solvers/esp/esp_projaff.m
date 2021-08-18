function [a,b] = esp_projaff(Ce,De,be,ExpectedDim)
%
% [a,b] = projaff(Ce,De,be,ExpectedDim)
%
%  Compute the set {x | exists y : Ce x + De y = be}
%
%  ExpectedDim : dimension that we expect the output to be. If
%  argument is included than this condition is checked and flagged if
%  not equal

global zerotol;

% Remove any zero columns from De
NZ = find(sum(abs(De),1) > zerotol);
De = De(:,NZ);

% If De is all zero, then we're already done
if(isempty(De))
  h = [Ce be];
else
  % Compute the nullspace of De
  [m,n] = size(De');
  [U,S,V] = svd(De',0);
  if     m  > 1, s = diag(S);
  elseif m == 1, s = S(1);
  else           s = 0;
  end;
  
  tol = max(m,n) * max(s) * eps;
  r = sum(s > tol);
  nDe = V(:,r+1:n);
  
  % Check if there are any singular values that fall between eps
  % accuracy and polytopic accuracy
  DodgyS = find((s>tol) & (s<zerotol));
  if(~isempty(DodgyS))
    fprintf('Matrix singular according to polytopic accurancy, but not eps accuracy\n');
  end;
  
  % Compute the projection
  h = nDe'*[Ce be];
end;

% Normalize and remove duplicates
h = H2normalize(h,1);
h = h.*repmat(sign(h(:,end)),1,size(h,2));
if(size(h,1) > 1)
  I = find(sum(abs(diff(h,1,1)),2) > zerotol);
  if(sum(abs(diff(h(end-1:end,:),1,1)) < zerotol)) 
    I = [I size(h,1)]; 
  end;
  h = h(I,:);
end;

% Test to make sure what we've done makes sense
if(nargin > 3)
  if(rank(h) ~= ExpectedDim)
    error('Projection of the affine hull is the wrong dimension');
  end;
  
  if(size(h,1) > ExpectedDim)
    error('Something funny going on here...');
  end;
end;

a = h(:,1:end-1);
b = h(:,end);
a = a';
