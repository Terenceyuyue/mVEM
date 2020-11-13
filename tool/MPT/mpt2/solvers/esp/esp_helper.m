function [Q,E] = esp_helper(P,ax,testdim)
%
% [pP,E] = esp_helper(P,ax)
%
% Compute the projection of the polytope P onto axes ax
%
% Inputs:
%   P          : Polytope P = [A b]
%   ax [1 2]   : Axes to project to
%
% Outputs:
%   pP         : Projection of P (polytope)
%   E          : Equality sets of each facet
%
%  2004/04/01
%     Colin Jones, Cambridge Control Laboratory, Cambridge UK
%     cnj22@cam.ac.uk

global zerotol
global verbose

P = H2normalize(P);
C = P(:,ax);
D = P(:,esp_setdiff([1:polydim(P)],ax));
b = P(:,end);
d = length(ax);
k = size(D,2);

% Test if the polytope is full dimensional and contains the origin
[x,R] = esp_innerball(P);

if((R < zerotol) & (testdim == 1))
  [Ae,be,Ai,bi,E]=esp_affhull(P);
  
  % Compute an internal point now that we know the affine hull
  [x,R] = esp_innerball(P,E);
  if(R > 1/zerotol)
    error('Polytope is unbounded\n');
  end;
  
  % Shift by x
  P  = esp_shift(P,x);
  be = be-Ae*x;
  bi = bi-Ai*x;
  
  axc = esp_setdiff([1:polydim(P)],ax);
  Ce = Ae(:,ax);
  De = Ae(:,axc);
  be = be;
  NDe = esp_null(De',1);
  
  % Test if the projection is full dimensional
  r = rank(NDe'*Ce,zerotol);
  if(r == 0)
    % Projection is full-D - just go ahead
    [Q,E] = esp_full(P,ax);
  else
    % Map P into the affine hull of its projection
    [U,S,V] = svd(NDe'*Ce);
    S       = S(1:r,1:r);
    Uhat    = U(:,1:r);
    Vhat    = V(:,1:r);
    Vtld    = V(:,r+1:end);
    
    t  = Vhat*inv(S)*Uhat'*NDe'*be;
    nP = [C*Vtld D b-C*t];
    
    [nQ,E] = espfulldim(nP,[1:d-r]);
    [A,b] = a2s(nQ);
    Q = [A*Vtld' b;NDe'*[Ce be];-NDe'*[Ce be]];
  end;
else
    % Shift by x
    P = esp_shift(P,x);
    [Q,E] = esp_fulldim(P,ax);
end;

% Shift polytope back to origional location
Q = esp_shift(Q,-x(ax));
