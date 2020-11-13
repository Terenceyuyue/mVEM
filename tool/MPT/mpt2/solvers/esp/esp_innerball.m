function [x,R] = esp_innerball(P,eq,Rbnd);
%
%  [x,R]=polyinnerball(P,eq,Rbnd);
%
%  eq   = row indicies of equalities
%  Rbnd = maximum bound on R (for unbounded polyhedra)
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++
% Given the polyhedron P={Ax<=B}, returns the center x and the radius R of the
% largest ball inside P. Note that R indicates how 'flat' is P, while  x is a
% good 'center' for P (the norm between x and the facets of P is maximized).
%
% [x,R]=POLYINNERBALL(A,B,LPSOLVER) also specifies the LP solver to use
%
% (C) 2001 by M. Baotic, Zurich, 26.01.2001
% Version 1.0
%
% 12/02/2003 : CNJ : modified for use with polytopic control toolbox
% 26/06/2003 : CNJ : modified to handle equalities
% 01/04/2004 : CNJ : modified for use with esp

if(nargin < 2) eq = []; end;
if(nargin < 3) Rbnd = inf; end;

[A,B] = a2s(P);
[m,n] = size(A);
ineq  = esp_setdiff([1:m],eq);

Ae = A(eq,:);
Be = B(eq);
Ae = [Ae;-Ae];
Be = [Be;-Be];
Ai = A(ineq,:);
Bi = B(ineq);

Anorm=sqrt(sum(Ai.*Ai,2));

AAi = [Ai Anorm];
AAe = [Ae zeros(size(Ae,1),1)];
bbi = Bi;
bbe = Be;
c   = [zeros(1,n) -1];

if(Rbnd ~= inf)
  AAi = [AAi;zeros(1,size(Ai,2)) 1];
  bbi = [bbi;Rbnd];
end;

%[xopt,l,how,val] = tom_lp(c,[AAi bbi],[AAe bbe]);
[xopt,l,how,val] = esp_LP(c,[AAi bbi],[AAe bbe]);

if(how ~= 1)
  error('LP error in esp_innerball');
  x = zeros(n,1);
  R = 0;
else
  R = xopt(n+1); % This is radius of the ball
  x = xopt(1:n);  % This is center of the ball
end;

if(R < 0)
    R = 0;
    x = zeros(n,1);
end;
