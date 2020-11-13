function [prim,dual] = isdegen(c,Hi,He,x,l)
%
% [prim,dual] = isdegen(c,Hi,He,x,l)
%
% prim = 1 iff x is a primal degenerate solution to min c'x s.t. Ax<b
% dual = 1 iff lambda is a dual degenerate solution to 
%                 max lambda'b s.t. lambdaA=b, lambda<0
%
% Note l = [li;le] li - lagrange multipliers for inequalities
%
% From Murty, Theorems 4.14 and 4.15 (p.g. 209)


global zerotol

[D,d] = a2s(Hi);
D = -D; d = -d;
if(isempty(He))
  A = [];
  b = [];
else
  [A,b] =  a2s(He);
end;
mu = -l(1:size(Hi,1));
pi = -l(size(Hi,1)+1:end);

I = find(abs(D*x-d) < zerotol);
J = find(mu > zerotol);
i = abs(mu) < zerotol;
j = zeros(length(mu),1);
j(I) = 1;
L = find(i+j == 2);
%L = intersect(I,find(abs(mu) < zerotol));

nI = length(I);
nJ = length(J);
nL = length(L);

DI = D(I,:);
DJ = D(J,:);
DL = D(L,:);

dual = 0;
%if(rank([A;DI]) < max(size(DI,2),size(A,2)))
if(rank([A;DI]) < min(size(DI)))
    dual = 1;
else
    if(~isempty(L))
      Ae = [DJ;A];
      be = zeros(nJ+size(A,1),1);
      Ai = -DL;
      bi = zeros(nL,1);
%      [z,l,how,val] = tom_lp(-sum(DL,1),[Ai bi],[Ae be],1);
      [z,l,how,val] = esp_LP(-sum(DL,1),[Ai bi],[Ae be]);
      if(how ~= 1) 
        error('LP error'); 
      end;
      if(abs(val) > zerotol)
        dual = 1;
      end;
    end;
end;


prim = 0;
% $$$ if(rank([A;DJ]) < min(size(DJ)))
% $$$     prim = 1;
% $$$ else
% $$$     if(~isempty(L))
% $$$         Ae = [A' DJ' DL'];
% $$$         be = zeros(size(DJ,2),1);
% $$$         Ai = [zeros(nL,size(A,1)+nJ) -eye(nL)];
% $$$         bi = zeros(nL,1);
% $$$         AA = [Ae;-Ae;Ai];
% $$$         bb = [be;-be;bi];
% $$$         box = newpoly_box(-1000*ones(1,size(A,1)+nL+nJ),1000*ones(1,size(A,1)+nL+nJ));
% $$$         [boxA,boxb] = a2s(Hrep(box));
% $$$         AA = [AA;boxA];
% $$$         bb = [bb;boxb];
% $$$         [z,l,how,val] = cdd_lp_c(AA,bb,-[zeros(1,nJ) ones(1,nL)]);
% $$$         if(how ~= 1) 
% $$$             error('CDD error'); 
% $$$         end;
% $$$         if(abs(val) > zerotol)
% $$$             prim = 1;
% $$$         end;
% $$$     end;
% $$$ end;

