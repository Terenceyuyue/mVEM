function V = con2vert(A,b)
% CON2VERT - convert a convex set of constraint inequalities into the set
%            of vertices at the intersections of those inequalities; i.e.,
%            solve the "vertex enumeration" problem. 
%
% Copyright (C) Michael Kleder, modified by Terence Yu.

flag = 0;
c = A\b;
tol = 1e-07;
if ~all(abs(A*c - b) < tol)
    [c,~,ef] = fminsearch( @(x)obj(x, {A,b}),c);
    if ef ~= 1, flag = 1; end
end
if flag ==1
    V = [];
    %error('Unable to locate a point within the interior of a feasible region.')
else
    b = b-A*c;
    D = A ./ repmat(b,[1,size(A,2)]);
    k = convhulln(D);
    G  = zeros(size(k,1),size(D,2));
    for ix = 1:size(k,1)
        F = D(k(ix,:),:);
        G(ix,:) = F\ones(size(F,1),1);
    end
    V = G + repmat(c',[size(G,1),1]);
    [~,I] = unique(num2str(V,6),'rows');
    V = V(I,:);
end

function d = obj(c,param)
    A = param{1};  b = param{2};
    d = A*c-b;
    k = (d>=-1e-15);
    d(k) = d(k)+1;
    d = max([0;d]);