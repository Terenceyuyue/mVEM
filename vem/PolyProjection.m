function [phI,nodeI,elemI] = PolyProjection(node,elem,ph,pOrder)
% PolyProjecton returns the basic data structure of piecewise polynomial projection
% of p:
%       Proj(p) = ph(1)*m1 + ph(2)*m2 + ....
% pOrder = 0: constant
% pOrder = 1: linear
%
% Copyright (C)  Terence Yu.

%% nodeI, elemI
nodeI = node(horzcat(elem{:}),:);
elemLen = cellfun('length',elem);
elemI = mat2cell(1:sum(elemLen), 1, elemLen)';

%% ph
% the original ph = [p1; p2; p3]
% the new ph = [p1, p2, p3]
nBasis = (pOrder+2)*(pOrder+1)/2;
ph = reshape(ph, [], nBasis);

%% Get auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter;

%% uhI
NT = size(elem,1);
phI = cell(NT,1);
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    % scaled monomials
    m1 = @(x,y) 1+0*x;
    m2 = @(x,y) (x-xK)./hK;
    m3 = @(x,y) (y-yK)./hK;
    m = {m1,m2,m3};
    % L2 projection of ph
    pv = zeros(Nv,1);
    for j = 1:nBasis
        pv = pv + ph(iel,j)*m{j}(x,y);
    end
    phI{iel} = pv;
end
phI = vertcat(phI{:});