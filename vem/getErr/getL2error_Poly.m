function [ErrpL2,h] = getL2error_Poly(node,elem,ph,pde,pOrder)
% Proj(p) = ph(1)*m1 + ph(2)*m2 + ....

%% ph
% the original ph = [p1; p2; p3]
% the new ph = [p1, p2, p3]
nBasis = (pOrder+2)*(pOrder+1)/2;
ph = reshape(ph, [], nBasis);

%% Get auxiliary data
% exact solution
pe = pde.pexact;
% auxiliary data structure
aux = auxgeometry(node,elem); elem = aux.elem;
NT = size(elem,1);

%% Compute L2 error
ErrpL2 = 0;
%pL2 = 0;
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);    
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); hK = aux.diameter(iel); 
    
    % scaled monomials
    m1 = @(x,y) 1+0*x;
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK; 
    m = {m1,m2,m3};

    % integration 
    pf = @(x,y) 0*x;
    for j = 1:nBasis
        pf = @(x,y) pf(x,y) + ph(iel,j)*m{j}(x,y);
    end
    errp = @(x,y) (pe([x,y])-pf(x,y)).^2;
    %p = @(x,y) pe([x,y]).^2;
    % elementwise error
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    ErrpL2 = ErrpL2 + integralTri(errp,3,nodeT,elemT);
    %pL2 = pL2 + integralTri(p,3,nodeT,elemT);
end
h = mean(aux.diameter);
ErrpL2 = sqrt(ErrpL2); %sqrt(ErrpL2)/sqrt(pL2);