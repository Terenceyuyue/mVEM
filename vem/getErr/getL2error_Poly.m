function [ErrpL2,h] = getL2error_Poly(node,elem,ph,pde,pOrder)
% Proj(p) = ph(1)*m1 + ph(2)*m2 + ....
% Note: \int  p  d x = 0 or L_0^2 constraint is imposed.

%% ph
% the original ph = [p1; p2; p3]
% the new ph = [p1, p2, p3]
nBasis = (pOrder+2)*(pOrder+1)/2;
ph = reshape(ph, [], nBasis);

%% Get auxiliary data
% exact solution
pef = @(x,y) pde.pexact([x,y]) + 0*x;
% auxiliary data structure
aux = auxgeometry(node,elem); elem = aux.elem;  area = aux.area;
areaOmega = sum(area);
NT = size(elem,1);

%% Impose L_0^2 constraint
intp = 0;
intm = zeros(NT,nBasis);
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);    
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); 
    hK = aux.diameter(iel); 
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];

    % scaled monomials
    m1 = @(x,y) 1+0*x;
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;  
    mk1 = @(x,y) [m1(x,y), m2(x,y), m3(x,y)];

    % integral of p
    intp = intp + integralTri(pef,5,nodeT,elemT);
    
    % integral of m on each element
    intm(iel,:) = integralTri(mk1,5,nodeT,elemT);
end
pe = @(x,y) pef(x,y) - intp/areaOmega;


%% Compute L2 error
ErrpL2 = 0;
%pL2 = 0;
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);    
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); 
    hK = aux.diameter(iel); 
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    
    % scaled monomials
    m1 = @(x,y) 1+0*x;
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK; 
    m4 = @(x,y) (x-xK).^2/hK^2;            
    m5 = @(x,y) (x-xK).*(y-yK)/hK^2;       
    m6 = @(x,y) (y-yK).^2/hK^2;  
    m = {m1,m2,m3,m4,m5,m6};

    % integration 
    pf = @(x,y) 0*x;
    for j = 1:nBasis
        pf = @(x,y) pf(x,y) + ph(iel,j)*(m{j}(x,y) - intm(iel,j)/areaOmega);
    end
    errp = @(x,y) (pef(x,y)-pf(x,y)).^2;
    %p = @(x,y) pe([x,y]).^2;
    % elementwise error    
    ErrpL2 = ErrpL2 + integralTri(errp,3,nodeT,elemT);
    %pL2 = pL2 + integralTri(p,3,nodeT,elemT);
end
h = mean(aux.diameter);
ErrpL2 = sqrt(ErrpL2); %sqrt(ErrpL2)/sqrt(pL2);