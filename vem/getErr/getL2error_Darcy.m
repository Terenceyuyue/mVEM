function [ErruL2,ErrpL2,h] = getL2error_Darcy(node,elem,uh,ph,info,pde)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s 

K = pde.K; % coefficient matrix

%% Get Ph and chi
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s
if length(ph)>size(elem,1)  % lifting mixed vem
    ph = reshape(ph, [], 3);
end

%% Get auxiliary data
% exact solution
ue = pde.uexact;  pe = pde.pexact;
% auxiliary data structure
aux = auxgeometry(node,elem); elem = aux.elem;
NT = size(elem,1);

%% Compute L2 error
ErruL2 = 0; ErrpL2 = 0;
uL2 = 0; pL2 = 0;
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);    
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); hK = aux.diameter(iel); 
    % ------- scaled monomials --------
    m1 = @(x,y) 1+0*x;
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;
    % K(\nabla hK mj), j = 2,...
    mK11 = @(x,y) K(1,1)+0*x;
    mK21 = @(x,y) K(1,2)+0*x;
    mK31 = @(x,y) 2*K(1,1)*m2(x,y);
    mK41 = @(x,y) K(1,1)*m3(x,y)+K(1,2)*m2(x,y);
    mK51 = @(x,y) 2*K(1,2)*m3(x,y);
    mK12 = @(x,y) K(2,1)+0*x;
    mK22 = @(x,y) K(2,2)+0*x;
    mK32 = @(x,y) 2*K(2,1)*m2(x,y);
    mK42 = @(x,y) K(2,1)*m3(x,y)+K(2,2)*m2(x,y);
    mK52 = @(x,y) 2*K(2,2)*m3(x,y);
    mK1 = {mK11, mK21, mK31, mK41, mK51};
    mK2 = {mK12, mK22, mK32, mK42, mK52};
    
    % vector a
    Pis = Ph{iel}; a = Pis*chi{iel};
    % integration     
    amK1 = @(x,y) 0+0*x;  
    amK2 = @(x,y) 0+0*x; 
    for i = 1:5
        amK1 = @(x,y) amK1(x,y) + a(i)*mK1{i}(x,y);
        amK2 = @(x,y) amK2(x,y) + a(i)*mK2{i}(x,y);
    end
    amK = @(x,y) [amK1(x,y), amK2(x,y)];
    errf = @(x,y) (ue([x,y])-amK(x,y)).^2;
    errp = @(x,y) (pe([x,y])-ph(iel)).^2;
    if size(ph,2)>1 % lifting mixed vem
        pf = @(x,y) ph(iel,1)*m1(x,y)+ph(iel,2)*m2(x,y)+ph(iel,3)*m3(x,y);
        errp = @(x,y) (pe([x,y])-pf(x,y)).^2;
    end
    u = @(x,y) ue([x,y]).^2;  
    p = @(x,y) pe([x,y]).^2;
    % elementwise error
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    ErruL2 = ErruL2 + sum(integralTri(errf,3,nodeT,elemT));
    ErrpL2 = ErrpL2 + integralTri(errp,3,nodeT,elemT);
    uL2 = uL2 + sum(integralTri(u,3,nodeT,elemT));
    pL2 = pL2 + integralTri(p,3,nodeT,elemT);
end
h = mean(aux.diameter);
ErruL2 = sqrt(ErruL2)/sqrt(uL2);  
ErrpL2 = sqrt(ErrpL2)/sqrt(pL2);