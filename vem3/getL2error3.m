function [ErrL2,h] = getL2error3(node3,elem3,uh,info,pde,kOrder)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s
% kOrder: VE space containing k-th order polynomials
%
% Copyright (C)  Terence Yu.

%% Get Ph and chi
Ph = info.Ph;
indexDof = info.elem2dof; % elementwise global index
chi = cellfun(@(id) uh(id), indexDof, 'UniformOutput', false); % elementwise numerical d.o.f.s

%% Get auxiliary data
k = kOrder; % k = 1
Nm = 4;
% exact solution
ue = pde.uexact;
% auxiliary data structure
aux = auxgeometry3(node3,elem3); elem3 = aux.elem3;
NT = size(elem3,1);

%% Compute L2 error
ErrL2 = 0;
for iel = 1:NT
    % element information
    elemf = elem3{iel}; 
    xK = aux.centroid3(iel,1); yK = aux.centroid3(iel,2); zK = aux.centroid3(iel,3);
    hK = aux.diameter3(iel);
    % scaled monomials
    m1 = @(x,y,z)  1+0*x;
    m2 = @(x,y,z) (x-xK)./hK;
    m3 = @(x,y,z) (y-yK)./hK;
    m4 = @(x,y,z) (z-zK)./hK;
    m = {m1,m2,m3,m4};
    
    % vector a
    Pis = Ph{iel}; a = Pis*chi{iel};
    % integration
    if length(a)<=Nm    % scalar case
        am = @(x,y,z) 0+0*x;
        for i = 1:Nm
            am = @(x,y,z) am(x,y,z) + a(i)*m{i}(x,y,z);
        end
    else % vectorial case
        am1 = @(x,y,z) 0+0*x;  am2 = @(x,y,z) 0+0*x;
        for i = 1:Nm
            am1 = @(x,y,z) am1(x,y,z) + a(i)*m{i}(x,y,z);
            am2 = @(x,y,z) am2(x,y) + a(Nm+i)*m{i}(x,y,z);
        end
        am = @(x,y,z) [am1(x,y,z), am2(x,y,z)];
    end
    f = @(x,y,z) (ue([x,y,z])-am(x,y,z)).^2;
    % elementwise error
    ErrL2 = ErrL2 + sum(integralPolyhedron(f,k+1,node3,elemf));
end
h = mean(aux.diameter3);
ErrL2 = sqrt(ErrL2);