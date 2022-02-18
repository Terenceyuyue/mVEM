function [ErrH1,h] = getH1error3(node3,elem3,uh,info,pde,kOrder)
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
Due = pde.Du;
% auxiliary data structure
aux = auxgeometry3(node3,elem3); elem3 = aux.elem3;
NT = size(elem3,1);

%% Compute H1 error
ErrH1 = 0;
for iel = 1:NT
    % element information
    elemf = elem3{iel};  
    hK = aux.diameter3(iel);
    % scaled monomials
    gradm1 = @(x,y,z) [0+0*x, 0+0*x, 0+0*x];
    gradm2 = @(x,y,z) [1/hK+0*x, 0+0*x, 0+0*x];
    gradm3 = @(x,y,z) [0+0*x, 1/hK+0*x, 0+0*x];
    gradm4 = @(x,y,z) [0+0*x, 0+0*x, 1/hK+0*x];
    gradm = {gradm1,gradm2,gradm3,gradm4};
    
    % vector a
    Pis = Ph{iel}; a = Pis*chi{iel};
    % integration
    if length(a)<=Nm  % scalar case
        agradm = @(x,y,z) [0+0*x, 0+0*x, 0+0*x];
        for i = 1:Nm
            agradm = @(x,y,z) agradm(x,y,z) + a(i)*gradm{i}(x,y,z);
        end
    else  % vectorial case
        agradm1 = @(x,y,z) [0+0*x, 0+0*x, 0+0*x]; 
        agradm2 = @(x,y,z) [0+0*x, 0+0*x, 0+0*x];
        for i = 1:Nm
            agradm1 = @(x,y,z) agradm1(x,y,z) + a(i)*gradm{i}(x,y,z);
            agradm2 = @(x,y,z) agradm2(x,y,z) + a(Nm+i)*gradm{i}(x,y,z);
        end
        agradm = @(x,y,z) [agradm1(x,y,z), agradm2(x,y,z)];
    end
    Df = @(x,y,z) (Due([x,y,z])-agradm(x,y,z)).^2;
    % elementwise error
    ErrH1 = ErrH1 + sum(integralPolyhedron(Df,k+1,node3,elemf));
end
h = mean(aux.diameter3);
ErrH1 = sqrt(ErrH1);