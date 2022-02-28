function [ErrH1,ErrL2,h] = getError3(node3,elem3,uh,info,pde,kOrder)
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
% exact solution
ue = pde.uexact;
Due = pde.Du;
% auxiliary data structure
if isfield(info,'aux')
    aux = info.aux;
else
    aux = auxgeometry3(node3,elem3);
end
elem3 = aux.elem3;
NT = size(elem3,1);

%% Compute L2 error
ErrL2 = 0;  ErrH1 = 0;
for iel = 1:NT
    % element information
    elemf = elem3{iel};
    xK = aux.centroid3(iel,1); 
    yK = aux.centroid3(iel,2); 
    zK = aux.centroid3(iel,3);
    hK = aux.diameter3(iel);
    % vector a
    Pis = Ph{iel}; a = Pis*chi{iel};
    % integration
    am = @(x,y,z) a(1)*1 + a(2)*(x-xK)./hK ...
        + a(3)*(y-yK)./hK + a(4)*(z-zK)./hK;
    agradm = @(x,y,z) a(1)*[0+0*x, 0+0*x, 0+0*x] ...
        + a(2)*[1/hK+0*x, 0+0*x, 0+0*x] ...
        + a(3)*[0+0*x, 1/hK+0*x, 0+0*x] ...
        + a(4)*[0+0*x, 0+0*x, 1/hK+0*x];
    fun = @(x,y,z) [(ue([x,y,z])-am(x,y,z)).^2, ...  % f
        (Due([x,y,z])-agradm(x,y,z)).^2]; % 1-by-4   % Df
    % elementwise error
    Int = integralPolyhedron(fun,k+1,node3,elemf);
    ErrL2 = ErrL2 + Int(1);
    ErrH1 = ErrH1 + sum(Int(2:4));
end
h = mean(aux.diameter3);
ErrL2 = sqrt(ErrL2);   ErrH1 = sqrt(ErrH1);