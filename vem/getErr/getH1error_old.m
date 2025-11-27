function [ErrH1,h] = getH1error_old(node,elem,uh,info,pde,kOrder)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s
% kOrder: VE space containing k-th order polynomials

%% Get Ph and chi
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s

%% Get auxiliary data
k = kOrder; % Pk-VEM (k<=3)
Nm = (k+1)*(k+2)/2;
% exact solution
Due = pde.Du;
% auxiliary data structure
aux = auxgeometry(node,elem); elem = aux.elem;
NT = size(elem,1);

%% Compute H1 error
ErrH1 = 0;
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); hK = aux.diameter(iel);
    % scaled monomials
    gradm1 = @(x,y) [0+0*x, 0+0*x];
    gradm2 = @(x,y) [1/hK+0*x, 0+0*x];
    gradm3 = @(x,y) [0+0*x, 1/hK+0*x];
    gradm4 = @(x,y) [2*(x-xK)/hK^2, 0+0*x];
    gradm5 = @(x,y) [(y-yK)/hK^2,  (x-xK)/hK^2];
    gradm6 = @(x,y) [0+0*x, 2*(y-yK)/hK^2];
    gradm7 = @(x,y) [3*(x-xK).^2/hK^3, 0+0*x];
    gradm8 = @(x,y) [2*(x-xK).*(y-yK)/hK^3, (x-xK).^2/hK^3];
    gradm9 = @(x,y) [(y-yK).^2/hK^3, 2*(x-xK).*(y-yK)/hK^3];
    gradm10 = @(x,y) [0+0*x, 3*(y-yK).^2/hK^3];
    gradm = {gradm1,gradm2,gradm3,gradm4,gradm5,gradm6,gradm7,gradm8,gradm9,gradm10};
    
    % vector a
    Pis = Ph{iel}; a = Pis*chi{iel};
    % integration
    if length(a)<=Nm  % scalar case
        agradm = @(x,y) [0+0*x, 0+0*x];
        for i = 1:Nm
            agradm = @(x,y) agradm(x,y) + a(i)*gradm{i}(x,y);
        end
    else  % vectorial case
        agradm1 = @(x,y) [0+0*x, 0+0*x]; 
        agradm2 = @(x,y) [0+0*x, 0+0*x];
        for i = 1:Nm
            agradm1 = @(x,y) agradm1(x,y) + a(i)*gradm{i}(x,y);
            agradm2 = @(x,y) agradm2(x,y) + a(Nm+i)*gradm{i}(x,y);
        end
        agradm = @(x,y) [agradm1(x,y), agradm2(x,y)];
    end
    Df = @(x,y) (Due([x,y])-agradm(x,y)).^2;
    % elementwise error
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    ErrH1 = ErrH1 + sum(integralTri(Df,k+1,nodeT,elemT));
end
h = mean(aux.diameter);
ErrH1 = sqrt(ErrH1);