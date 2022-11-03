function [ErrL2,h] = getL2error_old(node,elem,uh,info,pde,kOrder)
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
ue = pde.uexact;
% auxiliary data structure
aux = auxgeometry(node,elem); elem = aux.elem;
NT = size(elem,1);

%% Compute L2 error
ErrL2 = 0;
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); hK = aux.diameter(iel);
    % scaled monomials
    m1 = @(x,y)  1+0*x;
    m2 = @(x,y) (x-xK)./hK;
    m3 = @(x,y) (y-yK)./hK;
    m4 = @(x,y) (x-xK).^2./hK^2;
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2;
    m6 = @(x,y) (y-yK).^2./hK^2;
    m7 = @(x,y) (x-xK).^3./hK^3;
    m8 = @(x,y) (x-xK).^2.*(y-yK)./hK^3;
    m9 = @(x,y) (x-xK).*(y-yK).^2./hK^3;
    m10= @(x,y) (y-yK).^3./hK^3;
    m = {m1,m2,m3,m4,m5,m6,m7,m8,m9,m10};
    
    % vector a
    Pis = Ph{iel}; a = Pis*chi{iel};
    % integration
    if length(a)<=Nm    % scalar case
        am = @(x,y) 0+0*x;
        for i = 1:Nm
            am = @(x,y) am(x,y) + a(i)*m{i}(x,y);
        end
    else % vectorial case
        am1 = @(x,y) 0+0*x;  am2 = @(x,y) 0+0*x;
        for i = 1:Nm
            am1 = @(x,y) am1(x,y) + a(i)*m{i}(x,y);
            am2 = @(x,y) am2(x,y) + a(Nm+i)*m{i}(x,y);
        end
        am = @(x,y) [am1(x,y), am2(x,y)];
    end
    f = @(x,y) (ue([x,y])-am(x,y)).^2;
    % elementwise error
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    ErrL2 = ErrL2 + sum(integralTri(f,k+1,nodeT,elemT));
end
h = mean(aux.diameter);
ErrL2 = sqrt(ErrL2);