function [ErrL2,h] = getL2errorTri(node,elem,uh,info,pde,kOrder)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s
% kOrder: VE space containing k-th order polynomials

%% Get Ph, chi and a
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s
a = cell2mat(cellfun(@(A,b) (A*b)', Ph, chi, 'UniformOutput', false)); % vector a

%% Get auxiliary data
% number
k = kOrder; % Pk-VEM (k<=3)
N = size(node,1); NT = size(elem,1);
Nm = (k+1)*(k+2)/2;
% auxiliary data
isTri = true;
aux = auxgeometry(node,elem,isTri);

%% Triangulation
nodeTri = aux.nodeTri;
elemTri = aux.elemTri;

%% Quadrature points
[lambda,weight] = quadpts(k+3);

%% Compute L2 error
ErrL2 = zeros(NT,1);  % elementwise errors
for iel = 1:NT
    % --- element info ---
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2);
    hK = aux.diameter(iel);
    aK = a(iel,:);  % a1, a2, a3
    nodeT = nodeTri{iel};
    elemT = elemTri{iel};
    area = abs(simplexvolume(nodeT,elemT));

    % --- compute errors by looping over the triangles of each polygon ---
    Ntri = size(elemT,1);
    err = zeros(Ntri,1);
    for p = 1:length(weight)
        pxy = lambda(p,1)*nodeT(elemT(:,1),:) ...
            + lambda(p,2)*nodeT(elemT(:,2),:) ...
            + lambda(p,3)*nodeT(elemT(:,3),:);
        x = pxy(:,1);  y = pxy(:,2);
        switch k
            case 1
                m = [1+0*x, (x-xK)/hK, (y-yK)/hK];
            case 2
                m = [1+0*x, (x-xK)/hK, (y-yK)/hK, ...
                    (x-xK).^2/hK^2, (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2];
            case 3
                m = [1+0*x, (x-xK)./hK, (y-yK)./hK, ...
                    (x-xK).^2/hK^2, (x-xK).*(y-yK)/hK^2, (y-yK).^2/hK^2, ...
                    (x-xK).^3/hK.^3, (x-xK).^2.*(y-yK)/hK.^3, (x-xK).*(y-yK).^2/hK.^3,(y-yK).^3/hK.^3];
        end
        u = pde.uexact([x,y]);
        if size(a,2)<=Nm  % scalar case
            am = sum(repmat(aK,Ntri,1).*m,2);  
        else % vectorial case
            am = [sum(aK(:,1:Nm).*m,2), sum(aK(:,Nm+1:end).*m,2)]; % [am1, am2]
        end
        err = err + weight(p)*sum((u-am).^2,2);
    end
    ErrL2(iel) = dot(area,err);
end
ErrL2 = sum(ErrL2); % total errors
h = mean(aux.diameter);
ErrL2 = sqrt(ErrL2);