function [ErrH1,h] = getH1errorTri(node,elem,uh,info,pde,kOrder)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s
% kOrder: VE space containing k-th order polynomials

%% Get Ph, chi and a
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s
a = cell2mat(cellfun(@(A,b) (A*b)', Ph, chi, 'UniformOutput', false)); % vector a

%% Get auxiliary data
% numer
k = kOrder; % Pk-VEM (k<=3)
NT = size(elem,1);
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
ErrH1 = zeros(NT,1);  % elementwise errors
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
                mx = [0+0*x, 1/hK+0*x, 0+0*x];
                my = [0+0*x, 0+0*x, 1/hK+0*x];
            case 2
                mx = [0+0*x, 1/hK+0*x, 0+0*x, 2*(x-xK)/hK^2, (y-yK)/hK^2, 0+0*x];
                my = [0+0*x, 0+0*x, 1/hK+0*x, 0+0*x, (x-xK)/hK^2,2*(y-yK)/hK^2];
            case 3
                mx = [0+0*x, 1/hK+0*x, 0+0*x, ...
                    2*(x-xK)/hK^2, (y-yK)/hK^2, 0+0*x, ...
                    3*(x-xK).^2/hK^3, 2*(x-xK).*(y-yK)/hK^3,(y-yK).^2/hK^3, 0+0*x];
                my = [0+0*x, 0+0*x, 1/hK+0*x, ...
                    0+0*x, (x-xK)/hK^2,2*(y-yK)/hK^2, ...
                    0+0*x, (x-xK).^2/hK^3, 2*(x-xK).*(y-yK)/hK^3, 3*(y-yK).^2/hK^3];
        end
        Du = pde.Du([x,y]);
        if size(a,2)<=Nm  % scalar case
            amx = sum(repmat(aK,Ntri,1).*mx,2);  amy = sum(repmat(aK,Ntri,1).*my,2);
        else % vectorial case
            amx = [sum(aK(:,1:Nm).*mx,2), sum(aK(:,Nm+1:end).*mx,2)]; % [am1, am2]
            amy = [sum(aK(:,1:Nm).*my,2), sum(aK(:,Nm+1:end).*my,2)];
            Du = Du(:,[1,3,2,4]); % [u1x,u1y,u2x,u2y] --> [u1x,u2x,u1y,u2y]
        end
        agradm = [amx, amy];
        err = err + weight(p)*sum((Du-agradm).^2,2);
    end
    ErrH1(iel) = dot(area,err);
end
ErrH1 = sum(ErrH1); % total errors
h = mean(aux.diameter);
ErrH1 = sqrt(ErrH1);