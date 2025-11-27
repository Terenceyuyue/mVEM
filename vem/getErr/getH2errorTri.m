function [ErrH2,h] = getH2errorTri(node,elem,uh,info,pde,kOrder)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s 
% kOrder: VE space containing k-th order polynomials

%% Get Ph and chi
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s
a = cell2mat(cellfun(@(A,b) (A*b)', Ph, chi, 'UniformOutput', false)); % vector a

%% Get auxiliary data
k = kOrder; % Pk-VEM (k<=3)
Nm = (k+1)*(k+2)/2;
% auxiliary data structure
isTri = true;
aux = auxgeometry(node,elem,isTri); elem = aux.elem;
NT = size(elem,1);

%% Triangulation
nodeTri = aux.nodeTri;
elemTri = aux.elemTri;

%% Quadrature points
[lambda,weight] = quadpts(k+3);

%% Compute H2 error
ErrH2 = zeros(NT,1);  % elementwise errors
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
        mxx = [0*x, 0*x, 0*x, 2/hK^2+0*x, 0*x, 0*x, 6*(x-xK)/hK^3, 2*(y-yK)/hK^3, 0*x, 0*x];
        mxy = [0*x, 0*x, 0*x, 0*x, 1/hK^2+0*x, 0*x, 0*x, 2*(x-xK)/hK^3, 2*(y-yK)/hK^3, 0*x];
        myy = [0*x, 0*x, 0*x, 0*x, 0*x, 2/hK^2+0*x, 0*x, 0*x, 2*(x-xK)/hK^3, 6*(y-yK)/hK^3];
        mxx = mxx(:,1:Nm);
        mxy = mxy(:,1:Nm);
        myy = myy(:,1:Nm);
        amxx = sum(repmat(aK,Ntri,1).*mxx,2);  
        amxy = sum(repmat(aK,Ntri,1).*mxy,2);
        amyy = sum(repmat(aK,Ntri,1).*myy,2);
        DDam = [amxx, amxy, amxy, amyy]; % xy = yx
        DDue = pde.DDu([x,y]); % [Dxx,Dxy,Dyx,Dyy]
        err = err + weight(p)*sum((DDue-DDam).^2,2);
    end
    ErrH2(iel) = dot(area,err);
end
ErrH2 = sum(ErrH2); % total errors
h = mean(aux.diameter);
ErrH2 = sqrt(ErrH2);