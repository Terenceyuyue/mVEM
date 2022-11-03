function [ErrH1,h] = getH1error(node,elem,uh,info,pde,kOrder)
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
N = size(node,1); NT = size(elem,1);
Nm = (k+1)*(k+2)/2;
% auxiliary data
aux = auxgeometry(node,elem); elem = aux.elem;
elemLen = cellfun('length',elem); 
vertNum = unique(elemLen);

%% Triangulation
nodeTri = [node; aux.centroid];
% elemTri: [a_i, a_{i+1}, centroid]

%% Quadrature points
[lambda,weight] = quadpts(k+1);

%% Compute L2 error
ErrH1 = zeros(NT,1);  % elementwise errors
for Nv = vertNum(:)'
    % --- information of polygons with Nv vertices --- 
    idNv = find(elemLen == Nv); % find polygons with Nv vertices
    NTv = length(idNv);  % number of elements with Nv vertices
    elemNv = cell2mat(elem(idNv)); % elem
    xK = aux.centroid(idNv,1); yK = aux.centroid(idNv,2); 
    hK = aux.diameter(idNv); 
    aNv = a(idNv,:);  % a1, a2, a3
    
    % --- compute errors by looping over the triangles of each polygon ---
    v1 = 1:Nv;  v2 = [2:Nv,1];
    err = zeros(NTv,1);
    for s = 1:Nv  % loop over triangles of each element
        elemTri = [elemNv(:,v1(s)), elemNv(:,v2(s)), N+idNv]; 
        area = simplexvolume(nodeTri,elemTri);
        for p = 1:length(weight)
            pxy = lambda(p,1)*nodeTri(elemTri(:,1),:) ...
                + lambda(p,2)*nodeTri(elemTri(:,2),:) ...
                + lambda(p,3)*nodeTri(elemTri(:,3),:); 
            x = pxy(:,1);  y = pxy(:,2);
            switch k
                case 1
                    mx = [zeros(NTv,1), 1./hK, zeros(NTv,1)];
                    my = [zeros(NTv,1), zeros(NTv,1), 1./hK];
                case 2
                    mx = [zeros(NTv,1), 1./hK, zeros(NTv,1), ...
                          2*(x-xK)./hK.^2, (y-yK)./hK.^2, zeros(NTv,1)];
                    my = [zeros(NTv,1), zeros(NTv,1), 1./hK, ...
                          zeros(NTv,1), (x-xK)./hK.^2,2*(y-yK)./hK.^2];
                case 3
                    mx = [zeros(NTv,1), 1./hK, zeros(NTv,1), ...
                          2*(x-xK)./hK.^2, (y-yK)./hK.^2, zeros(NTv,1), ...
                          3*(x-xK).^2./hK.^3, 2*(x-xK).*(y-yK)./hK.^3,(y-yK).^2./hK.^3, zeros(NTv,1)];
                    my = [zeros(NTv,1), zeros(NTv,1), 1./hK, ...
                          zeros(NTv,1), (x-xK)./hK.^2,2*(y-yK)./hK.^2, ...
                          zeros(NTv,1), (x-xK).^2./hK.^3, 2*(x-xK).*(y-yK)./hK.^3, 3*(y-yK).^2./hK.^3];
            end
            Du = pde.Du([x,y]);  
            if size(a,2)<=Nm  % scalar case
                amx = sum(aNv.*mx,2);  amy = sum(aNv.*my,2);
            else % vectorial case
                amx = [sum(aNv(:,1:Nm).*mx,2), sum(aNv(:,Nm+1:end).*mx,2)]; % [am1, am2]
                amy = [sum(aNv(:,1:Nm).*my,2), sum(aNv(:,Nm+1:end).*my,2)];
                Du = Du(:,[1,3,2,4]); % [u1x,u1y,u2x,u2y] --> [u1x,u2x,u1y,u2y]
            end
            agradm = [amx, amy];
            err = err + weight(p)*sum((Du-agradm).^2,2);
        end
    end
    ErrH1(idNv) = area.*sum(err,2);
end
ErrH1 = sum(ErrH1); % total errors
h = mean(aux.diameter);
ErrH1 = sqrt(ErrH1);