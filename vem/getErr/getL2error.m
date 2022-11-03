function [ErrL2,h] = getL2error(node,elem,uh,info,pde,kOrder)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s
% kOrder: VE space containing k-th order polynomials

%% Get Ph, chi and a
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s
a = cell2mat(cellfun(@(A,b) (A*b)', Ph, chi, 'UniformOutput', false)); % vector a

%% Get auxiliary data
% auxiliary data
aux = auxgeometry(node,elem); elem = aux.elem;
% number
k = kOrder; % Pk-VEM (k<=3)
N = size(node,1); NT = size(elem,1);
Nm = (k+1)*(k+2)/2;
elemLen = cellfun('length',elem);  
vertNum = unique(elemLen);

%% Triangulation
nodeTri = [node; aux.centroid];
% elemTri: [a_i, a_{i+1}, centroid]

%% Quadrature points
[lambda,weight] = quadpts(k+1);

%% Compute L2 error
ErrL2 = zeros(NT,1);  % elementwise errors
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
                    m = [ones(NTv,1), (x-xK)./hK, (y-yK)./hK];
                case 2
                    m = [ones(NTv,1), (x-xK)./hK, (y-yK)./hK, ...
                        [(x-xK).^2, (x-xK).*(y-yK), (y-yK).^2]./repmat(hK.^2,1,3)
                        ];
                case 3
                    m = [ones(NTv,1), (x-xK)./hK, (y-yK)./hK, ...
                        [(x-xK).^2, (x-xK).*(y-yK), (y-yK).^2]./repmat(hK.^2,1,3), ...
                        [(x-xK).^3, (x-xK).^2.*(y-yK), (x-xK).*(y-yK).^2,(y-yK).^3]./repmat(hK.^3,1,4) 
                        ];
            end
            if size(a,2)<=Nm  % scalar case
                am = sum(aNv.*m,2);
            else % vectorial case
                am = [sum(aNv(:,1:Nm).*m,2), sum(aNv(:,Nm+1:end).*m,2)]; % [am1, am2]
            end
            ue = pde.uexact([x,y]);
            err = err + weight(p)*(ue-am).^2;
        end
    end
    ErrL2(idNv) = area.*sum(err,2);
end
ErrL2 = sum(ErrL2); % total errors
h = mean(aux.diameter);
ErrL2 = sqrt(ErrL2);