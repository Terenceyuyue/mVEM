function auxT = auxstructure3(node3,elem3)
%auxgeometry3 gets geometry data of polyhedral mesh
%
% Copyright (C) Terence Yu.

NT = size(elem3,1);  
if ~iscell(elem3) % tetrahderal mesh: transform to cell
    elemTet = cell(NT,1);
    for iel = 1:NT
        tet = elem3(iel,:);
        face = [tet(:,[2 4 3]);tet(:,[1 3 4]);tet(:,[1 4 2]);tet(:,[1 2 3])];
        elemTet{iel} = mat2cell(face, ones(4,1), 3);  
    end
    elem3 = elemTet; 
end

%% totalFace 
% note: padding by -1 not by NaN since unique(A,'rows') does not 
% remove the repeated rows
max_n_vertices = zeros(NT,1);
for iel = 1:NT
    elemf = elem3{iel};
    max_n_vertices(iel) = max(cellfun('length', elemf));
end
max_n_vertices = max(max_n_vertices);
allFace = vertcat(elem3{:});  % cell
padding_func = @(vertex_ind) [vertex_ind,...
    -ones(1,max_n_vertices-length(vertex_ind))];  % function to pad the vacancies
allFacePad = cellfun(padding_func, allFace, 'UniformOutput', false);
allFacePad = vertcat(allFacePad{:}); % matrix
% totalFace 
totalFace = sort(allFacePad, 2);

%% face
[~,i1,totalJ] = unique(totalFace,'rows'); % first occurence
face = allFace(i1);
% i2(totalJ)= 1:length(totalJ); i2 = i2(:); % second occurence by overlapping

%% elem2face
elemLen = cellfun('length',elem3); % length of each elem
elem2face = mat2cell(totalJ',1,elemLen)';

auxT.node3 = node3;  auxT.elem3 = elem3;
auxT.face = face;
auxT.elem2face = elem2face;