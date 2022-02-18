function bdStruct = setboundary3(node3,elem3,varargin)
% setboundary3 sets type of boundary faces and returns structure of boundary
% information.
% varargin: string for Neumann boundary
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

%% Find boundary faces
[~,i1,totalJ] = unique(totalFace,'rows');
i2(totalJ)= 1:length(totalJ); i2 = i2(:);
bdFace = allFace(i1(i1==i2));  %counterclockwise bdFace
bdFaceTri = allFacePad(i1(i1==i2),1:3); %triangles of bdFace

%% Set up boundary faces
nF = size(bdFace,1);
% initial as Dirichlet (true for Dirichlet, false for Neumann)
IdxD = true(nF,1);
nodebdFaceTri = (node3(bdFaceTri(:,1),:) + node3(bdFaceTri(:,2),:) + node3(bdFaceTri(:,3),:))/3;
x = nodebdFaceTri(:,1); y = nodebdFaceTri(:,2); z = nodebdFaceTri(:,3); %#ok<NASGU>
nvar = length(varargin); % 1 * size(varargin,2)
% note that length(varargin) = 1 for bdNeumann = [] or ''
if (nargin==2) || (~isempty(varargin{1}))
    for i = 1:nvar
        bdNeumann = varargin{i};
        if mycontains(bdNeumann,'==')    % transform 'x==1' to 'abs(x-1)<1e-4'
            bdStr = bdNeumann;
            bdStr = strrep(bdStr,'==','-');
            if mycontains(bdStr,'|')
                bdStr = strrep(bdStr, '|', ')<1e-4 | abs(');
            end
            bdNeumann = ['abs(',   bdStr,   ')<1e-4'];
        end
        id = eval(bdNeumann);
        IdxD(id) = false;
    end
end

bdFaceD = bdFace(IdxD,:);
bdFaceN = bdFace(~IdxD,:);
bdStruct.bdFace = bdFace;   % all boundary faces
bdStruct.bdFaceD = bdFaceD; % Dirichlet boundary faces
bdStruct.bdFaceN = bdFaceN; % Neumann boundary faces
bdFaceIdx = find(i1==i2);      % index of all boundary faces
bdFaceIdxD = bdFaceIdx(IdxD); % index of Dirichelt boundary faces
bdFaceIdxN = bdFaceIdx(~IdxD); % index of Neumann boundary faces
bdStruct.bdFaceIdx = bdFaceIdx;
bdStruct.bdFaceIdxD = bdFaceIdxD;
bdStruct.bdFaceIdxN = bdFaceIdxN;
bdStruct.bdNodeIdx = unique(horzcat(bdFace{:}));
bdStruct.bdNodeIdxD = unique(horzcat(bdFaceD{:}));
bdStruct.bdNodeIdxN = unique(horzcat(bdFaceN{:}));