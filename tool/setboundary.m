function bdStruct= setboundary(node,elem,varargin)
% varargin: string for Neumann boundary

NT = size(elem,1);
if ~iscell(elem) % transform to cell
    elem = mat2cell(elem,ones(NT,1),length(elem(1,:)));
end

%% totalEdge
shiftfun = @(verts) [verts(2:end),verts(1)];  % or shiftfun = @(verts) circshift(verts,-1);
T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
v0 = horzcat(elem{:})'; % the starting points of edges
v1 = horzcat(T1{:})'; % the ending points of edges
allEdge = [v0,v1];
totalEdge = sort(allEdge,2);

%% counterclockwise bdEdge
[~,~,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
[~, i1, ~] = unique(totalEdge,'rows');
bdEdge = allEdge(i1(s==1),:);

%% set up boundary
nE = size(bdEdge,1);
% initialized as Dirichlet (true for Dirichlet, false for Neumann)
Idx = true(nE,1);
midbdEdge = (node(bdEdge(:,1),:) + node(bdEdge(:,2),:))/2;
x = midbdEdge(:,1); y = midbdEdge(:,2); %#ok<NASGU>
nvar = length(varargin); % 1 * size(varargin,2)
% nvar==0: all boundary edges are Dirichlet
if nvar==0
    Idx = true(nE,1); % true for Dirichlet
end
if nvar==1
    bdNeumann = varargin{1};
    % case 1: all Dirichlet
    if isempty(bdNeumann)
        Idx = true(nE,1);  % true for Dirichlet
        % case 2: all Neumann
    elseif strcmpi(bdNeumann,'all')
        Idx = false(nE,1);
        % case 3:partial Neumann
    else
        if mycontains(bdNeumann,'==')    % transform 'x==1' to 'abs(x-1)<1e-4'
            bdStr = bdNeumann;
            bdStr = strrep(bdStr,'==','-');
            if mycontains(bdStr,'|')
                bdStr = strrep(bdStr, '|', ')<1e-4 | abs(');
            end
            bdNeumann = ['abs(',   bdStr,   ')<1e-4'];
        end
        id = eval(bdNeumann);
        Idx(id) = false;
    end
end
if nvar>1
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
        Idx(id) = false;
    end
end
bdStruct.bdEdge = bdEdge; % all boundary edges
bdStruct.bdEdgeD = bdEdge(Idx,:); % Dirichlet boundary edges
bdStruct.bdEdgeN = bdEdge(~Idx,:); % Neumann boundary edges
bdStruct.bdNodeIdx = unique(bdEdge); % index of all boundary nodes
bdStruct.bdNodeIdxD = unique(bdEdge(Idx,:)); % index of Dirichlet boundary nodes
bdStruct.bdNodeIdxN = unique(bdEdge(~Idx,:)); % index of Neumann boundary nodes
bdEdgeIdx = find(s==1);      % index of all boundary edges
bdStruct.bdEdgeIdx = bdEdgeIdx;
bdStruct.bdEdgeIdxD = bdEdgeIdx(Idx); % index of Dirichelt boundary edges
bdStruct.bdEdgeIdxN = bdEdgeIdx(~Idx); % index of Neumann boundary edges