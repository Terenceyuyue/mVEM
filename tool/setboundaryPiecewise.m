function bdStruct = setboundaryPiecewise(node,elem,varargin)
% varargin: string for piecewise boundaries
% The following usages are equivalent:
%
%    setboundaryPiecewise(node,elem, {'x==1', 'y==1'} )
%    setboundaryPiecewise(node,elem, 'x==1', 'y==1' )
%
% Copyright (C)  Terence Yu.

NT = size(elem,1);
if ~iscell(elem) % transform to cell
    elem = mat2cell(elem,ones(NT,1),length(elem(1,:)));
end

%% counterclockwise bdEdge
shiftfun = @(verts) [verts(2:end),verts(1)];  % or shiftfun = @(verts) circshift(verts,-1);
T1 = cellfun(shiftfun, elem, 'UniformOutput', false);
v0 = horzcat(elem{:})'; % the starting points of edges
v1 = horzcat(T1{:})'; % the ending points of edges
allEdge = [v0,v1];
totalEdge = sort(allEdge,2);
[~,~,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
[~, i1, ~] = unique(totalEdge,'rows');
bdEdge = allEdge(i1(s==1),:);
bdEdgeIdx = find(s==1);      % index of all boundary edges

%% empty string for the whole boundary
nvar = length(varargin);
if nvar < 1
    bdEdgeType = {bdEdge};  % the whole boundary
    bdEdgeIdxType = {bdEdgeIdx};
    bdStruct.bdEdgeType = bdEdgeType;
    bdStruct.bdEdgeIdxType = bdEdgeIdxType;
    return;
end

%% determine the strings for piecewise boundaries
bdString = varargin;
if nvar == 1
    bdString = horzcat(bdString{:}); 
end
if ~iscell(bdString), bdString = {bdString}; end
nP = length(bdString);
bdEdgeType = cell(1,nP + 1); % the +1 corresponds to the remaining part
bdEdgeIdxType = cell(1,nP + 1);

%% set up boundary
midbdEdge = (node(bdEdge(:,1),:) + node(bdEdge(:,2),:))/2;
x = midbdEdge(:,1); y = midbdEdge(:,2); %#ok<NASGU>
nE = size(bdEdge,1);
IdxTotal = false(nE,1);
% parts for given strings
for i = 1:nP
    Idx = false(nE,1);
    bdStr1 = bdString{i};
    if mycontains(bdStr1,'==') 
        bdStr1 = strrep(bdStr1,'==','-');
        if mycontains(bdStr1,'|')
            bdStr1 = strrep(bdStr1, '|', ')<1e-4 | abs(');
        end
        bdStr = ['abs(',   bdStr1,   ')<1e-4'];
    else
        bdStr = bdStr1;
    end
    id = eval(bdStr);
    Idx(id) = true;
    IdxTotal(id) = true;
    bdEdgeType{i} = bdEdge(Idx,:);
    bdEdgeIdxType{i} = bdEdgeIdx(Idx,:);
end
% the remaining part
if sum(IdxTotal)<nE
    bdEdgeType{end} = bdEdge(~IdxTotal,:);
    bdEdgeIdxType{end} = bdEdgeIdx(~IdxTotal,:);
else
    bdEdgeType = bdEdgeType(1:nP);
    bdEdgeIdxType = bdEdgeIdxType(1:nP);
end

nbd = length(bdEdgeType);
bdNodeIdxType = cell(1,nbd);
for i = 1:nbd
    bdNodeIdxType{i} = unique(bdEdgeType{i});
end
bdStruct.bdEdgeType = bdEdgeType;
bdStruct.bdEdgeIdxType = bdEdgeIdxType;
bdStruct.bdNodeIdxType = bdNodeIdxType;
bdStruct.bdEdge = bdEdge;
bdStruct.bdEdgeIdx = bdEdgeIdx;