function findedgeTr(node,elem,bdInd)
%FindedgeTr highlights edges for triangulation
% bdInd = 1; % boundary edge;
% other cases: all edges

hold on
% -------- edge matrix -------
totalEdge = sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2);
[i,j,s] = find(sparse(totalEdge(:,2),totalEdge(:,1),1));
edge = [j,i];
% bdEdge = edge(s==1,:);

% ------- range ---------
if nargin==2 || bdInd~=1  
    range = (1:size(edge,1))'; % all edges
else 
    range = find(s==1); % boundary edges
end

% ------ edge index ------
midEdge = (node(edge(range,1),:)+node(edge(range,2),:))/2;
plot(midEdge(:,1),midEdge(:,2),'s','LineWidth',1,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.6 0.5 0.8],'MarkerSize',20);
text(midEdge(:,1)-0.025,midEdge(:,2),int2str(range), ...
    'FontSize',12,'FontWeight','bold','Color','k');

