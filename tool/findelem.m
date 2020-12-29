function findelem(node,elem,range)
%Findelem highlights some elements
%
% Copyright (C) Terence Yu

hold on

NT = size(elem,1);
if ~iscell(elem) % transform to cell
    elem = mat2cell(elem,ones(NT,1),length(elem(1,:)));
end

if nargin==2
    range = (1:NT)';
end

center = zeros(length(range),2);
s = 1;
for iel = range(1):range(end)
    index = elem{iel};
    V = node(index, :); 
    center(s,:) = polycentroid(V);
    s = s+1;
end

plot(center(:,1),center(:,2),'o','LineWidth',1,'MarkerEdgeColor','k',...
    'MarkerFaceColor','y','MarkerSize',18);
text(center(:,1)-0.02,center(:,2),int2str(range),'FontSize',12,...
    'FontWeight','bold','Color','k');

hold off