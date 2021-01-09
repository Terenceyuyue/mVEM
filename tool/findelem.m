function findelem(node,elem,varargin)
%Findelem highlights some elements
%
% Copyright (C) Terence Yu

hold on

NT = size(elem,1);
if ~iscell(elem) % transform to cell
    elem = mat2cell(elem,ones(NT,1),length(elem(1,:)));
end

range = 1:NT;
if nargin==3, range = unique(varargin{1}); end
range = range(:);

center = zeros(length(range),2);
s = 1;
for iel = range(:)' % only valid for row vector
    index = elem{iel};
    V = node(index, :); 
    center(s,:) = polycentroid(V);
    s = s+1;
end
plot(center(:,1),center(:,2),'o','LineWidth',1,'MarkerEdgeColor','k',...
    'MarkerFaceColor','y','MarkerSize',18);
text(center(:,1)-0.01,center(:,2),int2str(range),'FontSize',12,...
    'FontWeight','bold','Color','k');

hold off