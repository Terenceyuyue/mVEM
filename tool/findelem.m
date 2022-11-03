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

d = size(node,2);

center = zeros(length(range),d);
s = 1;
for iel = range(:)' % only valid for row vector
    if d==2
        index = elem{iel};
        V = node(index, :);
        center(s,:) = polycentroid(V);
    end
    if d==3
        elemf = elem{iel};
        center(s,:) = polycentroid3(node,elemf);
    end
    s = s+1;
end

if d==2
    plot(center(:,1),center(:,2),'o','LineWidth',1,'MarkerEdgeColor','k',...
        'MarkerFaceColor','y','MarkerSize',15);
    text(center(:,1)-0.04,center(:,2),int2str(range),'FontSize',10,...
        'FontWeight','bold','Color','k');
end

if d==3
    plot3(center(:,1),center(:,2),center(:,3),'o','LineWidth',1,'MarkerEdgeColor','k',...
        'MarkerFaceColor','w','MarkerSize',18);
    text(center(:,1)-0.05,center(:,2),center(:,3),int2str(range),'FontSize',12,...
        'FontWeight','bold','Color','r');
end

hold off