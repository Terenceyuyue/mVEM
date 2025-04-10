function findnode(node,varargin)
%Findnode highlights nodes in certain range.
%
%    findnode(node)
%    findnode(node,range)
%    findnode(node,'noindex')
%    findnode(node,range,'noindex')
%
% Copyright (C) Terence Yu.

range = (1:size(node,1))';  shownum = true;

if nargin==2
    var = varargin{1};
    if isnumeric(var), range = var; end
    if iscell(var)
        if size(var{1},1) == 1   % row vector
            range = horzcat(var{:});
        else
            range = vertcat(var{:});
        end
    end
    if strcmpi(var,'noindex'), shownum = false; end
end
if nargin==3
    range = varargin{1}; shownum = false;
end
range = range(:);


dim = size(node,2);

hold on
dotColor = 'k.';

if dim==2
    plot(node(range,1),node(range,2),dotColor, 'MarkerSize', 5);
    if shownum
        shift = [0.015 0.015];
        text(node(range,1)+shift(1),node(range,2)+shift(2),int2str(range), ...
            'FontSize',10,'FontWeight','bold'); % show index number
    end
end

if dim==3
    plot3(node(range,1),node(range,2),node(range,3),dotColor, 'MarkerSize', 8);
    if shownum
        shift = [0.015 0.015 0.015];
        text(node(range,1)+shift(1),node(range,2)+shift(2),node(range,3)+shift(3),int2str(range), ...
            'FontSize',6,'FontWeight','bold'); % show index number
    end
end
hold off