function findnode(node,range)
%Findnode highlights nodes in certain range.

dim = size(node,2);

hold on
dotColor = 'k.';
if nargin==1
    range = (1:size(node,1))';
end

if dim==2
    plot(node(range,1),node(range,2),dotColor, 'MarkerSize', 15);
    shift = [0.015 0.015];
    text(node(range,1)+shift(1),node(range,2)+shift(2),int2str(range), ...
        'FontSize',8,'FontWeight','bold'); % show index number
end

if dim==3
    plot3(node(range,1),node(range,2),node(range,3),dotColor, 'MarkerSize', 15);
    shift = [0.015 0.015 0.015];
    text(node(range,1)+shift(1),node(range,2)+shift(2),node(range,3)+shift(3),int2str(range), ...
        'FontSize',8,'FontWeight','bold'); % show index number
end
hold off