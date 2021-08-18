function [pv,ev] = removeNodes(pv,ev,np)
% remove un-referenced nodes from set and re-index.

% mark referenced nodes
ki = zeros(size(pv,1),+1);
ki(+1:np) = +1;
ki(ev(:,1:2)) = +1;
%  keep referenced nodes
pv = pv(ki==+1,:);
% init re-indexing
ki(ki==+1) = 1:size(pv,1);
% edge re-indexing
ev(:,1:2) = ki(ev(:,1:2));
