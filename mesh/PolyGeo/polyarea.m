function area = polyarea(V)

% the vertices are in cyclic order

V1 = circshift(V,-1);
area = 0.5*abs(sum(V(:,1).*V1(:,2) - V1(:,1).*V(:,2)));