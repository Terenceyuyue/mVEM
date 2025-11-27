function g = decompGeo(P1,P2)
%decompGeo Generates the decomposed geometry matrix according to starting points 
% and ending points of boundary edges
% The line segments must be in the conventional directions.
%

% Copyright (C)  Terence Yu

N = size(P1,1);
g = zeros(7,N);
g(1,:) = 2;  % 2 for line 
g(2,:) = P1(:,1); % x: starting points
g(3,:) = P2(:,1); % x: ending points
g(4,:) = P1(:,2); % y: starting points
g(5,:) = P2(:,2); % y: ending points
g(6,:) = 1; % label of subdomain on the left
g(7,:) = 0; % label of subdomain on the right