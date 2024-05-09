function [uhI,nodeI,elemI] = EllipticProjection(node,elem,uh,info,kOrder)
%EllipticProjection returns the piecewise basic data structure of elliptic projection
% of uh for the conforming and nonconforming VEMs (scalar or vectorial problems).
%
% Copyright (C)  Terence Yu.

if nargin==4, kOrder = 1; end

%% nodeI, elemI
nodeI = node(horzcat(elem{:}),:);
elemLen = cellfun('length',elem);
elemI = mat2cell(1:sum(elemLen), 1, elemLen)';

%% Get Ph and chi
Ph = info.Ph;
index = info.elem2dof; % elementwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s

%% Get auxiliary data
k = kOrder; % Pk-VEM (k<=3)
Nm = (k+1)*(k+2)/2;
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter;

%% uhI
NT = size(elem,1);
uhI = cell(NT,1);
for iel = 1:NT
    % element information
    index = elem{iel};    Nv = length(index);
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);
    % scaled monomials
    m1 = @(x,y)  1+0*x;
    m2 = @(x,y) (x-xK)./hK;
    m3 = @(x,y) (y-yK)./hK;
    m4 = @(x,y) (x-xK).^2./hK^2;
    m5 = @(x,y) (x-xK).*(y-yK)./hK^2;
    m6 = @(x,y) (y-yK).^2./hK^2;
    m7 = @(x,y) (x-xK).^3./hK^3;
    m8 = @(x,y) (x-xK).^2.*(y-yK)./hK^3;
    m9 = @(x,y) (x-xK).*(y-yK).^2./hK^3;
    m10= @(x,y) (y-yK).^3./hK^3;
    m = {m1,m2,m3,m4,m5,m6,m7,m8,m9,m10};
    % coefficients of elliptic projection
    a = Ph{iel}*chi{iel};
    % elliptic projection of uh
    if length(a)<=Nm
        am = zeros(Nv,1);
        for i = 1:Nm
            am = am + a(i)*m{i}(x,y);
        end
    else % vectorial case
        am = zeros(Nv,2);
        for i = 1:Nm
            am(:,1) = am(:,1) + a(i)*m{i}(x,y);
            am(:,2) = am(:,2) + a(i+Nm)*m{i}(x,y);
        end
    end
    uhI{iel} = am;
end
uhI = vertcat(uhI{:});