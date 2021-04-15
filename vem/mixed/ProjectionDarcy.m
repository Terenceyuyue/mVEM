function [uhI,phI,nodeI,elemI] = ProjectionDarcy(node,elem,uh,ph,info,pde)
% ProjectionDarcy returns the piecewise basic data structure of elliptic projection 
% of uh for mixed VEM of the Darcy problem.
%
% Copyright (C)  Terence Yu. 

%% nodeI, elemI
nodeI = node(horzcat(elem{:}),:);
elemLen = cellfun('length',elem);
elemI = mat2cell(1:sum(elemLen), 1, elemLen)';

%% Get Ph and chi
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s

%% Get auxiliary data
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; diameter = aux.diameter; 

%% uhI
NT = size(elem,1); K = pde.K;  
uhI = cell(NT,1);  phI = cell(NT,1);
for iel = 1:NT    
    % element information
    index = elem{iel};  Nv = length(index);   
    xK = centroid(iel,1); yK = centroid(iel,2); hK = diameter(iel);
    x = node(index,1); y = node(index,2);     
    % scaled monomials    
    m2 = @(x,y) (x-xK)./hK; 
    m3 = @(x,y) (y-yK)./hK;    
    % \hat{m}_a = K*grad(hK*m_{a+1})
    mh1 = @(x,y) [K(1,1)+0*x,K(2,1)+0*x];
    mh2 = @(x,y) [K(1,2)+0*x,K(2,2)+0*x];
    mh3 = @(x,y) [2*K(1,1)*m2(x,y),2*K(2,1)*m2(x,y)];
    mh4 = @(x,y) [K(1,1)*m3(x,y)+K(1,2)*m2(x,y),K(2,1)*m3(x,y)+K(2,2)*m2(x,y)];
    mh5 = @(x,y) [2*K(1,2)*m3(x,y),2*K(2,2)*m3(x,y)];  
    % coefficients of elliptic projection
    a = Ph{iel}*chi{iel};    
    % elliptic projection of uh
    uhI{iel} = a(1)*mh1(x,y)+a(2)*mh2(x,y)+a(3)*mh3(x,y)+a(4)*mh4(x,y)+a(5)*mh5(x,y);
    phI{iel} = ph(iel)*ones(Nv,1); % piecewise constant
end
uhI = vertcat(uhI{:});  % uh = [u1,u2]
phI = vertcat(phI{:});