function eta = Strain_gradient_elasticity_C0VEM_indicator(node,elem,uh,info,pde)
% This function returns the local error indicator of Poisson equation in 2-D.
% 
% Copyright (C) Terence Yu.

NT = size(elem,1);
uexact = pde.uexact;

%% Get Ph and chi
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s

eta = zeros(NT,1);
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);  
    nodeK = node(index,:);
    % error    
    ue = uexact(nodeK); ue = ue(:,1);
    unum = chi{iel}; unum = unum(1:Nv,1);
    eta(iel) = norm(unum-ue);
end
