function [ErrH2,h] = getH2error(node,elem,uh,info,pde,kOrder)
% info.Ph:  elementwise Pis
% info.chi: elementwise numerical d.o.f.s 
% kOrder: VE space containing k-th order polynomials

%% Get Ph and chi
Ph = info.Ph;
index = info.elem2dof; % elemenwise global index
chi = cellfun(@(id) uh(id), index, 'UniformOutput', false); % elementwise numerical d.o.f.s

%% Get auxiliary data
k = kOrder; % Pk-VEM (k<=3)
Nm = (k+1)*(k+2)/2;
% exact solution
DDue = pde.DDu; % [Dxx,Dxy,Dyx,Dyy]
% auxiliary data structure
aux = auxgeometry(node,elem); elem = aux.elem;
NT = size(elem,1);

%% Compute H2 error
ErrH2 = 0;
for iel = 1:NT
    % element information
    index = elem{iel};  Nv = length(index);    
    xK = aux.centroid(iel,1); yK = aux.centroid(iel,2); hK = aux.diameter(iel); 
    % scaled monomials
    m1xx = @(x,y) 0*x; m1xy = @(x,y) 0*x; m1yy = @(x,y) 0*x;
    m2xx = @(x,y) 0*x; m2xy = @(x,y) 0*x; m2yy = @(x,y) 0*x;
    m3xx = @(x,y) 0*x; m3xy = @(x,y) 0*x; m3yy = @(x,y) 0*x;
    m4xx = @(x,y) 2/hK^2+0*x; m4xy = @(x,y) 0*x; m4yy = @(x,y) 0*x;
    m5xx = @(x,y) 0*x; m5xy = @(x,y) 1/hK^2+0*x; m5yy = @(x,y) 0*x;
    m6xx = @(x,y) 0*x; m6xy = @(x,y) 0*x; m6yy = @(x,y) 2/hK^2+0*x;
    m7xx = @(x,y) 6*(x-xK)/hK^3; m7xy = @(x,y) 0*x; m7yy = @(x,y) 0*x;
    m8xx = @(x,y) 2*(y-yK)/hK^3; m8xy = @(x,y) 2*(x-xK)/hK^3; m8yy = @(x,y) 0*x;
    m9xx = @(x,y) 0*x; m9xy = @(x,y) 2*(y-yK)/hK^3; m9yy = @(x,y) 2*(x-xK)/hK^3;
    m10xx = @(x,y) 0*x; m10xy = @(x,y) 0*x; m10yy = @(x,y) 6*(y-yK)/hK^3;
    mxx = {m1xx,m2xx,m3xx,m4xx,m5xx,m6xx,m7xx,m8xx,m9xx,m10xx};
    mxy = {m1xy,m2xy,m3xy,m4xy,m5xy,m6xy,m7xy,m8xy,m9xy,m10xy};
    myy = {m1yy,m2yy,m3yy,m4yy,m5yy,m6yy,m7yy,m8yy,m9yy,m10yy};    
    % vector a
    Pis = Ph{iel}; a = Pis*chi{iel};
    % integration     
    amxx = @(x,y) 0+0*x; amxy = @(x,y) 0+0*x; amyy = @(x,y) 0+0*x;
    for i = 1:Nm
        amxx = @(x,y) amxx(x,y) + a(i)*mxx{i}(x,y);
        amxy = @(x,y) amxy(x,y) + a(i)*mxy{i}(x,y);
        amyy = @(x,y) amyy(x,y) + a(i)*myy{i}(x,y);
    end
    DDam = @(x,y) [amxx(x,y), amxy(x,y), amxy(x,y), amyy(x,y)]; % xy = yx
    DDf = @(x,y) (DDue([x,y])-DDam(x,y)).^2; % xy = yx
    % elementwise error
    nodeT = [node(index,:);aux.centroid(iel,:)];
    elemT = [(Nv+1)*ones(Nv,1),(1:Nv)',[2:Nv,1]'];
    ErrH2 = ErrH2 + sum(integralTri(DDf,k,nodeT,elemT));
end
h = mean(aux.diameter);
ErrH2 = sqrt(ErrH2);  