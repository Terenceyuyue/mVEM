function chie = getDof_C0VEM(node,elem,pde)

uexact = pde.uexact;

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge;

%% Compute the d.o.f.s
% chi1: evaluation at all vertices
chi1 = uexact(node);
% chi2: evaluation at all mid-edge vertices
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:); zc = (z1+z2)/2;
chi2 = uexact(zc);
% chi3: moments of \partial_n v on edges with given orientation
e = z1-z2;  % e = z2-z1
Ne = [-e(:,2),e(:,1)];  % scaled ne
Dw = pde.Du;
Dw1 = Dw(z1); Dwc = Dw(zc); Dw2 = Dw(z2);
wnD = sum(1/6*(Dw1+4*Dwc+Dw2).*Ne,2);
chi3 = wnD;
% chi
chie = [chi1; chi2; chi3]; 
chie = chie(:);