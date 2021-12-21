function chie = getDof_C0VEM(node,elem,pde)

uexact = pde.uexact;

%% Get auxiliary data
% auxgeometry
aux = auxgeometry(node,elem);
node = aux.node; elem = aux.elem;
centroid = aux.centroid; area = aux.area; diameter = aux.diameter;
% auxstructure
auxT = auxstructure(node,elem);
edge = auxT.edge;
% numbers
NT = size(elem,1);

%% Compute the d.o.f.s
% chi1: evaluation at all vertices
chi1 = uexact(node);
% chi2: evaluation at all mid-edge vertices
z1 = node(edge(:,1),:); z2 = node(edge(:,2),:); zc = (z1+z2)/2;
chi2 = uexact(zc);
% chi3: moments of \partial_n v on edges with given orientation
e = z1-z2;  % e = z2-z1
Ne = [-e(:,2),e(:,1)];  % scaled ne
Du = pde.Du;
Dw1 = Du(z1); Dwc = Du(zc); Dw2 = Du(z2); % Dw = [u1x,u1y,u2x,u2y]
wnD1 = sum(1/6*(Dw1(:,1:2)+4*Dwc(:,1:2)+Dw2(:,1:2)).*Ne,2);
wnD2 = sum(1/6*(Dw1(:,3:4)+4*Dwc(:,3:4)+Dw2(:,3:4)).*Ne,2);
wnD = [wnD1,wnD2];
chi3 = wnD;
% chi
chie = [chi1; chi2; chi3]; 
chie = chie(:);