function R_P = PolyMesher_Reflect(P,Domain)

% ---------- Compute the reflection pointset ---------------
eps = 1e-8; eta = 0.9; c = 1.5; NT = size(P,1);
BdBox = Domain('BdBox'); 
Area = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
Alpha = c*sqrt(Area/NT);

d = Domain('Dist',P);  
NBdrySegs = size(d,2)-1;          %Number of constituent bdry segments
n1 = (Domain('Dist',P+[eps,0])-d)/eps;
n2 = (Domain('Dist',P+[0,eps])-d)/eps;
I = abs(d(:,1:NBdrySegs))<Alpha;  %Logical index of seeds near the bdry
P1 = repmat(P(:,1),1,NBdrySegs);  %[NT x NBdrySegs] extension of P(:,1)
P2 = repmat(P(:,2),1,NBdrySegs);  %[NT x NBdrySegs] extension of P(:,2)
R_P = [P1(I), P2(I)] - 2*[n1(I), n2(I)].*d(I);
d_R_P = Domain('Dist',R_P);
J = d_R_P(:,end)>0 & abs(d_R_P(:,end))>=eta*abs(d(I));
R_P = R_P(J,:); R_P = unique(R_P,'rows');