function RP = PolyMesher_Reflect(P,Domain,Area)

% ---------- Compute the reflection pointset ---------------
eps = 1e-8; eta = 0.9; c = 1.5; NT = size(P,1);
Alpha = c*sqrt(Area/NT);

d = Domain('Dist',P);  
NBdrySegs = size(d,2)-1;          %Number of constituent bdry segments
n1 = (Domain('Dist',P+repmat([eps,0],NT,1))-d)/eps;
n2 = (Domain('Dist',P+repmat([0,eps],NT,1))-d)/eps;
I = abs(d(:,1:NBdrySegs))<Alpha;  %Logical index of seeds near the bdry
P1 = repmat(P(:,1),1,NBdrySegs);  %[NT x NBdrySegs] extension of P(:,1)
P2 = repmat(P(:,2),1,NBdrySegs);  %[NT x NBdrySegs] extension of P(:,2)
RP = [P1(I), P2(I)] - 2*[n1(I), n2(I)].*repmat(d(I),1,2);
d_RP = Domain('Dist',RP);
J = (d_RP(:,end)>0) & (abs(d_RP(:,end))>=eta*abs(d(I)));
RP = RP(J,:); 
RP = unique(RP,'rows');