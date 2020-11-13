function P = PolyMesher_init_Pointset3D(Domain,NT)

BdBox = Domain('BdBox');
% ------------ Generate random pointset -----------
P = zeros(NT,3);  s = 0;
while s < NT
    p(:,1) = (BdBox(2)-BdBox(1))*rand(NT,1)+BdBox(1);
    p(:,2) = (BdBox(4)-BdBox(3))*rand(NT,1)+BdBox(3);
    p(:,3) = (BdBox(6)-BdBox(5))*rand(NT,1)+BdBox(5);
    d = Domain('Dist',p);
    I = find(d(:,end)<0);       % index of seeds inside the domain
    NumAdded = min(NT-s,length(I));  % number of seeds that can be added
    P(s+1:s+NumAdded,:) = p(I(1:NumAdded),:);
    s = s+NumAdded;
end