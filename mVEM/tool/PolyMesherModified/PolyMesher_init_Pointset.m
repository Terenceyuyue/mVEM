function P = PolyMesher_init_Pointset(Domain,NT)

BdBox = Domain('BdBox');
% ------------ Generate random pointset -----------
P = zeros(NT,2);  s = 0;
while s < NT
    xy(:,1) = (BdBox(2)-BdBox(1))*rand(NT,1)+BdBox(1);
    xy(:,2) = (BdBox(4)-BdBox(3))*rand(NT,1)+BdBox(3);
    d = Domain('Dist',xy); d = d(:,end);
    I = find(d<0);       % index of seeds inside the domain
    NumAdded = min(NT-s,length(I));  % number of seeds that can be added
    P(s+1:s+NumAdded,:) = xy(I(1:NumAdded),:);
    s = s+NumAdded;
end

