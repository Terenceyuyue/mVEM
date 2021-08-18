function P = PolyMesher_init_Pointset3D(Domain,varargin)

BdBox = Domain('BdBox');
nvar = length(varargin);

% ------------ Generate random pointset -----------
if nvar==1
    NT = varargin{1}; P = zeros(NT,3);  s = 0;
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
    return;
end


% ------------ Generate uniform pointset -----------
% if nvar==2   
    Nx = varargin{1}; Ny = varargin{2}; Nz = varargin{3};
    x = linspace(BdBox(1),BdBox(2),Nx+1)';
    y = linspace(BdBox(3),BdBox(4),Ny+1)';
    z = linspace(BdBox(5),BdBox(6),Nz+1)';
    xc = (x(1:end-1)+x(2:end))/2; 
    yc = (y(1:end-1)+y(2:end))/2;
    zc = (z(1:end-1)+z(2:end))/2;
    [X,Y,Z] = ndgrid(xc,yc,zc); P = [X(:),Y(:),Z(:)];
    d = Domain('Dist',P);  P = P(d(:,end)<0,:);
% end