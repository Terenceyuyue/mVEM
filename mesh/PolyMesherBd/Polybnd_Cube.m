function Pb = Polybnd_Cube(varargin)

if nargin==1
    BdBox = varargin{1};
    if length(BdBox)==6    % [x1,x2,y1,y2,z1,z2]
        x1 = BdBox(1); x2 = BdBox(2);
        y1 = BdBox(3); y2 = BdBox(4);
        z1 = BdBox(5); z2 = BdBox(6);
    end
    if length(BdBox)==2  % [x1,x2], x1=y1=z1, x2=y2=z2
        x1 = BdBox(1); x2 = BdBox(2);
        y1 = x1;       y2 = x2;
        z1 = x1;       z2 = x2;
    end
end

if nargin==2    % (x1,x2), x1=y1=z1, x2=y2=z2
    x1 = varargin{1}; x2 = varargin{2};
    y1 = x1;          y2 = x2;
    z1 = x1;          z2 = x2;
end

if nargin==6  % (x1,x2,y1,y2,z1,z2)
    x1 = varargin{1}; x2 = varargin{2};
    y1 = varargin{3}; y2 = varargin{4};
    z1 = varargin{5}; z2 = varargin{6};
end

[x,y,z] = meshgrid([x1,x2],[y1,y2],[z1,z2]);
Pb = [x(:), y(:), z(:)];